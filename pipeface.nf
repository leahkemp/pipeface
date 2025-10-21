nextflow.enable.dsl=2

// tag pipeface version
def pipeface_version = "0.9.4"

// create dummy NONE file for optional pipeface inputs
new File("NONE").text = "Dummy file for optional pipeface inputs. Don't delete during a pipeline run unless you want a bad time.\n"

// set defaults for optional params (set default optional input files to dummy NONE file)
// secret sauce second outdir
params.outdir2 = ""
params.tandem_repeat = "NONE"
params.annotate_override = ""
params.in_data_format_override = ""

// default values for chrX and chrY contig names. Used when sex = 'XY' and 
// haploidaware = 'yes'. Don't change unless your reference genome has different names. 
params.chrXseq = "chrX"
params.chrYseq = "chrY"

process scrape_settings {

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$filename" }, pattern: '*pipeface_settings.txt'

    input:
        tuple val(sample_id), val(family_id), val(files), val(data_type), val(regions_of_interest), val(clair3_model), val(family_position)
        val pipeface_version
        val in_data
        val in_data_format
        val ref
        val ref_index
        val tandem_repeat
        val mode
        val snp_indel_caller
        val sv_caller
        val annotate
        val calculate_depth
        val analyse_base_mods
        val check_relatedness
        val sites
        val outdir
        val outdir2
        val haploidaware
        val sex
        val parbed

    output:
        tuple val(sample_id), val(family_id), path('pipeface_settings.txt')

    script:
        // conditionally define reported SV caller
        if (sv_caller == 'both') {
            reported_sv_caller = 'cutesv,sniffles'
        }
        else if (sv_caller == 'sniffles') {
            reported_sv_caller = 'sniffles'
        }
        else if (sv_caller == 'cutesv') {
            reported_sv_caller = 'cutesv'
        }
        else if (sv_caller == 'NONE') {
            reported_sv_caller = 'NONE'
        }
        // conditionally define reported in data format
        if (in_data_format == 'ubam_fastq') {
            reported_in_data_format = 'unaligned BAM or FASTQ'
        }
        else if (in_data_format == 'aligned_bam') {
            reported_in_data_format = 'aligned BAM'
        }
        else if (in_data_format == 'snp_indel_vcf') {
            reported_in_data_format = 'SNP/indel VCF'
        }
        else if (in_data_format == 'sv_vcf') {
            reported_in_data_format = 'SV VCF'
        }
        if (in_data_format in ['ubam_fastq', 'aligned_bam'])
        """
        F="pipeface_settings.txt"
        echo "Pipeface version: $pipeface_version" >> \${F}
        echo "Sample ID: $sample_id" >> \${F}
        echo "Family ID: $family_id" >> \${F}
        echo "Family position: $family_position" >> \${F}
        echo "In data format: $reported_in_data_format" >> \${F}
        echo "Input data file/files: $files" >> \${F}
        echo "Data type: $data_type" >> \${F}
        echo "Regions of interest file: $regions_of_interest" >> \${F}
        echo "Clair3 model: $clair3_model" >> \${F}
        echo "In data csv path: $in_data" >> \${F}
        echo "Reference genome: $ref" >> \${F}
        echo "Reference genome index: $ref_index" >> \${F}
        echo "Haploid-aware: $haploidaware" >> \${F}
        echo "Sex: $sex" >> \${F}
        echo "PAR regions: $parbed" >> \${F}
        echo "Tandem repeat file: $tandem_repeat" >> \${F}
        echo "Mode: $mode" >> \${F}
        echo "SNP/indel caller: $snp_indel_caller" >> \${F}
        echo "SV caller: $reported_sv_caller" >> \${F}
        echo "Annotate: $annotate" >> \${F}
        echo "Calculate depth: $calculate_depth" >> \${F}
        echo "Analyse base modifications: $analyse_base_mods" >> \${F}
        echo "Check relatedness: $check_relatedness" >> \${F}
        echo "Somalier sites file: $sites" >> \${F}
        echo "Outdir: $outdir" >> \${F}
        """
        else if (in_data_format in ['snp_indel_vcf', 'sv_vcf'])
        """
        F="pipeface_settings.txt"
        echo "Pipeface version: $pipeface_version" >> \${F}
        echo "Sample ID: $sample_id" >> \${F}
        echo "Family ID: $family_id" >> \${F}
        echo "Family position: $family_position" >> \${F}
        echo "In data format: $reported_in_data_format" >> \${F}
        echo "Input data file/files: $files" >> \${F}
        echo "In data csv path: $in_data" >> \${F}
        echo "Annotate: $annotate" >> \${F}
        echo "Outdir: $outdir" >> \${F}
        """

    stub:
        """
        touch pipeface_settings.txt
        """

}

process merge_runs {

    input:
        tuple val(sample_id), val(family_id), val(extension), path(files)

    output:
        tuple val(sample_id), val(family_id), path('merged')

    script:
        if (extension == 'gz')
        """
        ${
            if (files instanceof List && files.size() > 1) {
                "cat ${files} > merged"
            } else {
                "ln -s ${files} merged"
            }
        }
        """
        else if (extension == 'fastq')
        """
        ${
            if (files instanceof List && files.size() > 1) {
                "cat ${files} | bgzip -@ ${task.cpus} > merged"
            } else {
                "ln -s ${files} merged"
            }
        }
        """
        else if (extension == 'bam')
        """
        ${
            if (files instanceof List && files.size() > 1) {
                "samtools merge -@ ${task.cpus} ${files} -o merged"
            } else {
                "ln -s ${files} merged"
            }
        }
        """

    stub:
        """
        touch merged
        """

}

process minimap2 {

    input:
        tuple val(sample_id), val(family_id), path(merged), val(extension), val(data_type)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(family_id), path('sorted.bam'), path('sorted.bam.bai')

    script:
        // conditionally define preset
        if (data_type == 'ont') {
            preset = 'lr:hq'
        }
        else if (data_type == 'pacbio') {
            preset = 'map-hifi'
        }
        if (extension == 'bam')
        """
        # run minimap
        samtools fastq -@ ${task.cpus} -T '*' $merged | minimap2 -R '@RG\\tID:${sample_id}\\tSM:${sample_id}' -y -Y --secondary=no --MD -a -x $preset -t ${task.cpus} $ref - | samtools sort -@ ${task.cpus} -o sorted.bam -
        # index bam
        samtools index -@ ${task.cpus} sorted.bam
        """
        else if (extension in ['gz', 'fastq'])
        """
        # run minimap
        minimap2 -R '@RG\\tID:${sample_id}\\tSM:${sample_id}' -Y --secondary=no --MD -a -x $preset -t ${task.cpus} $ref $merged | samtools sort -@ ${task.cpus} -o sorted.bam -
        # index bam
        samtools index -@ ${task.cpus} sorted.bam
        """

    stub:
        """
        touch sorted.bam
        touch sorted.bam.bai
        """

}

process mosdepth {

    def depth_software = "mosdepth"

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$depth_software.$filename" }, pattern: 'depth.txt'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(regions_of_interest)
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path('depth.txt')

    script:
        // define a string to optionally pass regions of interest bed file
        def regions_of_interest_optional = file(regions_of_interest).name != 'NONE' ? "-b $regions_of_interest" : ''
        """
        # run mosdepth
        mosdepth depth $bam $regions_of_interest_optional --no-per-base -t ${task.cpus}
        # rename file
        ln -s depth.mosdepth.summary.txt depth.txt
        """

    stub:
        """
        touch depth.txt
        """

}

process clair3 {

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$snp_indel_caller.$filename" }, pattern: 'snp_indel.g.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(data_type), val(regions_of_interest), val(clair3_model)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path('snp_indel.vcf.gz'), path('snp_indel.vcf.gz.tbi')
        tuple val(sample_id), path('snp_indel.g.vcf.gz'), path('snp_indel.g.vcf.gz.tbi')

    script:
        // define a string to optionally pass regions of interest bed file
        def regions_of_interest_optional = file(regions_of_interest).name != 'NONE' ? "--bed_fn=$regions_of_interest" : ''
        // conditionally define platform
        if (data_type == 'ont') {
            platform = 'ont'
        }
        else if (data_type == 'pacbio') {
            platform = 'hifi'
        }
        """
        # run clair3
        run_clair3.sh --bam_fn=$bam --ref_fn=$ref --output=./ --threads=${task.cpus} --platform=$platform --model_path=$clair3_model --sample_name=$sample_id --gvcf --include_all_ctgs $regions_of_interest_optional
        # rename files
        ln -s merge_output.vcf.gz snp_indel.vcf.gz
        ln -s merge_output.vcf.gz.tbi snp_indel.vcf.gz.tbi
        ln -s merge_output.gvcf.gz snp_indel.g.vcf.gz
        ln -s merge_output.gvcf.gz.tbi snp_indel.g.vcf.gz.tbi
        """

    stub:
        """
        touch snp_indel.vcf.gz
        touch snp_indel.vcf.gz.tbi
        """

}

process clair3_pre_processing {

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(regions_of_interest)
        val ref
        val ref_index
        path parbed

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path('diploid.bed'), path('haploid.bed')

    script:
        """
        # if regions were given, we use them to subset the genome index to construct the genome.bed file
        if [[ -f "${regions_of_interest}" ]]; then
            awk 'NR==FNR { len[\$1]=\$2; next } \$1 in len { print \$1 "\\t0\\t" len[\$1] }' \
                "${ref}.fai" "${regions_of_interest}" > genome.bed
        else
        # fallback: use full genome .fai
            awk '{ print \$1 "\\t0\\t" \$2 }' "${ref}.fai" > genome.bed
        fi
        genomebed="genome.bed"
        # separate out genome bed into a haploid and diploid bed each
        grep -v -E "^chrX|^chrY" \${genomebed} > autosomes.bed
        grep "^chrX" ${parbed} > xpar.bed
        grep "^chrY" ${parbed} > ypar.bed
        grep "^chrX" \${genomebed} > chrX_full.bed
        grep "^chrY" \${genomebed} > chrY_full.bed
        bedtools subtract -a chrX_full.bed -b xpar.bed > xnonpar.bed
        bedtools subtract -a chrY_full.bed -b ypar.bed > ynonpar.bed
        cat autosomes.bed xpar.bed ypar.bed | sort -k1,1V -k2,2n > diploid.bed
        cat xnonpar.bed ynonpar.bed | sort -k1,1 -k2,2n > haploid.bed
        """

    stub:
        """
        touch diploid.bed
        touch haploid.bed
        """

}

process clair3_haploid_aware {

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$snp_indel_caller.$filename" }, pattern: 'snp_indel.g.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path(diploid_bed), path(haploid_bed), val(data_type), val(clair3_model)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path('haploid_snp_indel.vcf.gz'), path('haploid_snp_indel.vcf.gz.tbi'), path('diploid_snp_indel.vcf.gz'), path('diploid_snp_indel.vcf.gz.tbi'), path('haploid_snp_indel.g.vcf.gz'), path('haploid_snp_indel.g.vcf.gz.tbi'), path('diploid_snp_indel.g.vcf.gz'), path('diploid_snp_indel.g.vcf.gz.tbi')

    script:
        // conditionally define platform
        if (data_type == 'ont') {
            platform = 'ont'
        }
        else if (data_type == 'pacbio') {
            platform = 'hifi'
        }
        """
        mkdir -p {haploid,diploid}
        # run clair3 (haploid)
        run_clair3.sh --bam_fn="${bam}" --ref_fn="${ref}" --output="haploid" --threads=${task.cpus} --platform=${platform} --model_path="${clair3_model}" --sample_name="${sample_id}" --gvcf --bed_fn="${haploid_bed}" --haploid_precise
        # run clair3 (diploid)
        run_clair3.sh --bam_fn="${bam}" --ref_fn="${ref}" --output="diploid" --threads=${task.cpus} --platform=${platform} --model_path="${clair3_model}" --sample_name="${sample_id}" --gvcf --bed_fn="${diploid_bed}"
        # rename files
        ln -s haploid/merge_output.vcf.gz haploid_snp_indel.vcf.gz
        ln -s haploid/merge_output.vcf.gz.tbi haploid_snp_indel.vcf.gz.tbi
        ln -s diploid/merge_output.vcf.gz diploid_snp_indel.vcf.gz
        ln -s diploid/merge_output.vcf.gz.tbi diploid_snp_indel.vcf.gz.tbi
        ln -s haploid/merge_output.gvcf.gz haploid_snp_indel.g.vcf.gz
        ln -s haploid/merge_output.gvcf.gz.tbi haploid_snp_indel.g.vcf.gz.tbi
        ln -s diploid/merge_output.gvcf.gz diploid_snp_indel.g.vcf.gz
        ln -s diploid/merge_output.gvcf.gz.tbi diploid_snp_indel.g.vcf.gz.tbi
        """

    stub:
        """
        touch haploid_merge_output.vcf.gz
        touch diploid_merge_output.vcf.gz
        touch haploid_snp_indel.g.vcf.gz
        touch diploid_snp_indel.g.vcf.gz
        """

}

process clair3_post_processing {

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$snp_indel_caller.$filename" }, pattern: 'snp_indel.g.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path(haploid_snp_indel_vcf), path(haploid_snp_indel_vcf_index), path(diploid_snp_indel_vcf), path(diploid_snp_indel_vcf_index), path(haploid_snp_indel_gvcf), path(haploid_snp_indel_gvcf_index), path(diploid_snp_indel_gvcf), path(diploid_snp_indel_gvcf_index)
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path('snp_indel.vcf.gz'), path('snp_indel.vcf.gz.tbi')
        tuple val(sample_id), path('snp_indel.g.vcf.gz'), path('snp_indel.g.vcf.gz.tbi')

    script:
        """
        # merge diploid and haploid vcf & gvcf files
        bcftools +fixploidy $haploid_snp_indel_vcf -- -f 2 -t GT > haploid_ploidyfixed.vcf
        bcftools +fixploidy $haploid_snp_indel_gvcf -- -f 2 -t GT > haploid_ploidyfixed.g.vcf
        bgzip -@ ${task.cpus} haploid_ploidyfixed.vcf
        bgzip -@ ${task.cpus} haploid_ploidyfixed.g.vcf
        tabix -p vcf haploid_ploidyfixed.vcf.gz
        tabix -p vcf haploid_ploidyfixed.g.vcf.gz
        bcftools concat -a -Oz -o unsorted.vcf.gz $diploid_snp_indel_vcf haploid_ploidyfixed.vcf.gz
        bcftools sort -Oz -o snp_indel.vcf.gz unsorted.vcf.gz
        tabix -p vcf snp_indel.vcf.gz
        bcftools concat -a -Oz -o unsorted.g.vcf.gz $diploid_snp_indel_gvcf haploid_ploidyfixed.g.vcf.gz
        bcftools sort -Oz -o snp_indel.g.vcf.gz unsorted.g.vcf.gz
        tabix -p vcf snp_indel.g.vcf.gz
        """

    stub:
        """
        touch snp_indel.g.vcf.gz
        touch snp_indel.g.vcf.gz.tbi
        """

}

process deepvariant_dry_run {

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(data_type)
        val ref
        val ref_index
        val sex
        val haploidaware
        val parbed

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), env(make_examples_args), env(call_variants_args)

    script:
        // conditionally define model type
        if (data_type == 'ont') {
            model = 'ONT_R104'
        }
        else if (data_type == 'pacbio') {
            model = 'PACBIO'
        }
        // conditionally define haploid contigs and par regions
        if (haploidaware == 'yes') {
            haploidparameter = (sex == "XX") ? "" : "--haploid_contigs ${params.chrXseq},${params.chrYseq}"
            parbedparameter = (sex == "XX") ? "" : "--par_regions_bed $parbed"
        }
        else {
            haploidparameter = ""
            parbedparameter = ""
        }
        """
        # do a dry-run of deepvariant
        run_deepvariant --reads=$bam --ref=$ref --sample_name=$sample_id --output_vcf=snp_indel.vcf.gz --output_gvcf=snp_indel.g.vcf.gz --model_type=$model --dry_run=true \
        ${haploidparameter} ${parbedparameter} > commands.txt
        # extract arguments for make_examples and call_variants stages
        make_examples_args=\$(grep "/opt/deepvariant/bin/make_examples" commands.txt | awk -F'/opt/deepvariant/bin/make_examples' '{print \$2}' | sed 's/--mode calling//g' | sed 's/--ref "[^"]*"//g' | sed 's/--reads "[^"]*"//g' | sed 's/--sample_name "[^"]*"//g' | sed 's/--examples "[^"]*"//g' | sed 's/--gvcf "[^"]*"//g')
        call_variants_args=\$(grep "/opt/deepvariant/bin/call_variants" commands.txt | awk -F'/opt/deepvariant/bin/call_variants' '{print \$2}' | sed 's/--outfile "[^"]*"//g' | sed 's/--examples "[^"]*"//g')
        """

    stub:
        """
        make_examples_args=""
        call_variants_args=""
        touch sorted.bam
        touch sorted.bam.bai
        """

}

process deepvariant_make_examples {

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(make_examples_args), val(call_variants_args), val(regions_of_interest)
        val ref
        val ref_index
        val parbed

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(call_variants_args), path('make_examples.*.gz{,.example_info.json}'), path('make_examples_call_variant_outputs.*.gz'), path('gvcf.*.gz')

    script:
        // define an optional string to pass regions of interest bed file
        def regions_of_interest_optional = file(regions_of_interest).name != 'NONE' ? "--regions $regions_of_interest" : ''
        """
        seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer make_examples \\
            --mode calling --ref "${ref}" --reads "${bam}" --sample_name "${sample_id}" ${regions_of_interest_optional} --examples "make_examples.tfrecord@${task.cpus}.gz" --gvcf "gvcf.tfrecord@${task.cpus}.gz" ${make_examples_args}
        """

    stub:
        """
        touch make_examples.tfrecord-00000-of-00104.gz
        touch make_examples.tfrecord-00000-of-00104.gz.example_info.json
        touch make_examples_call_variant_outputs.tfrecord-00000-of-00104.gz
        touch gvcf.tfrecord-00000-of-00104.gz
        """

}

process deepvariant_call_variants {

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(call_variants_args), path(make_examples_out), val(make_examples_call_variant_out), val(gvcf)

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(make_examples_call_variant_out), val(gvcf), path('call_variants_output*.gz')

    script:
        def matcher = make_examples_out[0].baseName =~ /^(.+)-\d{5}-of-(\d{5})$/
        def num_shards = matcher[0][2] as int
        """
        call_variants --outfile "call_variants_output.tfrecord.gz" --examples "make_examples.tfrecord@${num_shards}.gz" ${call_variants_args}
        """

    stub:
        """
        touch call_variants_output-00000-of-00016.tfrecord.gz
        """

}

process deepvariant_post_processing {

    // for when deeptrio selected as the snp/indel caller
    def reported_snp_indel_caller = "deepvariant"

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename ->
        if (params.snp_indel_caller != 'deeptrio') {
            return "$sample_id.$ref_name.$snp_indel_caller.$filename"
        } else {
            return "$sample_id.$ref_name.$reported_snp_indel_caller.$filename"
        }
    }, pattern: 'snp_indel.g.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path(make_examples_call_variant_out), path(gvcf), path(call_variants_out), val(family_position)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path('snp_indel.vcf.gz'), path('snp_indel.vcf.gz.tbi')
        tuple val(sample_id), val(family_id), val(family_position), path("${family_position}.sorted.bam"), path("${family_position}.sorted.bam.bai"), path("${family_position}_snp_indel.g.vcf.gz")
        tuple val(sample_id), path('snp_indel.g.vcf.gz'), path('snp_indel.g.vcf.gz.tbi')

    script:
        def matcher = gvcf[0].baseName =~ /^(.+)-\d{5}-of-(\d{5})$/
        def num_shards = matcher[0][2] as int
        """
        # postprocess_variants and vcf_stats_report stages in deepvariant
        postprocess_variants --ref "${ref}" --infile "call_variants_output.tfrecord.gz" --outfile "snp_indel.vcf.gz" --nonvariant_site_tfrecord_path "gvcf.tfrecord@${num_shards}.gz" --gvcf_outfile "snp_indel.g.vcf.gz" --cpus "${task.cpus}" --small_model_cvo_records "make_examples_call_variant_outputs.tfrecord@${num_shards}.gz" --sample_name "${sample_id}"
        vcf_stats_report --input_vcf "snp_indel.vcf.gz" --outfile_base "snp_indel"
        # tag bam and gvcf with family_position for downstream glnexus
        ln -s snp_indel.g.vcf.gz ${family_position}_snp_indel.g.vcf.gz
        ln -s $bam ${family_position}.sorted.bam
        ln -s $bam_index ${family_position}.sorted.bam.bai
        """

    stub:
        """
        touch snp_indel.vcf.gz
        touch snp_indel.vcf.gz.tbi
        touch snp_indel.g.vcf.gz
        touch snp_indel.g.vcf.gz.tbi
        """

}

process filter_ref_call {

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path(snp_indel_vcf), path(snp_indel_vcf_index)

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path('snp_indel.filtered.vcf.gz')

    script:
        """
        # filter out refcall variants
        bcftools view -f 'PASS' $snp_indel_vcf -o snp_indel.filtered.vcf.gz
        """
    stub:
        """
        touch snp_indel.filtered.vcf.gz
        """

}

process split_multiallele {

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path(snp_indel_vcf)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path('snp_indel.split.vcf.gz'), path('snp_indel.split.vcf.gz.tbi')

    script:
        """
        # run bcftools norm
        bcftools norm --threads ${task.cpus} -m -any -f $ref $snp_indel_vcf > snp_indel.split.unsorted.vcf
        # sort
        bcftools sort -o snp_indel.split.vcf snp_indel.split.unsorted.vcf
        # compress and index vcf
        bgzip -@ ${task.cpus} snp_indel.split.vcf
        tabix snp_indel.split.vcf.gz
        """

    stub:
        """
        touch snp_indel.split.vcf.gz
        touch snp_indel.split.vcf.gz.tbi
        """

}

process whatshap_phase {

    // for when deeptrio selected as the snp/indel caller
    def reported_snp_indel_caller = "deepvariant"

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename ->
        if (params.snp_indel_caller != 'deeptrio') {
            return "$sample_id.$ref_name.$snp_indel_caller.$filename"
        } else {
            return "$sample_id.$ref_name.$reported_snp_indel_caller.$filename"
        }
    }, pattern: 'snp_indel.phased.*'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path(snp_indel_split_vcf), path(snp_indel_split_vcf_index)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path('snp_indel.phased.vcf.gz'), path('snp_indel.phased.vcf.gz.tbi')
        tuple val(sample_id), val(family_id), path('snp_indel.phased.vcf.gz')
        tuple val(sample_id), val(family_id), path('snp_indel.phased.vcf.gz'), path('snp_indel.phased.vcf.gz.tbi'), path('snp_indel.phased.read_list.txt'), path('snp_indel.phased.stats.gtf')

    script:
        """
        # rename file for publishing purposes
        ln -sf $snp_indel_split_vcf snp_indel.vcf.gz
        ln -sf $snp_indel_split_vcf_index snp_indel.vcf.gz.tbi
        # run whatshap phase
        whatshap phase --reference $ref --output snp_indel.phased.vcf.gz --output-read-list snp_indel.phased.read_list.txt --sample $sample_id --ignore-read-groups snp_indel.vcf.gz $bam 
        # index vcf
        tabix snp_indel.phased.vcf.gz
        # run whatshap stats
        whatshap stats snp_indel.phased.vcf.gz --gtf snp_indel.phased.stats.gtf --sample $sample_id 
        """

    stub:
        """
        touch snp_indel.phased.vcf.gz
        touch snp_indel.phased.vcf.gz.tbi
        touch snp_indel.phased.read_list.txt
        touch snp_indel.phased.stats.gtf
        """

}

process whatshap_haplotag {

    def mapper_phaser = "minimap2.whatshap"

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$mapper_phaser.$filename" }, pattern: 'sorted.haplotagged.*'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path(snp_indel_split_vcf), path(snp_indel_split_vcf_index), val(family_position)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path('sorted.haplotagged.bam'), path('sorted.haplotagged.bam.bai')
        tuple val(sample_id), val(family_id), val(family_position), path("${family_position}.sorted.haplotagged.bam"), path("${family_position}.sorted.haplotagged.bam.bai")
        tuple val(sample_id), val(family_id), path('sorted.haplotagged.tsv')

    script:
        """
        # run whatshap haplotag
        whatshap haplotag --reference $ref --output sorted.haplotagged.bam --sample $sample_id --tag-supplementary --ignore-read-groups --output-threads ${task.cpus} --output-haplotag-list sorted.haplotagged.tsv $snp_indel_split_vcf $bam
        # index bam
        samtools index -@ ${task.cpus} sorted.haplotagged.bam
        # tag bam with family_position for downstream deeptrio
        ln -s sorted.haplotagged.bam ${family_position}.sorted.haplotagged.bam
        ln -s sorted.haplotagged.bam.bai ${family_position}.sorted.haplotagged.bam.bai
        """

    stub:
        """
        touch sorted.haplotagged.bam
        touch sorted.haplotagged.bam.bai
        touch sorted.haplotagged.tsv
        """

}

process deeptrio_dry_run {

    input:
        tuple val(proband_sample_id), val(proband_family_id), val(proband_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), val(proband_data_type)
        tuple val(father_sample_id), val(father_family_id), val(father_family_position), path(father_haplotagged_bam), path(father_haplotagged_bam_index), val(father_data_type)
        tuple val(mother_sample_id), val(mother_family_id), val(mother_family_position), path(mother_haplotagged_bam), path(mother_haplotagged_bam_index), val(mother_data_type)
        val ref
        val ref_index

    output:
        tuple val(proband_sample_id), val(proband_family_id), val(proband_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), val(proband_data_type)
        tuple val(father_sample_id), val(father_family_id), val(father_family_position), path(father_haplotagged_bam), path(father_haplotagged_bam_index), val(father_data_type)
        tuple val(mother_sample_id), val(mother_family_id), val(mother_family_position), path(mother_haplotagged_bam), path(mother_haplotagged_bam_index), val(mother_data_type)
        tuple env(make_examples_cs_args), env(make_examples_calling_args), env(call_variants_proband_args), env(call_variants_father_args), env(call_variants_mother_args)
        
    script:
        // conditionally define model type
        if (proband_data_type == 'ont') {
            model = 'ONT'
        }
        else if (proband_data_type == 'pacbio') {
            model = 'PACBIO'
        }
        """
        run_deeptrio --model_type=$model --ref=$ref --sample_name_child=$proband_sample_id --sample_name_parent1=$father_sample_id --sample_name_parent2=$mother_sample_id --reads_child=$proband_haplotagged_bam --reads_parent1=$father_haplotagged_bam --reads_parent2=$mother_haplotagged_bam --output_vcf_child=child.vcf.gz --output_vcf_parent1=parent1.vcf.gz --output_vcf_parent2=parent2.vcf.gz --output_gvcf_child=child.g.vcf.gz --output_gvcf_parent1=parent1.g.vcf.gz --output_gvcf_parent2=parent2.g.vcf.gz --dry_run=true > commands.txt
        make_examples_cs_args=\$(grep "/opt/deepvariant/bin/deeptrio/make_examples --mode candidate_sweep" commands.txt | awk -F'/opt/deepvariant/bin/deeptrio/make_examples' '{print \$2}' | sed 's/--ref "[^"]*"//g' | sed 's/--sample_name "[^"]*"//g' | sed 's/--reads "[^"]*"//g' | sed 's/--sample_name_parent1 "[^"]*"//g' | sed 's/--reads_parent1 "[^"]*"//g' | sed 's/--sample_name_parent2 "[^"]*"//g' | sed 's/--reads_parent2 "[^"]*"//g' | sed 's/--examples "[^"]*"//g' | sed 's/--candidate_positions "[^"]*"//g' | sed 's/--gvcf "[^"]*"//g')
        make_examples_calling_args=\$(grep "/opt/deepvariant/bin/deeptrio/make_examples --mode calling" commands.txt | awk -F'/opt/deepvariant/bin/deeptrio/make_examples' '{print \$2}' | sed 's/--ref "[^"]*"//g' | sed 's/--sample_name "[^"]*"//g' | sed 's/--reads "[^"]*"//g' | sed 's/--sample_name_parent1 "[^"]*"//g' | sed 's/--reads_parent1 "[^"]*"//g' | sed 's/--sample_name_parent2 "[^"]*"//g' | sed 's/--reads_parent2 "[^"]*"//g' | sed 's/--examples "[^"]*"//g' | sed 's/--candidate_positions "[^"]*"//g' | sed 's/--gvcf "[^"]*"//g')
        call_variants_proband_args=\$(grep "/opt/deepvariant/bin/call_variants" commands.txt | grep "child" | awk -F'/opt/deepvariant/bin/call_variants' '{print \$2}' | sed 's/--outfile "[^"]*"//g' | sed 's/--examples "[^"]*"//g')
        call_variants_father_args=\$(grep "/opt/deepvariant/bin/call_variants" commands.txt | grep "parent1" | awk -F'/opt/deepvariant/bin/call_variants' '{print \$2}' | sed 's/--outfile "[^"]*"//g' | sed 's/--examples "[^"]*"//g')
        call_variants_mother_args=\$(grep "/opt/deepvariant/bin/call_variants" commands.txt | grep "parent2" | awk -F'/opt/deepvariant/bin/call_variants' '{print \$2}' | sed 's/--outfile "[^"]*"//g' | sed 's/--examples "[^"]*"//g')
        """

    stub:
        """
        make_examples_cs_args=""
        make_examples_calling_args=""
        call_variants_proband_args=""
        call_variants_father_args=""
        call_variants_mother_args=""
        """

}

process deeptrio_make_examples {

    input:
        tuple val(proband_sample_id), val(proband_family_id), val(proband_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), val(proband_data_type)
        tuple val(father_sample_id), val(father_family_id), val(father_family_position), path(father_haplotagged_bam), path(father_haplotagged_bam_index), val(father_data_type)
        tuple val(mother_sample_id), val(mother_family_id), val(mother_family_position), path(mother_haplotagged_bam), path(mother_haplotagged_bam_index), val(mother_data_type)
        tuple val(make_examples_cs_args), val(make_examples_calling_args), val(call_variants_proband_args), val(call_variants_father_args), val(call_variants_mother_args)
        val ref
        val ref_index

    output:
        tuple val(proband_sample_id), val(proband_family_id), val(proband_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path('make_examples_child.*.gz'), path('gvcf_child.*.gz'), path('*.example_info.json'), val(call_variants_proband_args), emit: proband
        tuple val(father_sample_id), val(father_family_id), val(father_family_position), path(father_haplotagged_bam), path(father_haplotagged_bam_index), path('make_examples_parent1.*.gz'), path('gvcf_parent1.*.gz'), path('*.example_info.json'), val(call_variants_father_args), emit: father
        tuple val(mother_sample_id), val(mother_family_id), val(mother_family_position), path(mother_haplotagged_bam), path(mother_haplotagged_bam_index), path('make_examples_parent2.*.gz'), path('gvcf_parent2.*.gz'), path('*.example_info.json'), val(call_variants_mother_args), emit: mother

    script:
        """
        seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer make_examples \\
            --ref "${ref}" --sample_name "${proband_sample_id}" --reads "${proband_haplotagged_bam}" --sample_name_parent1 "${father_sample_id}" --reads_parent1 "${father_haplotagged_bam}" \\
            --sample_name_parent2 "${mother_sample_id}" --reads_parent2 "${mother_haplotagged_bam}" --examples "make_examples.tfrecord@${task.cpus}.gz" --gvcf "gvcf.tfrecord@${task.cpus}.gz" --candidate_positions "candidate_positions@${task.cpus}.gz" ${make_examples_cs_args}
        
        seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer make_examples \\
            --ref "${ref}" --sample_name "${proband_sample_id}" --reads "${proband_haplotagged_bam}" --sample_name_parent1 "${father_sample_id}" --reads_parent1 "${father_haplotagged_bam}" \\
            --sample_name_parent2 "${mother_sample_id}" --reads_parent2 "${mother_haplotagged_bam}" --examples "make_examples.tfrecord@${task.cpus}.gz" --gvcf "gvcf.tfrecord@${task.cpus}.gz" --candidate_positions "candidate_positions@${task.cpus}.gz" ${make_examples_calling_args}
        """

    stub:
        """
        touch make_examples_child.tfrecord-00000-of-00104.gz
        touch make_examples_parent1.tfrecord-00000-of-00104.gz
        touch make_examples_parent2.tfrecord-00000-of-00104.gz
        touch make_examples.tfrecord-00000-of-00104.gz.example_info.json
        touch gvcf_child.tfrecord-00000-of-00104.gz
        touch gvcf_parent1.tfrecord-00000-of-00104.gz
        touch gvcf_parent2.tfrecord-00000-of-00104.gz
        """

}

process deeptrio_call_variants {

    input:
        tuple val(sample_id), val(proband_family_id), val(family_position), path(haplotagged_bam), path(haplotagged_bam_index), path(make_examples), val(gvcf), path(example_info), val(call_variants_args)

    output:
        tuple val(sample_id), val(proband_family_id), val(family_position), path(haplotagged_bam), path(haplotagged_bam_index), path('*.gz'), val(gvcf)

    script:
        def matcher = make_examples[0].baseName =~ /^(.+)-\d{5}-of-(\d{5})$/
        def make_examples_name = matcher[0][1]
        def make_examples_num_shards = matcher[0][2] as int
        """
        call_variants --outfile "call_variants_output.tfrecord.gz" --examples "${make_examples_name}@${make_examples_num_shards}.gz" ${call_variants_args}
        """

    stub:
        """
        touch call_variants_output-00000-of-00016.tfrecord.gz
        """

}

process deeptrio_postprocessing {
    
    input:
        tuple val(sample_id), val(proband_family_id), val(family_position), path(haplotagged_bam), path(haplotagged_bam_index), path(call_variants), path(gvcf)
        val ref
        val ref_index
        
    output:
        tuple val(sample_id), val(proband_family_id), val(family_position), path(haplotagged_bam), path(haplotagged_bam_index), path("${family_position}_snp_indel.g.vcf.gz")

    script:
        def matcher = gvcf[0].baseName =~ /^(.+)-\d{5}-of-(\d{5})$/
        def gvcf_name = matcher[0][1]
        def gvcf_num_shards = matcher[0][2] as int
        """
        postprocess_variants --ref "${ref}" --sample_name "${sample_id}" --infile "call_variants_output.tfrecord.gz" --nonvariant_site_tfrecord_path "${gvcf_name}@${gvcf_num_shards}.gz" --cpus "${task.cpus}" --outfile "${family_position}_snp_indel.vcf.gz" --gvcf_outfile "${family_position}_snp_indel.g.vcf.gz"
        vcf_stats_report --input_vcf "${family_position}_snp_indel.vcf.gz" --outfile_base "${family_position}_snp_indel"
        """

    stub:
        """
        touch ${family_position}_snp_indel.vcf.gz
        touch ${family_position}_snp_indel.vcf.gz.tbi
        touch ${family_position}_snp_indel.g.vcf
        """

}

process somalier_duo {

    publishDir "$outdir/$proband_family_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$proband_family_id.$ref_name.$filename" }, pattern: 'somalier*'

    input:
        tuple val(proband_sample_id), val(proband_family_id), val(proband_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(proband_gvcf)
        tuple val(parent_sample_id), val(parent_family_id), val(parent_family_position), path(parent_haplotagged_bam), path(parent_haplotagged_bam_index), path(parent_gvcf)
        val ref
        val ref_index
        val sites
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(proband_sample_id), path('somalier.samples.tsv'), path('somalier.pairs.tsv'), path('somalier.html')

    script:
        """
        # run somalier extract
        SOMALIER_SAMPLE_NAME=$proband_sample_id somalier extract -d extracted --sites $sites -f $ref $proband_haplotagged_bam
        SOMALIER_SAMPLE_NAME=$parent_sample_id somalier extract -d extracted --sites $sites -f $ref $parent_haplotagged_bam
        # run somalier relate
        somalier relate extracted/*.somalier
        """

}

process somalier_trio {

    publishDir "$outdir/$proband_family_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$proband_family_id.$ref_name.$filename" }, pattern: 'somalier*'

    input:
        tuple val(proband_sample_id), val(proband_family_id), val(proband_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(proband_gvcf)
        tuple val(father_sample_id), val(father_family_id), val(father_family_position), path(father_haplotagged_bam), path(father_haplotagged_bam_index), path(father_gvcf)
        tuple val(mother_sample_id), val(mother_family_id), val(mother_family_position), path(mother_haplotagged_bam), path(mother_haplotagged_bam_index), path(mother_gvcf)
        val ref
        val ref_index
        val sites
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(proband_sample_id), path('somalier.samples.tsv'), path('somalier.pairs.tsv'), path('somalier.html')

    script:
        """
        # run somalier extract
        SOMALIER_SAMPLE_NAME=$proband_sample_id somalier extract -d extracted --sites $sites -f $ref $proband_haplotagged_bam
        SOMALIER_SAMPLE_NAME=$mother_sample_id somalier extract -d extracted --sites $sites -f $ref $mother_haplotagged_bam
        SOMALIER_SAMPLE_NAME=$father_sample_id somalier extract -d extracted --sites $sites -f $ref $father_haplotagged_bam
        # run somalier relate
        somalier relate extracted/*.somalier
        """

}

process glnexus_duo {

    input:
        tuple val(proband_sample_id), val(proband_family_id), val(proband_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(proband_gvcf)
        tuple val(parent_sample_id), val(parent_family_id), val(parent_family_position), path(parent_haplotagged_bam), path(parent_haplotagged_bam_index), path(parent_gvcf)

    output:
        tuple val(proband_sample_id), val(parent_sample_id), val(proband_family_id), val(parent_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(parent_haplotagged_bam), path(parent_haplotagged_bam_index), path('snp_indel.bcf')

    script:
        """
        glnexus_cli --config DeepVariant $proband_gvcf $parent_gvcf > snp_indel.bcf
        """

    stub:
        """
        touch snp_indel.bcf
        """

}

process glnexus_trio {

    input:
        tuple val(proband_sample_id), val(proband_family_id), val(proband_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(proband_gvcf)
        tuple val(father_sample_id), val(father_family_id), val(father_family_position), path(father_haplotagged_bam), path(father_haplotagged_bam_index), path(father_gvcf)
        tuple val(mother_sample_id), val(mother_family_id), val(mother_family_position), path(mother_haplotagged_bam), path(mother_haplotagged_bam_index), path(mother_gvcf)

    output:
        tuple val(proband_sample_id), val(father_sample_id), val(mother_sample_id), val(proband_family_id), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(father_haplotagged_bam), path(father_haplotagged_bam_index), path(mother_haplotagged_bam), path(mother_haplotagged_bam_index), path('snp_indel.bcf')

    script:
        """
        glnexus_cli --config DeepVariant $proband_gvcf $father_gvcf $mother_gvcf > snp_indel.bcf
        """

    stub:
        """
        touch snp_indel.bcf
        """

}

process glnexus_post_processing_duo {

    input:
        tuple val(proband_sample_id), val(parent_sample_id), val(proband_family_id), val(parent_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(parent_haplotagged_bam), path(parent_haplotagged_bam_index), path(snp_indel_bcf)

    output:
        tuple val(proband_sample_id), val(parent_sample_id), val(proband_family_id), val(parent_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(parent_haplotagged_bam), path(parent_haplotagged_bam_index), path('snp_indel.vcf.gz'), path('snp_indel.vcf.gz.tbi')

    script:
        """
        bcftools view snp_indel.bcf | bgzip -@ ${task.cpus} -c > snp_indel.vcf.gz
        tabix snp_indel.vcf.gz
        """

    stub:
        """
        touch snp_indel.vcf.gz
        touch snp_indel.vcf.gz.tbi
        """

}

process glnexus_post_processing_trio {

    input:
        tuple val(proband_sample_id), val(father_sample_id), val(mother_sample_id), val(proband_family_id), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(father_haplotagged_bam), path(father_haplotagged_bam_index), path(mother_haplotagged_bam), path(mother_haplotagged_bam_index), path(snp_indel_bcf)

    output:
        tuple val(proband_sample_id), val(father_sample_id), val(mother_sample_id), val(proband_family_id), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(father_haplotagged_bam), path(father_haplotagged_bam_index), path(mother_haplotagged_bam), path(mother_haplotagged_bam_index), path('snp_indel.vcf.gz'), path('snp_indel.vcf.gz.tbi')

    script:
        """
        bcftools view snp_indel.bcf | bgzip -@ ${task.cpus} -c > snp_indel.vcf.gz
        tabix snp_indel.vcf.gz
        """

    stub:
        """
        touch snp_indel.vcf.gz
        touch snp_indel.vcf.gz.tbi
        """

}

process split_multiallele_duo {

    input:
        tuple val(proband_sample_id), val(parent_sample_id), val(proband_family_id), val(parent_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(parent_haplotagged_bam), path(parent_haplotagged_bam_index), path(snp_indel_vcf), path(snp_indel_vcf_index)
        val ref
        val ref_index

    output:
        tuple val(proband_sample_id), val(parent_sample_id), val(proband_family_id), val(parent_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(parent_haplotagged_bam), path(parent_haplotagged_bam_index), path('snp_indel.split.vcf.gz'), path('snp_indel.split.vcf.gz.tbi')

    script:
        """
        # run bcftools norm
        bcftools norm --threads ${task.cpus} -m -any -f $ref snp_indel.vcf.gz > snp_indel.split.vcf
        # compress and index vcf
        bgzip -@ ${task.cpus} snp_indel.split.vcf
        tabix snp_indel.split.vcf.gz
        """

    stub:
        """
        touch snp_indel.split.vcf.gz
        touch snp_indel.split.vcf.gz.tbi
        """

}

process split_multiallele_trio {

    input:
        tuple val(proband_sample_id), val(father_sample_id), val(mother_sample_id), val(proband_family_id), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(father_haplotagged_bam), path(father_haplotagged_bam_index), path(mother_haplotagged_bam), path(mother_haplotagged_bam_index), path(snp_indel_vcf), path(snp_indel_vcf_index)
        val ref
        val ref_index

    output:
        tuple val(proband_sample_id), val(father_sample_id), val(mother_sample_id), val(proband_family_id), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(father_haplotagged_bam), path(father_haplotagged_bam_index), path(mother_haplotagged_bam), path(mother_haplotagged_bam_index), path('snp_indel.split.vcf.gz'), path('snp_indel.split.vcf.gz.tbi')

    script:
        """
        # run bcftools norm
        bcftools norm --threads ${task.cpus} -m -any -f $ref snp_indel.vcf.gz > snp_indel.split.vcf
        # compress and index vcf
        bgzip -@ ${task.cpus} snp_indel.split.vcf
        tabix snp_indel.split.vcf.gz
        """

    stub:
        """
        touch snp_indel.split.vcf.gz
        touch snp_indel.split.vcf.gz.tbi
        """

}

process whatshap_phase_duo {

    publishDir "$outdir/$proband_family_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$proband_family_id.$ref_name.$snp_indel_caller.$filename" }, pattern: 'snp_indel.phased.*'

    input:
        tuple val(proband_sample_id), val(parent_sample_id), val(proband_family_id), val(parent_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(parent_haplotagged_bam), path(parent_haplotagged_bam_index), path(snp_indel_split_vcf), path(snp_indel_split_vcf_index)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller

    output:
        tuple val(proband_sample_id), val(proband_family_id), path('snp_indel.phased.vcf.gz')
        tuple val(proband_family_id), path('snp_indel.phased.vcf.gz'), path('snp_indel.phased.vcf.gz.tbi'), path('snp_indel.phased.read_list.txt'), path('snp_indel.phased.stats.gtf')

    script:
        """
        # handle missing parent
        if [[ $parent_family_position == 'father' ]]; then
            FATHER_SAMPLE_ID=$parent_sample_id
            MOTHER_SAMPLE_ID="0"
        elif [[ $parent_family_position == 'mother' ]]; then
            FATHER_SAMPLE_ID="0"
            MOTHER_SAMPLE_ID=$parent_sample_id
        fi
        # rename file for publishing purposes
        ln -sf $snp_indel_split_vcf snp_indel.vcf.gz
        printf "$proband_family_id\t$proband_sample_id\t\${FATHER_SAMPLE_ID}\t\${MOTHER_SAMPLE_ID}\t0\t1\n" > pedigree.ped
        # run whatshap phase
        whatshap phase --reference $ref --output snp_indel.phased.vcf.gz --output-read-list snp_indel.phased.read_list.txt --ped pedigree.ped snp_indel.vcf.gz $proband_haplotagged_bam $parent_haplotagged_bam
        # index vcf
        tabix snp_indel.phased.vcf.gz
        # run whatshap stats
        whatshap stats snp_indel.phased.vcf.gz --gtf snp_indel.phased.stats.gtf
        """

    stub:
        """
        touch snp_indel.phased.vcf.gz
        touch snp_indel.phased.vcf.gz.tbi
        touch snp_indel.phased.read_list.txt
        touch snp_indel.phased.stats.gtf
        """

}

process whatshap_phase_trio {

    publishDir "$outdir/$proband_family_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$proband_family_id.$ref_name.$snp_indel_caller.$filename" }, pattern: 'snp_indel.phased.*'

    input:
        tuple val(proband_sample_id), val(father_sample_id), val(mother_sample_id), val(proband_family_id), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path(father_haplotagged_bam), path(father_haplotagged_bam_index), path(mother_haplotagged_bam), path(mother_haplotagged_bam_index), path(snp_indel_split_vcf), path(snp_indel_split_vcf_index)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller

    output:
        tuple val(proband_sample_id), val(proband_family_id), path('snp_indel.phased.vcf.gz')
        tuple val(proband_family_id), path('snp_indel.phased.vcf.gz'), path('snp_indel.phased.vcf.gz.tbi'), path('snp_indel.phased.read_list.txt'), path('snp_indel.phased.stats.gtf')

    script:
        """
        # rename file for publishing purposes
        ln -sf $snp_indel_split_vcf snp_indel.vcf.gz
        # create pedigree file
        printf "$proband_family_id\t$proband_sample_id\t$father_sample_id\t$mother_sample_id\t0\t1\n" > pedigree.ped
        # run whatshap phase
        whatshap phase --reference $ref --output snp_indel.phased.vcf.gz --output-read-list snp_indel.phased.read_list.txt --ped pedigree.ped snp_indel.vcf.gz $proband_haplotagged_bam $father_haplotagged_bam $mother_haplotagged_bam
        # index vcf
        tabix snp_indel.phased.vcf.gz
        # run whatshap stats
        whatshap stats snp_indel.phased.vcf.gz --gtf snp_indel.phased.stats.gtf
        """

    stub:
        """
        touch snp_indel.phased.vcf.gz
        touch snp_indel.phased.vcf.gz.tbi
        touch snp_indel.phased.read_list.txt
        touch snp_indel.phased.stats.gtf
        """

}

process vep_snp_indel {

    publishDir { task ->
        if (mode in ['singleton', 'NONE']) {
            return "$outdir/$family_id/$outdir2/$sample_id"
        } else {
            return "$outdir/$family_id/$outdir2"
        }
    }, mode: 'copy', overwrite: true, saveAs: { filename ->
        if (mode in ['singleton', 'NONE']) {
            return "$sample_id.$ref_name.$snp_indel_caller.$filename"
        } else {
            return "$family_id.$ref_name.$snp_indel_caller.$filename"
        }
    }, pattern: 'snp_indel.phased.annotated.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(snp_indel_split_phased_vcf)
        val ref
        val ref_index
        val vep_db
        val revel_db
        val gnomad_db
        val clinvar_db
        val cadd_snv_db
        val cadd_indel_db
        val spliceai_snv_db
        val spliceai_indel_db
        val alphamissense_db
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller
        val mode

    output:
        tuple val(sample_id), val(family_id), path('snp_indel.phased.annotated.vcf.gz'), path('snp_indel.phased.annotated.vcf.gz.tbi')

    script:
        """
        # run vep
        vep -i $snp_indel_split_phased_vcf -o snp_indel.phased.annotated.vcf.gz --format vcf --vcf --fasta $ref --dir $vep_db --assembly GRCh38 --species homo_sapiens --cache --offline --merged --sift b --polyphen b --symbol --hgvs --hgvsg --uploaded_allele --check_existing --filter_common --distance 0 --nearest gene --canonical --mane --pick --fork ${task.cpus} --no_stats --compress_output bgzip --dont_skip \
        --plugin REVEL,file=$revel_db --custom file=$gnomad_db,short_name=gnomAD,format=vcf,type=exact,fields=AF_joint%AF_exomes%AF_genomes%nhomalt_joint%nhomalt_exomes%nhomalt_genomes \
        --custom file=$clinvar_db,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG \
        --plugin CADD,snv=$cadd_snv_db,indels=$cadd_indel_db \
        --plugin SpliceAI,snv=$spliceai_snv_db,indel=$spliceai_indel_db \
        --plugin AlphaMissense,file=$alphamissense_db
        # index vcf
        tabix snp_indel.phased.annotated.vcf.gz
        """

    stub:
        """
        touch snp_indel.phased.annotated.vcf.gz
        touch snp_indel.phased.annotated.vcf.gz.tbi
        """

}

process minimod {

    def software = "minimod"

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$software.$filename" }, pattern: 'modfreqs_*.b*'

    input:
        tuple val(sample_id), val(family_id), path(haplotagged_bam), path(haplotagged_bam_index), val(data_type)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path('modfreqs_hap1.bw'), path('modfreqs_hap1.bed'), path('modfreqs_hap2.bw'), path('modfreqs_hap2.bed'), path('modfreqs_combined.bw'), path('modfreqs_combined.bed'), path('modfreqs_unphased.bed'), optional: true

    script:
        """
        # run minimod
        minimod mod-freq $ref $haplotagged_bam -t ${task.cpus} --haplotypes -o modfreqs.tmp.bed
        # sort
        awk 'NR > 1 { print }' modfreqs.tmp.bed | sort -k1,1 -k2,2n > modfreqs.bed
        if [ -s modfreqs.bed ]; then
            # seperate haplotypes
            for FILE in modfreqs_hap1 modfreqs_hap2 modfreqs_combined modfreqs_unphased; do head -n1 modfreqs.tmp.bed > \${FILE}.bed; done
            awk '\$9==1' modfreqs.bed >> modfreqs_hap1.bed
            awk '\$9==2' modfreqs.bed >> modfreqs_hap2.bed
            awk '\$9=="*"' modfreqs.bed >> modfreqs_combined.bed
            awk '\$9==0' modfreqs.bed >> modfreqs_unphased.bed
            # generate bigwig
            cut -f1,2 $ref_index > chrom.sizes
            for FILE in modfreqs_hap1 modfreqs_hap2 modfreqs_combined; do awk 'NR > 1 { print }' \${FILE}.bed | cut -f1-3,7 > \${FILE}.formatted.bed; done
            for FILE in modfreqs_hap1 modfreqs_hap2 modfreqs_combined; do bedGraphToBigWig \${FILE}.formatted.bed chrom.sizes \${FILE}.bw; done
        fi
        """

    stub:
        """
        touch modfreqs_hap1.bed
        touch modfreqs_hap1.bw
        touch modfreqs_hap2.bed
        touch modfreqs_hap2.bw
        touch modfreqs_combined.bed
        touch modfreqs_combined.bw
        """

}

process sniffles {

    def sv_caller = "sniffles"

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$sv_caller.$filename" }, pattern: 'sv.phased.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(haplotagged_bam), path(haplotagged_bam_index), val(family_position)
        val ref
        val ref_index
        val tandem_repeat
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path('sv.phased.vcf.gz')
        tuple val(sample_id), val(family_id), path('sv.phased.vcf.gz'), path('sv.phased.vcf.gz.tbi')
        tuple val(sample_id), val(family_id), val(family_position), path("${family_position}.sv.phased.vcf.gz"), path("${family_position}.sorted.haplotagged.bam")

    script:
        // define a string to optionally pass tandem repeat bed file
        def tandem_repeat_optional = file(tandem_repeat).name != 'NONE' ? "--tandem-repeats $tandem_repeat" : ''
        """
        # run sniffles
        sniffles \
        --reference $ref --input $haplotagged_bam --threads ${task.cpus} --sample-id $sample_id --vcf sv.phased.vcf.gz --output-rnames --minsvlen 50 --phase $tandem_repeat_optional
        # tag vcf and bam with family_position for downstream jasmine
        ln -s sv.phased.vcf.gz ${family_position}.sv.phased.vcf.gz
        ln -s sorted.haplotagged.bam ${family_position}.sorted.haplotagged.bam
        """

    stub:
        """
        touch sv.phased.vcf.gz
        touch sv.phased.vcf.gz.tbi
        """

}

process cutesv {

    def sv_caller = "cutesv"

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$sv_caller.$filename" }, pattern: 'sv.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(haplotagged_bam), path(haplotagged_bam_index), val(data_type), val(family_position)
        val ref
        val ref_index
        val tandem_repeat
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path('sv.vcf.gz')
        tuple val(sample_id), val(family_id), path('sv.vcf.gz'), path('sv.vcf.gz.tbi')
        tuple val(sample_id), val(family_id), val(family_position), path("${family_position}.sv.vcf.gz"), path("${family_position}.sorted.haplotagged.bam")

    script:
        // define platform specific settings
        if (data_type == 'ont') {
            settings = '--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3'
        }
        else if (data_type == 'pacbio') {
            settings = '--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5'
        }
        """
        # run cuteSV
        cuteSV $haplotagged_bam $ref sv.vcf ./ --sample ${sample_id} -t ${task.cpus} --genotype --report_readid --min_size 50 $settings
        # compress and index vcf
        bgzip -@ ${task.cpus} sv.vcf
        tabix sv.vcf.gz
        # tag vcf and bam with family_position for downstream jasmine
        ln -s sv.vcf.gz ${family_position}.sv.vcf.gz
        ln -s sorted.haplotagged.bam ${family_position}.sorted.haplotagged.bam
        """

    stub:
        """
        touch sv.vcf.gz
        touch sv.vcf.gz.tbi
        """

}

process jasmine_sniffles_duo {

    def sv_caller_merger = "sniffles.jasmine"

    publishDir "$outdir/$proband_family_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$proband_family_id.$ref_name.$sv_caller_merger.$filename" }, pattern: 'sv.phased.vcf*'

    input:
        tuple val(proband_sample_id), val(proband_family_id), val(proband_family_position), path(proband_sv_phased_vcf), path(proband_haplotagged_bam), val(proband_data_type)
        tuple val(parent_sample_id), val(parent_family_id), val(parent_family_position), path(parent_sv_phased_vcf), path(parent_haplotagged_bam)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(proband_sample_id), val(proband_family_id), path('sv.phased.vcf.gz')
        tuple val(proband_family_id), path('sv.phased.vcf.gz'), path('sv.phased.vcf.gz.tbi')

    script:
        // conditionally define iris arguments
        // as default, iris will pass minimap -x map-ont
        // the --pacbio flag passed to iris will pass minimap -x map-pb
        // these are the only two options iris makes available for minimaps -x argument, so I can't use lr:hq and map-hifi boo
        if (proband_data_type == 'ont') {
            iris_args = '--run_iris iris_args=min_ins_length=20,--rerunracon,--keep_long_variants'
        }
        else if (proband_data_type == 'pacbio') {
            iris_args = '--run_iris iris_args=min_ins_length=20,--rerunracon,--keep_long_variants,--pacbio'
        }
        """
        # unzip vcfs
        gunzip -c $proband_sv_phased_vcf > proband.sv.phased.vcf
        gunzip -c $parent_sv_phased_vcf > parent.sv.phased.vcf
        # create file lists
        realpath proband.sv.phased.vcf >> vcfs.txt
        realpath parent.sv.phased.vcf >> vcfs.txt
        realpath $proband_haplotagged_bam >> bams.txt
        realpath $parent_haplotagged_bam >> bams.txt
        # run jasmine
        jasmine threads=${task.cpus} out_dir=./ genome_file=$ref file_list=vcfs.txt bam_list=bams.txt out_file=sv.phased.tmp.vcf min_support=1 --mark_specific spec_reads=7 spec_len=20 --pre_normalize --output_genotypes --clique_merging --dup_to_ins --normalize_type --require_first_sample --default_zero_genotype $iris_args
        # fix vcf header (remove prefix to sample names that jasmine adds)
        grep '##' sv.phased.tmp.vcf > sv.phased.vcf
        grep '#CHROM' sv.phased.tmp.vcf | sed 's/0_//g' | sed 's/1_//g' >> sv.phased.vcf
        grep -v '#' sv.phased.tmp.vcf >> sv.phased.vcf
        bcftools sort sv.phased.vcf -o sv.phased.vcf
        # compress and index vcf
        bgzip -@ ${task.cpus} sv.phased.vcf
        tabix sv.phased.vcf.gz
        """

    stub:
        """
        touch sv.phased.vcf.gz
        touch sv.phased.vcf.gz.tbi
        """

}

process jasmine_cutesv_duo {

    def sv_caller_merger = "cutesv.jasmine"

    publishDir "$outdir/$proband_family_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$proband_family_id.$ref_name.$sv_caller_merger.$filename" }, pattern: 'sv.vcf*'

    input:
        tuple val(proband_sample_id), val(proband_family_id), val(proband_family_position), path(proband_sv_vcf), path(proband_haplotagged_bam), val(proband_data_type)
        tuple val(parent_sample_id), val(parent_family_id), val(parent_family_position), path(parent_sv_vcf), path(parent_haplotagged_bam)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(proband_sample_id), val(proband_family_id), path('sv.vcf.gz')
        tuple val(proband_family_id), path('sv.vcf.gz'), path('sv.vcf.gz.tbi')

    script:
        // conditionally define iris arguments
        // as default, iris will pass minimap -x map-ont
        // the --pacbio flag passed to iris will pass minimap -x map-pb
        // these are the only two options iris makes available for minimaps -x argument, so I can't use lr:hq and map-hifi boo
        if (proband_data_type == 'ont') {
            iris_args = '--run_iris iris_args=min_ins_length=20,--rerunracon,--keep_long_variants'
        }
        else if (proband_data_type == 'pacbio') {
            iris_args = '--run_iris iris_args=min_ins_length=20,--rerunracon,--keep_long_variants,--pacbio'
        }
        """
        # unzip vcfs
        gunzip -c $proband_sv_vcf > proband.sv.vcf
        gunzip -c $parent_sv_vcf > parent.sv.vcf
        # create file lists
        realpath proband.sv.vcf >> vcfs.txt
        realpath parent.sv.vcf >> vcfs.txt
        realpath $proband_haplotagged_bam >> bams.txt
        realpath $parent_haplotagged_bam >> bams.txt
        # run jasmine
        jasmine threads=${task.cpus} out_dir=./ genome_file=$ref file_list=vcfs.txt bam_list=bams.txt out_file=sv.tmp.vcf min_support=1 --mark_specific spec_reads=7 spec_len=20 --pre_normalize --output_genotypes --clique_merging --dup_to_ins --normalize_type --require_first_sample --default_zero_genotype $iris_args
        # fix vcf header (remove prefix to sample names that jasmine adds)
        grep '##' sv.tmp.vcf > sv.vcf
        grep '#CHROM' sv.tmp.vcf | sed 's/0_//g' | sed 's/1_//g' >> sv.vcf
        grep -v '#' sv.tmp.vcf >> sv.vcf
        bcftools sort sv.vcf -o sv.vcf
        # compress and index vcf
        bgzip -@ ${task.cpus} sv.vcf
        tabix sv.vcf.gz
        """

    stub:
        """
        touch sv.vcf.gz
        touch sv.vcf.gz.tbi
        """

}

process jasmine_sniffles_trio {

    def sv_caller_merger = "sniffles.jasmine"

    publishDir "$outdir/$proband_family_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$proband_family_id.$ref_name.$sv_caller_merger.$filename" }, pattern: 'sv.phased.vcf*'

    input:
        tuple val(proband_sample_id), val(proband_family_id), val(proband_family_position), path(proband_sv_phased_vcf), path(proband_haplotagged_bam), val(proband_data_type)
        tuple val(father_sample_id), val(father_family_id), val(father_family_position), path(father_sv_phased_vcf), path(father_haplotagged_bam)
        tuple val(mother_sample_id), val(mother_family_id), val(mother_family_position), path(mother_sv_phased_vcf), path(mother_haplotagged_bam)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(proband_sample_id), val(proband_family_id), path('sv.phased.vcf.gz')
        tuple val(proband_family_id), path('sv.phased.vcf.gz'), path('sv.phased.vcf.gz.tbi')

    script:
        // conditionally define iris arguments
        // as default, iris will pass minimap -x map-ont
        // the --pacbio flag passed to iris will pass minimap -x map-pb
        // these are the only two options iris makes available for minimaps -x argument, so I can't use lr:hq and map-hifi boo
        if (proband_data_type == 'ont') {
            iris_args = '--run_iris iris_args=min_ins_length=20,--rerunracon,--keep_long_variants'
        }
        else if (proband_data_type == 'pacbio') {
            iris_args = '--run_iris iris_args=min_ins_length=20,--rerunracon,--keep_long_variants,--pacbio'
        }
        """
        # unzip vcfs
        gunzip -c $proband_sv_phased_vcf > proband.sv.phased.vcf
        gunzip -c $father_sv_phased_vcf > father.sv.phased.vcf
        gunzip -c $mother_sv_phased_vcf > mother.sv.phased.vcf
        # create file lists
        realpath proband.sv.phased.vcf >> vcfs.txt
        realpath father.sv.phased.vcf >> vcfs.txt
        realpath mother.sv.phased.vcf >> vcfs.txt
        realpath $proband_haplotagged_bam >> bams.txt
        realpath $father_haplotagged_bam >> bams.txt
        realpath $mother_haplotagged_bam >> bams.txt
        # run jasmine
        jasmine threads=${task.cpus} out_dir=./ genome_file=$ref file_list=vcfs.txt bam_list=bams.txt out_file=sv.phased.tmp.vcf min_support=1 --mark_specific spec_reads=7 spec_len=20 --pre_normalize --output_genotypes --clique_merging --dup_to_ins --normalize_type --require_first_sample --default_zero_genotype $iris_args
        # fix vcf header (remove prefix to sample names that jasmine adds)
        grep '##' sv.phased.tmp.vcf > sv.phased.vcf
        grep '#CHROM' sv.phased.tmp.vcf | sed 's/0_//g' | sed 's/1_//g' | sed 's/2_//g' >> sv.phased.vcf
        grep -v '#' sv.phased.tmp.vcf >> sv.phased.vcf
        bcftools sort sv.phased.vcf -o sv.phased.vcf
        # compress and index vcf
        bgzip -@ ${task.cpus} sv.phased.vcf
        tabix sv.phased.vcf.gz
        """

    stub:
        """
        touch sv.phased.vcf.gz
        touch sv.phased.vcf.gz.tbi
        """

}

process jasmine_cutesv_trio {

    def sv_caller_merger = "cutesv.jasmine"

    publishDir "$outdir/$proband_family_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$proband_family_id.$ref_name.$sv_caller_merger.$filename" }, pattern: 'sv.vcf*'

    input:
        tuple val(proband_sample_id), val(proband_family_id), val(proband_family_position), path(proband_sv_vcf), path(proband_haplotagged_bam), val(proband_data_type)
        tuple val(father_sample_id), val(father_family_id), val(father_family_position), path(father_sv_vcf), path(father_haplotagged_bam)
        tuple val(mother_sample_id), val(mother_family_id), val(mother_family_position), path(mother_sv_vcf), path(mother_haplotagged_bam)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(proband_sample_id), val(proband_family_id), path('sv.vcf.gz')
        tuple val(proband_family_id), path('sv.vcf.gz'), path('sv.vcf.gz.tbi')

    script:
        // conditionally define iris arguments
        // as default, iris will pass minimap -x map-ont
        // the --pacbio flag passed to iris will pass minimap -x map-pb
        // these are the only two options iris makes available for minimaps -x argument, so I can't use lr:hq and map-hifi boo
        if (proband_data_type == 'ont') {
            iris_args = '--run_iris iris_args=min_ins_length=20,--rerunracon,--keep_long_variants'
        }
        else if (proband_data_type == 'pacbio') {
            iris_args = '--run_iris iris_args=min_ins_length=20,--rerunracon,--keep_long_variants,--pacbio'
        }
        """
        # unzip vcfs
        gunzip -c $proband_sv_vcf > proband.sv.vcf
        gunzip -c $father_sv_vcf > father.sv.vcf
        gunzip -c $mother_sv_vcf > mother.sv.vcf
        # create file lists
        realpath proband.sv.vcf >> vcfs.txt
        realpath father.sv.vcf >> vcfs.txt
        realpath mother.sv.vcf >> vcfs.txt
        realpath $proband_haplotagged_bam >> bams.txt
        realpath $father_haplotagged_bam >> bams.txt
        realpath $mother_haplotagged_bam >> bams.txt
        # run jasmine
        jasmine threads=${task.cpus} out_dir=./ genome_file=$ref file_list=vcfs.txt bam_list=bams.txt out_file=sv.tmp.vcf min_support=1 --mark_specific spec_reads=7 spec_len=20 --pre_normalize --output_genotypes --clique_merging --dup_to_ins --normalize_type --require_first_sample --default_zero_genotype $iris_args
        # fix vcf header (remove prefix to sample names that jasmine adds)
        grep '##' sv.tmp.vcf > sv.vcf
        grep '#CHROM' sv.tmp.vcf | sed 's/0_//g' | sed 's/1_//g' | sed 's/2_//g' >> sv.vcf
        grep -v '#' sv.tmp.vcf >> sv.vcf
        bcftools sort sv.vcf -o sv.vcf
        # compress and index vcf
        bgzip -@ ${task.cpus} sv.vcf
        tabix sv.vcf.gz
        """

    stub:
        """
        touch sv.vcf.gz
        touch sv.vcf.gz.tbi
        """

}

process vep_sniffles_sv {

    def sv_caller = "sniffles"
    def sv_caller_merger = "sniffles.jasmine"

    publishDir { task ->
        if (mode in ['singleton', 'NONE']) {
            return "$outdir/$family_id/$outdir2/$sample_id"
        } else {
            return "$outdir/$family_id/$outdir2"
        }
    }, mode: 'copy', overwrite: true, saveAs: { filename ->
        if (mode in ['singleton', 'NONE']) {
            return "$sample_id.$ref_name.$sv_caller.$filename"
        } else {
            return "$family_id.$ref_name.$sv_caller_merger.$filename"
        }
    }, pattern: 'sv.phased.annotated.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(sv_phased_vcf)
        val ref
        val ref_index
        val vep_db
        val gnomad_db
        val cadd_sv_db
        val outdir
        val outdir2
        val ref_name
        val mode

    output:
        tuple val(sample_id), val(family_id), path('sv.phased.annotated.vcf.gz'), path('sv.phased.annotated.vcf.gz.tbi')

    script:
        """
        # run vep
        vep -i $sv_phased_vcf -o sv.phased.annotated.vcf.gz --format vcf --vcf --fasta $ref --dir $vep_db --assembly GRCh38 --species homo_sapiens --cache --offline --merged --sift b --polyphen b --symbol --hgvs --hgvsg  --uploaded_allele --check_existing --filter_common --distance 0 --nearest gene --canonical --mane --pick --fork ${task.cpus} --no_stats --compress_output bgzip --dont_skip --plugin CADD,sv=$cadd_sv_db
        # index vcf
        tabix sv.phased.annotated.vcf.gz
        """

    stub:
        """
        touch sv.phased.annotated.vcf.gz
        touch sv.phased.annotated.vcf.gz.tbi
        """

}

process vep_cutesv_sv {

    def sv_caller = "cutesv"
    def sv_caller_merger = "cutesv.jasmine"

    publishDir { task ->
        if (mode in ['singleton', 'NONE']) {
            return "$outdir/$family_id/$outdir2/$sample_id"
        } else {
            return "$outdir/$family_id/$outdir2"
        }
    }, mode: 'copy', overwrite: true, saveAs: { filename ->
        if (mode in ['singleton', 'NONE']) {
            return "$sample_id.$ref_name.$sv_caller.$filename"
        } else {
            return "$family_id.$ref_name.$sv_caller_merger.$filename"
        }
    }, pattern: 'sv.annotated.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(sv_vcf)
        val ref
        val ref_index
        val vep_db
        val gnomad_db
        val cadd_sv_db
        val outdir
        val outdir2
        val ref_name
        val mode

    output:
        tuple val(sample_id), val(family_id), path('sv.annotated.vcf.gz'), path('sv.annotated.vcf.gz.tbi')

    script:
        """
        # run vep
        vep -i $sv_vcf -o sv.annotated.vcf.gz --format vcf --vcf --fasta $ref --dir $vep_db --assembly GRCh38 --species homo_sapiens --cache --offline --merged --sift b --polyphen b --symbol --hgvs --hgvsg --uploaded_allele --check_existing --filter_common --distance 0 --nearest gene --canonical --mane --pick --fork ${task.cpus} --no_stats --compress_output bgzip --dont_skip --plugin CADD,sv=$cadd_sv_db
        # index vcf
        tabix sv.annotated.vcf.gz
        """

    stub:
        """
        touch sv.annotated.vcf.gz
        touch sv.annotated.vcf.gz.tbi
        """

}

workflow {

    // grab parameters
    in_data = "${params.in_data}".trim()
    sex = "${params.sex}".trim()
    parbed = "${params.parbed}".trim()
    in_data_format = "${params.in_data_format}".trim()
    in_data_format_override = "${params.in_data_format_override}".trim()
    ref = "${params.ref}".trim()
    ref_index = "${params.ref_index}".trim()
    tandem_repeat = "${params.tandem_repeat}".trim()
    mode = "${params.mode}".trim()
    snp_indel_caller = "${params.snp_indel_caller}".trim()
    sv_caller = "${params.sv_caller}".trim()
    annotate = "${params.annotate}".trim()
    haploidaware = "${params.haploidaware}".trim()
    annotate_override = "${params.annotate_override}".trim()
    calculate_depth = "${params.calculate_depth}".trim()
    analyse_base_mods = "${params.analyse_base_mods}".trim()
    check_relatedness = "${params.check_relatedness}".trim()
    sites = "${params.sites}".trim()
    outdir = "${params.outdir}".trim()
    outdir2 = "${params.outdir2}".trim()
    vep_db = "${params.vep_db}".trim()
    revel_db = "${params.revel_db}".trim()
    gnomad_db = "${params.gnomad_db}".trim()
    clinvar_db = "${params.clinvar_db}".trim()
    cadd_snv_db = "${params.cadd_snv_db}".trim()
    cadd_indel_db = "${params.cadd_indel_db}".trim()
    cadd_sv_db = "${params.cadd_sv_db}".trim()
    spliceai_snv_db = "${params.spliceai_snv_db}".trim()
    spliceai_indel_db = "${params.spliceai_indel_db}".trim()
    alphamissense_db = "${params.alphamissense_db}".trim()

    // check user provided parameters
    if (!file(in_data).exists()) {
        exit 1, "In data csv file does not exist, 'in_data = ${in_data}' provided."
    }
    if (!(in_data_format in ['ubam_fastq', 'aligned_bam', 'snp_indel_vcf', 'sv_vcf'])) {
        exit 1, "In data format should be 'ubam_fastq', 'aligned_bam', 'snp_indel_vcf' or 'sv_vcf', in_data_format = '${in_data_format}' provided."
    }
    if (!file(ref).exists()) {
        exit 1, "Reference genome file does not exist, ref = '${ref}' provided."
    }
    if (!file(ref_index).exists()) {
        exit 1, "Reference genome index file does not exist, ref_index = '${ref_index}' provided."
    }
    if (!(haploidaware in ['yes', 'no'])) {
        exit 1, "haploidaware must be either 'yes' or 'no', haploidaware = '${haploidaware}' provided."
    }
    if (haploidaware == "yes") {
        if (sex != "XY") {
            exit 1, "Haploid-aware mode is only supported when sex = 'XY', sex = '${sex}' provided."
        }
        if (parbed == "NONE") {
            exit 1, "In haploid-aware mode, provide a valid PAR BED file, parbed = 'NONE' provided."
        }
        if (!file(parbed).exists()) {
            exit 1, "PAR BED file does not exist, parbed = '${parbed}' provided."
        }
    }
    else if (haploidaware == "no") {
        // sex can be anything, including unset
        if (parbed != "NONE") {
            exit 1, "In haploid-aware mode, set the PAR BED file to 'NONE', haploidaware = '${haploidaware}' and parbed = '${parbed}' provided."
        }
    }
    if (!file(tandem_repeat).exists()) {
        exit 1, "Tandem repeat bed file does not exist, tandem_repeat = '${tandem_repeat}' provided. Set to 'NONE' if not required."
    }
    if (in_data_format in ['ubam_fastq', 'aligned_bam']) {
        if (!(mode in ['singleton', 'duo', 'trio'])) {
            exit 1, "Mode should be 'singleton', 'duo' or 'trio', mode = '${mode}' provided."
        }
        if (mode == 'singleton' && !(snp_indel_caller in ['clair3', 'deepvariant'])) {
            exit 1, "When in singleton mode, the SNP/indel caller must be either 'clair3' or 'deepvariant', mode = '${mode}' and snp_indel_caller = '${snp_indel_caller}' provided."
        }
        if (mode == 'duo' && snp_indel_caller != 'deepvariant') {
            exit 1, "When in duo mode, the SNP/indel caller must be 'deepvariant', mode = '${mode}' and snp_indel_caller = '${snp_indel_caller}' provided."
        }
        if (mode == 'trio' && snp_indel_caller != 'deeptrio') {
            exit 1, "When in trio mode, the SNP/indel caller must be 'deeptrio', mode = '${mode}' and snp_indel_caller = '${snp_indel_caller}' provided."
        }
        if (!(sv_caller in ['cutesv', 'sniffles', 'both'])) {
            exit 1, "SV calling software should be 'sniffles', 'cutesv', or 'both', sv_caller = '${sv_caller}' provided."
        }
    }
    if (!(annotate in ['yes', 'no'])) {
        exit 1, "Choice to annotate should be either 'yes' or 'no', annotate = '${annotate}' provided."
    }
    if (annotate == 'yes') {
        if (!file(vep_db).exists()) {
            exit 1, "VEP cache directory does not exist, vep_db = '${vep_db}' provided."
        }
        if (!file(revel_db).exists()) {
            exit 1, "REVEL database file does not exist, revel_db = '${revel_db}' provided."
        }
        if (!file(gnomad_db).exists()) {
            exit 1, "gnomAD database file does not exist, gnomad_db = '${gnomad_db}' provided."
        }
        if (!file(clinvar_db).exists()) {
            exit 1, "ClinVar database file does not exist, clinvar_db = '${clinvar_db}' provided."
        }
        if (!file(cadd_snv_db).exists()) {
            exit 1, "CADD SNV database file does not exist, cadd_snv_db = '${cadd_snv_db}' provided."
        }
        if (!file(cadd_indel_db).exists()) {
            exit 1, "CADD indel database file does not exist, cadd_indel_db = '${cadd_indel_db}' provided."
        }
        if (!file(cadd_sv_db).exists()) {
            exit 1, "CADD SV database file does not exist, cadd_sv_db = '${cadd_sv_db}' provided."
        }
        if (!file(spliceai_snv_db).exists()) {
            exit 1, "SpliceAI SNV database file does not exist, spliceai_snv_db = '${spliceai_snv_db}' provided."
        }
        if (!file(spliceai_indel_db).exists()) {
            exit 1, "SpliceAI indel database file does not exist, spliceai_indel_db = '${spliceai_indel_db}' provided."
        }
        if (!file(alphamissense_db).exists()) {
            exit 1, "AlphaMissense database file does not exist, alphamissense_db = '${alphamissense_db}' provided."
        }
        if (!ref.toLowerCase().contains('hg38') && !ref.toLowerCase().contains('grch38') && annotate_override != 'yes') {
            exit 1, "Only hg38/GRCh38 is supported for annotation. It looks like you may not be passing a hg38/GRCh38 reference genome based on the filename of the reference genome. ref = '${ref}' provided. Pass '--annotate_override yes' on the command line to override this error."
        }
    }
    if (!(calculate_depth in ['yes', 'no'])) {
        exit 1, "Choice to calculate depth should be either 'yes', or 'no', calculate_depth = '${calculate_depth}' provided."
    }
    if (!(analyse_base_mods in ['yes', 'no'])) {
        exit 1, "Choice to analyse base modifications should be either 'yes', or 'no', analyse_base_mods = '${analyse_base_mods}' provided."
    }
    if (!(check_relatedness in ['yes', 'no'])) {
        exit 1, "Choice to check relatedness should be either 'yes', or 'no', check_relatedness = '${check_relatedness}' provided."
    }
    if (check_relatedness == 'yes') {
        if (!(mode in ['duo', 'trio'])) {
            exit 1, "Checking relatedness is only available when running in duo or trio mode, mode = '${mode}' and check_relatedness = '${check_relatedness}' provided."
        }
        if (sites == 'NONE') {
            exit 1, "When checking relatedness, set an appropriate sites file. sites = '${sites}' provided."
        }
        if (!file(sites).exists()) {
            exit 1, "Sites file does not exist, sites = '${sites}' provided."
        }
    }
    else if (check_relatedness == 'no') {
        if (sites != 'NONE') {
            exit 1, "When not checking relatedness, set sites file to 'NONE', sites = '${sites}' provided."
        }
    }
    if (!outdir) {
        exit 1, "No output directory (outdir) provided."
    }
    if (in_data_format in ['snp_indel_vcf', 'sv_vcf']) {
        if (mode != 'NONE') {
            exit 1, "When the in data format is SNP/indel VCF or SV VCF, set the mode to 'NONE', in_data_format = '${in_data_format}' and mode = '${mode}' provided."
        }
        if (tandem_repeat != 'NONE') {
            exit 1, "When the in data format is SNP/indel VCF or SV VCF, set the tandem repeat file to 'NONE', in_data_format = '${in_data_format}' and tandem_repeat = '${tandem_repeat}' provided."
        }
        if (annotate == 'no') {
            exit 1, "In data format is SNP/indel VCF or SV VCF, but you've chosen not to annotate, in_data_format = '${in_data_format}' and annotate = '${annotate}' provided. Nothing for pipeface to do."
        }
        if (calculate_depth == 'yes') {
            exit 1, "When the in data format is SNP/indel VCF or SV VCF, set calculate depth to 'no', in_data_format = '${in_data_format}' and calculate_depth = '${calculate_depth}' provided."
        }
        if (check_relatedness != 'no') {
            exit 1, "When the in data format is SNP/indel VCF or SV VCF, set checking relatedness to 'no', in_data_format = '${in_data_format}' and check_relatedness = '${check_relatedness}' provided."
        }
    }
    if (in_data_format == 'snp_indel_vcf') {
        if (sv_caller != 'NONE') {
            exit 1, "When the in data format is SNP/indel VCF, set the SV calling software to 'NONE', sv_caller = '${sv_caller}' provided."
        }
        if (snp_indel_caller == 'NONE') {
            exit 1, "When the input data format is 'snp_indel_vcf', pass the SNP/indel calling software which was used to generate the input data (not 'NONE'), snp_indel_caller = '${snp_indel_caller}' provided."
        }
        if (ref == 'NONE' || ref_index == 'NONE') {
            exit 1, "When the input data format is SNP/indel VCF, pass the reference genome and it's index which was used to generate the input data (not 'NONE'), ref = '${ref}' and ref_index = '${ref_index}' provided."
        }
    }
    else if (in_data_format == 'sv_vcf') {
        if (snp_indel_caller != 'NONE') {
            exit 1, "When the in data format is SV VCF, set the SNP/indel calling software to 'NONE', snp_indel_caller = '${snp_indel_caller}' provided."
        }
        if (sv_caller == 'NONE') {
            exit 1, "When the input data format is SV VCF, please pass the SV calling software which was used to generate the input data (not 'NONE'), sv_caller = '${sv_caller}' provided."
        }
    }

    // build variable
    ref_name = file(ref).getSimpleName()

    // build input tuple
    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            def sample_id = row.sample_id
            def family_id = row.family_id
            def files = row.file
            def data_type = row.data_type
            def regions_of_interest = row.regions_of_interest
            def clair3_model = row.clair3_model

            // only apply the haploidaware check here if it's on
            if (params.haploidaware == 'yes') {
                def check_file = (regions_of_interest != 'NONE' && file(regions_of_interest).exists())
                                    ? file(regions_of_interest)
                                    : file(params.ref_index)
                def fileContent = check_file.text
                def chrX_found = fileContent.contains(params.chrXseq)
                def chrY_found = fileContent.contains(params.chrYseq)
                if (!chrX_found || !chrY_found) {
                    exit 1, "Haploid-aware mode requires both chrX and chrY to be present in ${check_file}."
                }
            }
            return tuple(sample_id, family_id, files, data_type, regions_of_interest, clair3_model)
        }
        .groupTuple(by: [0,1,3,4,5])
        .set { in_data_tuple }

    // build channels from in_data.csv file for describing each metadata associated with samples/families
    // I can use these downstream to join to input tuples throughout the pipeline as needed to avoid needing to have all sample metadata in all process input/output tuples
    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple(row.sample_id, row.family_id) }
        .groupTuple(by: [0,1])
        .set { id_tuple }

    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple(row.sample_id, row.family_id, row.family_position) }
        .groupTuple(by: [0,1,2])
        .set { family_position_tuple }

    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple(row.sample_id, row.family_id, file(row.file).getExtension()) }
        .groupTuple(by: [0,1,2])
        .set { extension_tuple }

    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple(row.sample_id, row.family_id, row.file) }
        .groupTuple(by: [0,1])
        .set { files_tuple }

    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple(row.sample_id, row.family_id, "${row.file}.bai") }
        .groupTuple(by: [0,1])
        .set { index_tuple }

    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple(row.sample_id, row.family_id, row.data_type) }
        .groupTuple(by: [0,1,2])
        .set { data_type_tuple }

    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple(row.sample_id, row.family_id, row.regions_of_interest) }
        .groupTuple(by: [0,1,2])
        .set { regions_of_interest_tuple }

    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple(row.sample_id, row.family_id, row.clair3_model) }
        .groupTuple(by: [0,1,2])
        .set { clair3_model_tuple }

    // check user provided parameters in in_data.csv file
    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            def sample_id = row.sample_id
            def family_id = row.family_id
            def family_position = row.family_position
            def in_file = row.file
            def data_type = row.data_type
            def regions_of_interest = row.regions_of_interest
            def clair3_model = row.clair3_model

            if (sample_id.isEmpty()) {
                exit 1, "There is an empty entry in the 'sample_id' column of '${in_data}'."
            }
            if (in_file.isEmpty()) {
                exit 1, "There is an empty entry in the 'file' column of '${in_data}'."
            }
            if (!file(in_file).exists()) {
                exit 1, "There is an entry in the 'file' column of '${in_data}' which doesn't exist. Check file '${in_file}'."
            }
            if (!(file(in_file).getExtension() in ['bam', 'gz', 'fastq'])) {
                exit 1, "There is an entry in the 'file' column of '$in_data' which doesn't have a 'bam', 'gz' or 'fastq' file extension. Check file '${in_file}'."
            }
            if (in_data_format == 'ubam_fastq') {
                if (snp_indel_caller == 'clair3' && clair3_model == 'NONE') {
                    exit 1, "When clair3 is selected as the SNP/indel calling software, provide a path to an appropriate clair3 model in the 'clair3_model' column of '${in_data}' rather than setting it to 'NONE'."
                }
            }
            else if (in_data_format == 'aligned_bam') {
                if (!file(in_file + '.bai').exists()) {
                    exit 1, "You've specified that the in data format is aligned BAM, but it looks like '${in_file}' defined in the file column of '${in_data}' isn't indexed (ie. '${in_file}.bai doesn't exist')."
                }
                if (!in_file.contains('bam') && in_data_format_override != 'yes') {
                    exit 1, "You've specified that the in data format is aligned BAM, but it looks like you may not be passing a BAM file based on the file names defined in the file column of '${in_data}'. Check file '${in_file}'. Pass '--in_data_format_override yes' on the command line to override this error."
                }
            }
            else if (in_data_format in ['snp_indel_vcf', 'sv_vcf']) {
                if (data_type != 'NONE') {
                    exit 1, "When the input data format is SNP/indel VCF or SV VCF, please set the data type (data_type column of '${in_data}') to 'NONE'."
                }
                if (!in_file.contains('vcf') && in_data_format_override != 'yes') {
                    exit 1, "You've specified that the in data format is SNP/indel VCF or SV VCF, but it looks like you may not be passing a VCF file based on the file names defined in the file column of '${in_data}'. ${in_file} passed. Pass '--in_data_format_override yes' on the command line to override this error."
                }
            }
            if (data_type.isEmpty()) {
                exit 1, "There is an empty entry in the 'data_type' column of '${in_data}'."
            }
            if (in_data_format in ['ubam_fastq', 'aligned_bam']) {
                if (!(data_type in ['ont', 'pacbio'])) {
                    exit 1, "Entries in the 'data_type' column of '${in_data}' should be 'ont' or 'pacbio', '${data_type}' provided."
                }
            }
            if (regions_of_interest.isEmpty()) {
                exit 1, "There is an empty entry in the 'regions_of_interest' column of '${in_data}'."
            }
            if (!file(regions_of_interest).exists()) {
                exit 1, "There is an entry in the 'regions_of_interest' column of '${in_data}' which doesn't exist. Check file '${regions_of_interest}'."
            }
            if (clair3_model.isEmpty()) {
                exit 1, "There is an empty entry in the 'clair3_model' column of '${in_data}'."
            }
            if (!file(clair3_model).exists()) {
                exit 1, "There is an entry in the 'clair3_model' column of '${in_data}' which doesn't exist. Check path '${clair3_model}'."
            }
            if (clair3_model != 'NONE') {
                if (snp_indel_caller != 'clair3') {
                    exit 1, "Pass 'NONE' in the 'clair3_model' column of '${in_data}' when clair3 is NOT selected as the SNP/indel calling software, '${clair3_model}' provided'."
                }
                if (in_data_format in ['snp_indel_vcf', 'sv_vcf']) {
                    exit 1, "When the input data format is SNP/indel VCF or SV VCF, please set the Clair3 model (clair3_model column of '${in_data}') to 'NONE'."
                }
            }
        }

    // check user provided parameters relating in in_data.csv file relating to cohorts
    if (mode in ['duo', 'trio']) {
        Channel
            .fromPath(in_data)
            .splitCsv(header: true, sep: ',', strip: true)
            .map { row -> tuple(row.sample_id, row.family_id, row.family_position, row.file, row.data_type, row.regions_of_interest, row.clair3_model) }
            .groupTuple(by: 1)
            .map { sample_ids, family_ids, family_positions, files, data_types, regions_of_interests, clair3_models ->
                if (!family_positions.every {it in ['proband', 'father', 'mother']}) {
                    exit 1, "Entries in the 'family_position' column of '$in_data' should contain 'proband', 'father' and 'mother' for every 'family_id' in duo/trio mode, '$family_positions' provided for family '$family_ids'. 'mode = ${mode}' provided."
                }
                if (mode == "duo") {
                    if (sample_ids.unique().size() != 2) {
                        exit 1, "Entries in the 'sample_id' column of '$in_data' should contain 2 unique values for every 'family_id' in duo mode, '$sample_ids' provided for family '$family_ids'. 'mode = ${mode}' provided."
                    }
                }
                if (mode == "trio") {
                    if (sample_ids.unique().size() != 3) {
                        exit 1, "Entries in the 'sample_id' column of '$in_data' should contain 3 unique values for every 'family_id' in trio mode, '$sample_ids' provided for family '$family_ids'. 'mode = ${mode}' provided."
                    }
                }
        }
    }

    // workflow
    // pre-process, alignment and qc
    scrape_settings(in_data_tuple.join(family_position_tuple, by: [0,1]), pipeface_version, in_data, in_data_format, ref, ref_index, tandem_repeat, mode, snp_indel_caller, sv_caller, annotate, calculate_depth, analyse_base_mods, check_relatedness, sites, outdir, outdir2, haploidaware, sex, parbed)
    if (in_data_format == 'ubam_fastq') {
        merged = merge_runs(id_tuple.join(extension_tuple, by: [0,1]).join(files_tuple, by: [0,1]))
        bam = minimap2(merged.join(extension_tuple, by: [0,1]).join(data_type_tuple, by: [0,1]), ref, ref_index)
    }
    if (in_data_format == 'aligned_bam') {
        bam = id_tuple.join(files_tuple, by: [0,1]).join(index_tuple, by: [0,1])
    }
    if (in_data_format in ['ubam_fastq', 'aligned_bam']) {
        if (calculate_depth == 'yes') {
            mosdepth(bam.join(regions_of_interest_tuple, by: [0,1]), outdir, outdir2, ref_name)
        }
        // snp/indel calling
        if (snp_indel_caller == 'clair3') {
            if (haploidaware == 'no' | sex == 'XX') {
                (snp_indel_vcf_bam, gvcf) = clair3(bam.join(data_type_tuple, by: [0,1]).join(regions_of_interest_tuple, by: [0,1]).join(clair3_model_tuple, by: [0,1]), ref, ref_index, outdir, outdir2, ref_name, snp_indel_caller)
            }
            if (haploidaware == 'yes' && sex == 'XY') {
                bam_diploid_haploid_bed = clair3_pre_processing(bam.join(regions_of_interest_tuple, by: [0,1]), ref, ref_index, parbed)
                haploid_diploid_vcf_gvcf = clair3_haploid_aware(bam_diploid_haploid_bed.join(data_type_tuple, by: [0,1]).join(clair3_model_tuple, by: [0,1]), ref, ref_index, outdir, outdir2, ref_name)
                (snp_indel_vcf_bam, gvcf) = clair3_post_processing(haploid_diploid_vcf_gvcf, outdir, outdir2, ref_name, snp_indel_caller)
            }
        }
        else if (snp_indel_caller in ['deepvariant', 'deeptrio']) {
            dv_commands = deepvariant_dry_run(bam.join(data_type_tuple, by: [0,1]), ref, ref_index, sex, haploidaware, parbed)
            dv_examples = deepvariant_make_examples(dv_commands.join(regions_of_interest_tuple, by: [0,1]), ref, ref_index, parbed)
            dv_calls = deepvariant_call_variants(dv_examples)
            (snp_indel_raw_vcf_bam, snp_indel_gvcf_bam, gvcf) = deepvariant_post_processing(dv_calls.join(family_position_tuple, by: [0,1]), ref, ref_index, outdir, outdir2, ref_name, snp_indel_caller)
            // filter refcall variants
            snp_indel_vcf_bam = filter_ref_call(snp_indel_raw_vcf_bam)
        }
        // split multiallelic variants
        snp_indel_split_vcf_bam = split_multiallele(snp_indel_vcf_bam, ref, ref_index)
        // phasing
        (snp_indel_split_phased_vcf_bam, snp_indel_split_phased_vcf, phased_read_list) = whatshap_phase(snp_indel_split_vcf_bam, ref, ref_index, outdir, outdir2, ref_name, snp_indel_caller)
        // haplotagging
        (haplotagged_bam, haplotagged_bam_fam, haplotagged_tsv) = whatshap_haplotag(snp_indel_split_phased_vcf_bam.join(family_position_tuple, by: [0,1]), ref, ref_index, outdir, outdir2, ref_name)
        // base mod analysis
        if (analyse_base_mods == 'yes') {
            minimod(haplotagged_bam.join(data_type_tuple, by: [0,1]), ref, ref_index, outdir, outdir2, ref_name)
        }
        // joint snp/indel calling
        if (mode == 'duo') {
            tmp = snp_indel_gvcf_bam.groupTuple(by: 1).transpose()
            proband_gvcf_bam = tmp.filter { tuple -> tuple[2].contains("proband") }
            parent_gvcf_bam = tmp.filter { tuple -> tuple[2].contains("father") || tuple[2].contains("mother") }
            // check relatedness
            if (check_relatedness == 'yes') {
                somalier_duo(proband_gvcf_bam, parent_gvcf_bam, ref, ref_index, sites, outdir, outdir2, ref_name)
            }
            // gvcf merging
            joint_snp_indel_bcf_bam = glnexus_duo(proband_gvcf_bam, parent_gvcf_bam)
            joint_snp_indel_vcf_bam = glnexus_post_processing_duo(joint_snp_indel_bcf_bam)
            // joint split multiallelic variants
            joint_snp_indel_split_vcf_bam = split_multiallele_duo(joint_snp_indel_vcf_bam, ref, ref_index)
            // joint phasing
            (joint_snp_indel_split_phased_vcf, joint_phased_read_list) = whatshap_phase_duo(joint_snp_indel_split_vcf_bam, ref, ref_index, outdir, outdir2, ref_name, snp_indel_caller)
        }
        if (mode == 'trio') {
            tmp = haplotagged_bam_fam.groupTuple(by: 1).transpose()
            proband_bam = tmp.filter { tuple -> tuple[2].contains("proband") }
            father_bam = tmp.filter { tuple -> tuple[2].contains("father") }
            mother_bam = tmp.filter { tuple -> tuple[2].contains("mother") }
            dt_commands = deeptrio_dry_run(proband_bam.join(data_type_tuple, by: [0,1]), father_bam.join(data_type_tuple, by: [0,1]), mother_bam.join(data_type_tuple, by: [0,1]), ref, ref_index)
            dt_examples = deeptrio_make_examples(dt_commands, ref, ref_index)
            dt_calls = deeptrio_call_variants(dt_examples.proband.mix(dt_examples.father, dt_examples.mother))
            snp_indel_gvcf_bam = deeptrio_postprocessing(dt_calls, ref, ref_index)
            tmp = snp_indel_gvcf_bam.groupTuple(by: 1).transpose()
            proband_gvcf_bam = tmp.filter { tuple -> tuple[2].contains("proband") }
            father_gvcf_bam = tmp.filter { tuple -> tuple[2].contains("father") }
            mother_gvcf_bam = tmp.filter { tuple -> tuple[2].contains("mother") }
            // check relatedness
            if (check_relatedness == 'yes') {
                somalier_trio(proband_gvcf_bam, father_gvcf_bam, mother_gvcf_bam, ref, ref_index, sites, outdir, outdir2, ref_name)
            }
            // gvcf merging
            joint_snp_indel_bcf_bam = glnexus_trio(proband_gvcf_bam, father_gvcf_bam, mother_gvcf_bam)
            joint_snp_indel_vcf_bam = glnexus_post_processing_trio(joint_snp_indel_bcf_bam)
            // joint split multiallelic variants
            joint_snp_indel_split_vcf_bam = split_multiallele_trio(joint_snp_indel_vcf_bam, ref, ref_index)
            // joint phasing
            (joint_snp_indel_split_phased_vcf, joint_phased_read_list) = whatshap_phase_trio(joint_snp_indel_split_vcf_bam, ref, ref_index, outdir, outdir2, ref_name, snp_indel_caller)
        }
        if (mode in ['duo', 'trio']) {
            // joint snp/indel annotation
            if (annotate == 'yes') {
                vep_snp_indel(joint_snp_indel_split_phased_vcf, ref, ref_index, vep_db, revel_db, gnomad_db, clinvar_db, cadd_snv_db, cadd_indel_db, spliceai_snv_db, spliceai_indel_db, alphamissense_db, outdir, outdir2, ref_name, snp_indel_caller, mode)
            }
        }
        // sv calling
        if (sv_caller in ['sniffles', 'both']) {
            (sv_vcf_sniffles, sv_vcf_sniffles_indexed, sv_vcf_haplotagged_bam_fam_sniffles) = sniffles(haplotagged_bam.join(family_position_tuple, by: [0,1]), ref, ref_index, tandem_repeat, outdir, outdir2, ref_name)
        }
        if (sv_caller in ['cutesv', 'both']) {
            (sv_vcf_cutesv, sv_vcf_cutesv_indexed, sv_vcf_haplotagged_bam_fam_cutesv) = cutesv(haplotagged_bam.join(data_type_tuple, by: [0,1]).join(family_position_tuple, by: [0,1]), ref, ref_index, tandem_repeat, outdir, outdir2, ref_name)
        }
    }
    if (in_data_format == 'snp_indel_vcf') {
        snp_indel_split_phased_vcf = id_tuple.join(files_tuple, by: [0,1])
    }
    if (in_data_format in ['ubam_fastq', 'aligned_bam', 'snp_indel_vcf']) {
        // annotation
        if (annotate == 'yes' && !(mode in ['duo', 'trio'])) {
            vep_snp_indel(snp_indel_split_phased_vcf, ref, ref_index, vep_db, revel_db, gnomad_db, clinvar_db, cadd_snv_db, cadd_indel_db, spliceai_snv_db, spliceai_indel_db, alphamissense_db, outdir, outdir2, ref_name, snp_indel_caller, mode)
        }
    }
    if (mode == 'duo') {
        // sv vcf merging
        if (sv_caller in ['sniffles', 'both']) {
            sniffles_tmp = sv_vcf_haplotagged_bam_fam_sniffles.groupTuple(by: 1).transpose()
            sniffles_proband_tuple = sniffles_tmp.filter { tuple -> tuple[2].contains("proband") }
            sniffles_parent_tuple = sniffles_tmp.filter { tuple -> tuple[2].contains("father") || tuple[2].contains("mother") }
            (joint_sv_vcf_sniffles, joint_sv_vcf_sniffles_indexed) = jasmine_sniffles_duo(sniffles_proband_tuple.join(data_type_tuple, by: [0,1]), sniffles_parent_tuple, ref, ref_index, outdir, outdir2, ref_name)
        }
        if (sv_caller in ['cutesv', 'both']) {
            cutesv_tmp = sv_vcf_haplotagged_bam_fam_cutesv.groupTuple(by: 1).transpose()
            cutesv_proband_tuple = cutesv_tmp.filter { tuple -> tuple[2].contains("proband") }
            cutesv_parent_tuple = cutesv_tmp.filter { tuple -> tuple[2].contains("father") || tuple[2].contains("mother") }
            (joint_sv_vcf_cutesv, joint_sv_vcf_cutesv_indexed) = jasmine_cutesv_duo(cutesv_proband_tuple.join(data_type_tuple, by: [0,1]), cutesv_parent_tuple, ref, ref_index, outdir, outdir2, ref_name)
        }
    }
    if (mode == 'trio') {
        // sv vcf merging
        if (sv_caller in ['sniffles', 'both']) {
            sniffles_tmp = sv_vcf_haplotagged_bam_fam_sniffles.groupTuple(by: 1).transpose()
            sniffles_proband_tuple = sniffles_tmp.filter { tuple -> tuple[2].contains("proband") }
            sniffles_father_tuple = sniffles_tmp.filter { tuple -> tuple[2].contains("father") }
            sniffles_mother_tuple = sniffles_tmp.filter { tuple -> tuple[2].contains("mother") }
            (joint_sv_vcf_sniffles, joint_sv_vcf_sniffles_indexed) = jasmine_sniffles_trio(sniffles_proband_tuple.join(data_type_tuple, by: [0,1]), sniffles_father_tuple, sniffles_mother_tuple, ref, ref_index, outdir, outdir2, ref_name)
        }
        if (sv_caller in ['cutesv', 'both']) {
            cutesv_tmp = sv_vcf_haplotagged_bam_fam_cutesv.groupTuple(by: 1).transpose()
            cutesv_proband_tuple = cutesv_tmp.filter { tuple -> tuple[2].contains("proband") }
            cutesv_father_tuple = cutesv_tmp.filter { tuple -> tuple[2].contains("father") }
            cutesv_mother_tuple = cutesv_tmp.filter { tuple -> tuple[2].contains("mother") }
            (joint_sv_vcf_cutesv, joint_sv_vcf_cutesv_indexed) = jasmine_cutesv_trio(cutesv_proband_tuple.join(data_type_tuple, by: [0,1]), cutesv_father_tuple, cutesv_mother_tuple, ref, ref_index, outdir, outdir2, ref_name)
        }
    }
    if (mode in ['duo', 'trio']) {
        // joint sv annotation
        if (annotate == 'yes') {
            if (sv_caller in ['sniffles', 'both']) {
                vep_sniffles_sv(joint_sv_vcf_sniffles, ref, ref_index, vep_db, gnomad_db, cadd_sv_db, outdir, outdir2, ref_name, mode)
            }
            if (sv_caller in ['cutesv', 'both']) {
                vep_cutesv_sv(joint_sv_vcf_cutesv, ref, ref_index, vep_db, gnomad_db, cadd_sv_db, outdir, outdir2, ref_name, mode)
            }
        }
    }
    if (in_data_format == 'sv_vcf' && !(mode in ['duo', 'trio'])) {
        sv_vcf_sniffles = id_tuple.join(files_tuple, by: [0,1])
        sv_vcf_cutesv = id_tuple.join(files_tuple, by: [0,1])
    }
    if (in_data_format in ['ubam_fastq', 'aligned_bam', 'sv_vcf']) {
        // annotation
        if (annotate == 'yes' && !(mode in ['duo', 'trio'])) {
            if (sv_caller in ['sniffles', 'both']) {
                vep_sniffles_sv(sv_vcf_sniffles, ref, ref_index, vep_db, gnomad_db, cadd_sv_db, outdir, outdir2, ref_name, mode)
            }
            if (sv_caller in ['cutesv', 'both']) {
                vep_cutesv_sv(sv_vcf_cutesv, ref, ref_index, vep_db, gnomad_db, cadd_sv_db, outdir, outdir2, ref_name, mode)
            }
        }
    }
}
