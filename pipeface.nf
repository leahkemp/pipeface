nextflow.enable.dsl=2

// tag pipeface version
def pipeface_version = "0.11.0"

// set defaults for undocumented params
params.outdir2 = ""
params.annotate_override = ""
params.in_data_format_override = ""

// default values for chrX and chrY contig names. Used when sex = 'XY' and 
// haploidaware = 'yes'. Don't change unless your reference genome has different names. 
params.chrXseq = "chrX"
params.chrYseq = "chrY"

process scrape_settings {

    publishDir "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$sample_id.$filename" }, pattern: '*pipeface_settings.txt'

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
        val tr_calling
        val tr_call_regions
        val check_relatedness
        val sites
        val somatic_calling
        val outdir
        val outdir2
        val haploidaware
        val sex
        val parbed
        val sv_mapq

    output:
        tuple val(sample_id), val(family_id), path("pipeface_settings.txt")

    script:
        // conditionally define reported sv caller and in data format
        reported_sv_caller = [both: 'cutesv,sniffles', sniffles: 'sniffles', cutesv: 'cutesv', NONE: 'NONE'][sv_caller]
        reported_in_data_format = [ubam_fastq: 'unaligned BAM or FASTQ', aligned_bam: 'aligned BAM', snp_indel_vcf: 'SNP/indel VCF', sv_vcf: 'SV VCF'][in_data_format]
        if (in_data_format in ['ubam_fastq', 'aligned_bam'])
        """
        {
            echo "Pipeface version: $pipeface_version"
            echo "Sample ID: $sample_id"
            echo "Family ID: $family_id"
            echo "Family position: $family_position"
            echo "In data format: $reported_in_data_format"
            echo "Input data file/files: $files"
            echo "Data type: $data_type"
            echo "Regions of interest file: $regions_of_interest"
            echo "Clair3 model: $clair3_model"
            echo "In data csv path: $in_data"
            echo "Reference genome: $ref"
            echo "Reference genome index: $ref_index"
            echo "Haploid-aware: $haploidaware"
            echo "Sex: $sex"
            echo "PAR regions: $parbed"
            echo "Tandem repeat file: $tandem_repeat"
            echo "Mode: $mode"
            echo "SNP/indel caller: $snp_indel_caller"
            echo "SV caller: $reported_sv_caller"
            echo "SV MAPQ filter threshold: $sv_mapq"
            echo "Annotate: $annotate"
            echo "Calculate depth: $calculate_depth"
            echo "Analyse base modifications: $analyse_base_mods"
            echo "Tandem repeat calling: $tr_calling"
            echo "Tandem repeat call regions: $tr_call_regions"
            echo "Check relatedness: $check_relatedness"
            echo "Somalier sites file: $sites"
            echo "Somatic calling: $somatic_calling"
            echo "Outdir: $outdir"
        } > pipeface_settings.txt
        """
        else if (in_data_format in ['snp_indel_vcf', 'sv_vcf'])
        """
        {
            echo "Pipeface version: $pipeface_version"
            echo "Sample ID: $sample_id"
            echo "Family ID: $family_id"
            echo "Family position: $family_position"
            echo "In data format: $reported_in_data_format"
            echo "Input data file/files: $files"
            echo "In data csv path: $in_data"
            echo "Annotate: $annotate"
            echo "Outdir: $outdir"
        } > pipeface_settings.txt
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
        tuple val(sample_id), val(family_id), path("merged")

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
        tuple val(sample_id), val(family_id), path("sorted.bam"), path("sorted.bam.bai")

    script:
        // conditionally define preset
        def preset = data_type == 'ont' ? 'lr:hq' : 'map-hifi'
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

    publishDir "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$sample_id.${ref_name}.mosdepth.$filename" }, pattern: 'depth.txt'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(regions_of_interest)
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path("depth.txt")

    script:
        // optionally pass regions of interest bed file
        def regions_of_interest_optional = regions_of_interest != 'NONE' ? "-b $regions_of_interest" : ''
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

process clair3_pre_processing {

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(regions_of_interest)
        val ref
        val ref_index
        path parbed

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path("diploid.bed"), path("haploid.bed")

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

process clair3 {

    publishDir "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$snp_indel_caller.$filename" }, pattern: 'snp_indel.g.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(data_type), val(regions_of_interest), val(clair3_model)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path("snp_indel.vcf.gz"), path("snp_indel.vcf.gz.tbi")
        tuple val(sample_id), val(family_id), path("snp_indel.g.vcf.gz"), path("snp_indel.g.vcf.gz.tbi")

    script:
        // optionally pass regions of interest bed file
        def regions_of_interest_optional = regions_of_interest != 'NONE' ? "--bed_fn=$regions_of_interest" : ''
        // conditionally define platform
        def platform = data_type == 'ont' ? 'ont' : 'hifi'
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

process clair3_haploid_aware {

    publishDir "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$snp_indel_caller.$filename" }, pattern: 'snp_indel.g.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path(diploid_bed), path(haploid_bed), val(data_type), val(clair3_model)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path("haploid_snp_indel.vcf.gz"), path("haploid_snp_indel.vcf.gz.tbi"), path("diploid_snp_indel.vcf.gz"), path("diploid_snp_indel.vcf.gz.tbi"), path("haploid_snp_indel.g.vcf.gz"), path("haploid_snp_indel.g.vcf.gz.tbi"), path("diploid_snp_indel.g.vcf.gz"), path("diploid_snp_indel.g.vcf.gz.tbi")

    script:
        // conditionally define platform
        def platform = data_type == 'ont' ? 'ont' : 'hifi'
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

    publishDir "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$snp_indel_caller.$filename" }, pattern: 'snp_indel.g.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path(haploid_snp_indel_vcf), path(haploid_snp_indel_vcf_index), path(diploid_snp_indel_vcf), path(diploid_snp_indel_vcf_index), path(haploid_snp_indel_gvcf), path(haploid_snp_indel_gvcf_index), path(diploid_snp_indel_gvcf), path(diploid_snp_indel_gvcf_index)
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path("snp_indel.vcf.gz"), path("snp_indel.vcf.gz.tbi")
        tuple val(sample_id), val(family_id), path("snp_indel.g.vcf.gz"), path("snp_indel.g.vcf.gz.tbi")

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

process clairs_to {

    publishDir "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$sample_id.${ref_name}.clairs_to.$filename" }, pattern: '{snv,indel}.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(clairs_to_platform), val(regions_of_interest)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path("snv.vcf.gz"), path("snv.vcf.gz.tbi"), path("indel.vcf.gz"), path("indel.vcf.gz.tbi")

    script:
        // optionally pass regions of interest bed file
        def regions_of_interest_optional = regions_of_interest != 'NONE' ? "--bed_fn $regions_of_interest" : ''
        """
        run_clairs_to --tumor_bam_fn $bam --ref_fn $ref --output_dir ./ --threads ${task.cpus} --platform $clairs_to_platform --sample_name $sample_id --include_all_ctgs $regions_of_interest_optional
        """

    stub:
        """
        touch snv.vcf.gz
        touch snv.vcf.gz.tbi
        touch indel.vcf.gz
        touch indel.vcf.gz.tbi
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
        def model = data_type == 'ont' ? 'ONT_R104' : 'PACBIO'
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
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(call_variants_args), path("make_examples.*.gz{,.example_info.json}"), path("make_examples_call_variant_outputs.*.gz"), path("gvcf.*.gz")

    script:
        // conditionally pass regions of interest bed file
        def regions_of_interest_optional = regions_of_interest != 'NONE' ? "--regions $regions_of_interest" : ''
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
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(make_examples_call_variant_out), val(gvcf), path("call_variants_output*.gz")

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

    publishDir "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id", mode: params.publish_mode, overwrite: true, saveAs: { filename ->
        // when deeptrio is selected, files are still published under the 'deepvariant' caller name
        def caller = params.snp_indel_caller != 'deeptrio' ? snp_indel_caller : 'deepvariant'
        return "$sample_id.$ref_name.$caller.$filename"
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
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path("snp_indel.vcf.gz"), path("snp_indel.vcf.gz.tbi")
        tuple val(sample_id), val(family_id), val(family_position), path("${family_position}.sorted.bam"), path("${family_position}.sorted.bam.bai"), path("${family_position}_snp_indel.g.vcf.gz")
        tuple val(sample_id), path("snp_indel.g.vcf.gz"), path("snp_indel.g.vcf.gz.tbi")

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
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path("snp_indel.filtered.vcf.gz"), path("snp_indel.filtered.vcf.gz.tbi")

    script:
        """
        # filter out refcall variants
        bcftools view -f 'PASS' $snp_indel_vcf -o snp_indel.filtered.vcf.gz
        # index vcf
        tabix snp_indel.filtered.vcf.gz
        """
    stub:
        """
        touch snp_indel.filtered.vcf.gz
        touch snp_indel.filtered.vcf.gz.tbi
        """

}

process split_multiallele {

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path(snp_indel_vcf), path(snp_indel_vcf_index)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path("snp_indel.split.vcf.gz"), path("snp_indel.split.vcf.gz.tbi")

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

    publishDir "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id", mode: params.publish_mode, overwrite: true, saveAs: { filename ->
        // when deeptrio is selected, files are still published under the 'deepvariant' caller name
        def caller = params.snp_indel_caller != 'deeptrio' ? snp_indel_caller : 'deepvariant'
        return "$sample_id.$ref_name.$caller.$filename"
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
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path("snp_indel.phased.vcf.gz"), path("snp_indel.phased.vcf.gz.tbi")
        tuple val(sample_id), val(family_id), path("snp_indel.phased.vcf.gz")
        tuple val(sample_id), val(family_id), path("snp_indel.phased.vcf.gz"), path("snp_indel.phased.vcf.gz.tbi"), path("snp_indel.phased.read_list.txt"), path("snp_indel.phased.stats.gtf")

    script:
        """
        # run whatshap phase
        whatshap phase --reference $ref --output snp_indel.phased.vcf.gz --output-read-list snp_indel.phased.read_list.txt --sample $sample_id --ignore-read-groups $snp_indel_split_vcf $bam
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

    publishDir "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$sample_id.${ref_name}.minimap2.whatshap.$filename" }, pattern: 'sorted.haplotagged.*'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path(snp_indel_split_vcf), path(snp_indel_split_vcf_index), val(family_position)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path("sorted.haplotagged.bam"), path("sorted.haplotagged.bam.bai")
        tuple val(sample_id), val(family_id), val(family_position), path("${family_position}.sorted.haplotagged.bam"), path("${family_position}.sorted.haplotagged.bam.bai")
        tuple val(sample_id), val(family_id), path("sorted.haplotagged.tsv")

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
        def model = proband_data_type == 'ont' ? 'ONT' : 'PACBIO'
        """
        run_deeptrio --model_type=$model --ref=$ref \
            --sample_name_child=$proband_sample_id --sample_name_parent1=$father_sample_id --sample_name_parent2=$mother_sample_id \
            --reads_child=$proband_haplotagged_bam --reads_parent1=$father_haplotagged_bam --reads_parent2=$mother_haplotagged_bam \
            --output_vcf_child=child.vcf.gz --output_vcf_parent1=parent1.vcf.gz --output_vcf_parent2=parent2.vcf.gz \
            --output_gvcf_child=child.g.vcf.gz --output_gvcf_parent1=parent1.g.vcf.gz --output_gvcf_parent2=parent2.g.vcf.gz \
            --dry_run=true > commands.txt
        make_examples_cs_args=\$(grep "/opt/deepvariant/bin/deeptrio/make_examples --mode candidate_sweep" commands.txt | awk -F'/opt/deepvariant/bin/deeptrio/make_examples' '{print \$2}' \
            | sed 's/--ref "[^"]*"//g' | sed 's/--sample_name "[^"]*"//g' | sed 's/--reads "[^"]*"//g' | sed 's/--sample_name_parent1 "[^"]*"//g' | sed 's/--reads_parent1 "[^"]*"//g' \
            | sed 's/--sample_name_parent2 "[^"]*"//g' | sed 's/--reads_parent2 "[^"]*"//g' | sed 's/--examples "[^"]*"//g' | sed 's/--candidate_positions "[^"]*"//g' | sed 's/--gvcf "[^"]*"//g')
        make_examples_calling_args=\$(grep "/opt/deepvariant/bin/deeptrio/make_examples --mode calling" commands.txt | awk -F'/opt/deepvariant/bin/deeptrio/make_examples' '{print \$2}' \
            | sed 's/--ref "[^"]*"//g' | sed 's/--sample_name "[^"]*"//g' | sed 's/--reads "[^"]*"//g' | sed 's/--sample_name_parent1 "[^"]*"//g' | sed 's/--reads_parent1 "[^"]*"//g' \
            | sed 's/--sample_name_parent2 "[^"]*"//g' | sed 's/--reads_parent2 "[^"]*"//g' | sed 's/--examples "[^"]*"//g' | sed 's/--candidate_positions "[^"]*"//g' | sed 's/--gvcf "[^"]*"//g')
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
        tuple val(proband_sample_id), val(proband_family_id), val(proband_family_position), path(proband_haplotagged_bam), path(proband_haplotagged_bam_index), path("make_examples_child.*.gz"), path("gvcf_child.*.gz"), path("*.example_info.json"), val(call_variants_proband_args), emit: proband
        tuple val(father_sample_id), val(father_family_id), val(father_family_position), path(father_haplotagged_bam), path(father_haplotagged_bam_index), path("make_examples_parent1.*.gz"), path("gvcf_parent1.*.gz"), path("*.example_info.json"), val(call_variants_father_args), emit: father
        tuple val(mother_sample_id), val(mother_family_id), val(mother_family_position), path(mother_haplotagged_bam), path(mother_haplotagged_bam_index), path("make_examples_parent2.*.gz"), path("gvcf_parent2.*.gz"), path("*.example_info.json"), val(call_variants_mother_args), emit: mother

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
        tuple val(sample_id), val(family_id), val(family_position), path(haplotagged_bam), path(haplotagged_bam_index), path(make_examples), val(gvcf), path(example_info), val(call_variants_args)

    output:
        tuple val(sample_id), val(family_id), val(family_position), path(haplotagged_bam), path(haplotagged_bam_index), path("*.gz"), val(gvcf)

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
        tuple val(sample_id), val(family_id), val(family_position), path(haplotagged_bam), path(haplotagged_bam_index), path(call_variants), path(gvcf)
        val ref
        val ref_index
        
    output:
        tuple val(sample_id), val(family_id), val(family_position), path(haplotagged_bam), path(haplotagged_bam_index), path("${family_position}_snp_indel.g.vcf.gz")

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

process somalier_extract {

    publishDir "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id", mode: params.publish_mode, overwrite: true, saveAs: { filename ->
        def clean = filename.replaceFirst(/^.*\.somalier$/, 'somalier')
        return "$sample_id.$ref_name.$clean"
    }, pattern: '*.somalier'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index)
        val ref
        val ref_index
        val sites
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path("${sample_id}.somalier")

    script:
        """
        SOMALIER_SAMPLE_NAME=$sample_id somalier extract -d extracted --sites $sites -f $ref $bam
        cp extracted/* .
        """

    stub:
        """
        touch ${sample_id}.somalier
        """

}

process somalier_relate {

    publishDir "$outdir/$family_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$family_id.$ref_name.$filename" }, pattern: 'somalier.*'

    input:
        tuple val(proband_sample_id), val(family_id), path(somalier_files)
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(proband_sample_id), path("somalier.samples.tsv"), path("somalier.pairs.tsv"), path("somalier.html")

    script:
        """
        somalier relate $somalier_files
        """

    stub:
        """
        touch somalier.samples.tsv
        touch somalier.pairs.tsv
        touch somalier.html
        """

}

process glnexus_pre_processing {

    input:
        tuple val(sample_id), val(family_id), path(gvcf)
        val ref_index

    output:
        tuple val(sample_id), val(family_id), path("*.amended.g.vcf.gz")

    script:
        """
        # reheader to include all contigs in gvcf header even if no variants were called in a contig to avoid this issue: https://github.com/HKU-BAL/Clair3/issues/371
        {
            bcftools view -h "$gvcf" | grep -v -E '^##contig=|^#CHROM'
            awk '{ printf "##contig=<ID=%s,length=%s>\\n", \$1, \$2 }' "$ref_index"
            bcftools view -h "$gvcf" | grep '#CHROM'
        } > ${sample_id}.amended.g.vcf
        # convert lower cases of soft-masked sequences to upper case to avoid this issue: https://github.com/HKU-BAL/Clair3/issues/359
        zgrep -v '#' $gvcf | awk -F'\\t' -v OFS='\\t' '{ if(\$0 !~ /^#/) { \$4=toupper(\$4); \$5=toupper(\$5) } print }' >> ${sample_id}.amended.g.vcf
        # compress vcf
        bgzip -@ ${task.cpus} ${sample_id}.amended.g.vcf
        """

    stub:
        """
        touch ${sample_id}.amended.g.vcf.gz
        """

}

process glnexus {

    input:
        tuple val(proband_sample_id), val(family_id), val(sample_ids), val(family_positions), path(bams, stageAs: 'bam?.bam'), path(bam_indices, stageAs: 'bam?.bam.bai'), path(gvcfs)
        val snp_indel_caller
        val clair3_config

    output:
        tuple val(proband_sample_id), val(family_id), val(sample_ids), val(family_positions), path(bams), path(bam_indices), path("snp_indel.bcf")

    script:
        if (snp_indel_caller == 'clair3')
        """
        glnexus_cli --config $clair3_config $gvcfs > snp_indel.bcf
        """
        else if (snp_indel_caller in ['deepvariant', 'deeptrio'])
        """
        glnexus_cli --config DeepVariant $gvcfs > snp_indel.bcf
        """

    stub:
        """
        touch snp_indel.bcf
        """

}

process glnexus_post_processing {

    input:
        tuple val(proband_sample_id), val(family_id), val(sample_ids), val(family_positions), path(bams, stageAs: 'bam?.bam'), path(bam_indices, stageAs: 'bam?.bam.bai'), path(snp_indel_bcf)

    output:
        tuple val(proband_sample_id), val(family_id), val(sample_ids), val(family_positions), path(bams), path(bam_indices), path("snp_indel.vcf.gz"), path("snp_indel.vcf.gz.tbi")

    script:
        """
        # convert bcf to vcf, reorder since glnexus sorts samples alphabetically
        bcftools view $snp_indel_bcf | bcftools view -s ${sample_ids.join(',')} | bgzip -@ ${task.cpus} -c > snp_indel.vcf.gz
        tabix snp_indel.vcf.gz
        """

    stub:
        """
        touch snp_indel.vcf.gz
        touch snp_indel.vcf.gz.tbi
        """

}

process split_multiallele_family {

    input:
        tuple val(proband_sample_id), val(family_id), val(sample_ids), val(family_positions), path(bams, stageAs: 'bam?.bam'), path(bam_indices, stageAs: 'bam?.bam.bai'), path(snp_indel_vcf), path(snp_indel_vcf_index)
        val ref
        val ref_index

    output:
        tuple val(proband_sample_id), val(family_id), val(sample_ids), val(family_positions), path(bams), path(bam_indices), path("snp_indel.split.vcf.gz"), path("snp_indel.split.vcf.gz.tbi")

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

process whatshap_phase_family {

    publishDir "$outdir/$family_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$family_id.$ref_name.$snp_indel_caller.$filename" }, pattern: 'snp_indel.phased.*'

    input:
        tuple val(proband_sample_id), val(family_id), val(sample_ids), val(family_positions), path(bams, stageAs: 'bam?.bam'), path(bam_indices, stageAs: 'bam?.bam.bai'), path(snp_indel_split_vcf), path(snp_indel_split_vcf_index)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller

    output:
        tuple val(proband_sample_id), val(family_id), path("snp_indel.phased.vcf.gz")
        tuple val(family_id), path("snp_indel.phased.vcf.gz"), path("snp_indel.phased.vcf.gz.tbi"), path("snp_indel.phased.read_list.txt"), path("snp_indel.phased.stats.gtf")

    script:
        """
        # build pedigree, find father and mother by family_position (use '0' if absent for duo)
        FATHER_ID="0"
        MOTHER_ID="0"
        SAMPLE_IDS=(${sample_ids.join(' ')})
        FAMILY_POSITIONS=(${family_positions.join(' ')})
        for i in \${!SAMPLE_IDS[@]}; do
            [[ \${FAMILY_POSITIONS[\$i]} == "father" ]] && FATHER_ID=\${SAMPLE_IDS[\$i]}
            [[ \${FAMILY_POSITIONS[\$i]} == "mother" ]] && MOTHER_ID=\${SAMPLE_IDS[\$i]}
        done
        printf "$family_id\t$proband_sample_id\t\${FATHER_ID}\t\${MOTHER_ID}\t0\t1\n" > pedigree.ped
        # run whatshap phase
        whatshap phase --reference $ref --output snp_indel.phased.vcf.gz --output-read-list snp_indel.phased.read_list.txt --ped pedigree.ped $snp_indel_split_vcf ${bams.join(' ')}
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
            return "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id"
        } else {
            return "$outdir/$family_id/$outdir2"
        }
    }, mode: params.publish_mode, overwrite: true, saveAs: { filename ->
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
        tuple val(sample_id), val(family_id), path("snp_indel.phased.annotated.vcf.gz"), path("snp_indel.phased.annotated.vcf.gz.tbi")

    script:
        """
        # run vep
        vep -i $snp_indel_split_phased_vcf -o snp_indel.phased.annotated.vcf.gz --format vcf --vcf --fasta $ref --dir $vep_db --assembly GRCh38 --species homo_sapiens --cache --offline --merged \
        --sift b --polyphen b --symbol --hgvs --hgvsg --uploaded_allele --check_existing --filter_common --distance 0 --nearest gene --canonical --mane --pick \
        --fork ${task.cpus} --no_stats --compress_output bgzip --dont_skip \
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

    publishDir "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$sample_id.${ref_name}.minimod.$filename" }, pattern: 'modfreqs_*.b*'

    input:
        tuple val(sample_id), val(family_id), path(haplotagged_bam), path(haplotagged_bam_index), val(data_type)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path("modfreqs_hap1.bw"), path("modfreqs_hap1.bed"), path("modfreqs_hap2.bw"), path("modfreqs_hap2.bed"), path("modfreqs_combined.bw"), path("modfreqs_combined.bed"), path("modfreqs_unphased.bed"), optional: true

    script:
        """
        # run minimod
        minimod mod-freq $ref $haplotagged_bam -t ${task.cpus} --haplotypes -o modfreqs.tmp.bed
        # sort
        awk 'NR > 1 { print }' modfreqs.tmp.bed | sort -k1,1 -k2,2n > modfreqs.bed
        if [ -s modfreqs.bed ]; then
            # separate haplotypes
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

process longtr_pre_processing {

    input:
        val tr_call_regions

    output:
        path("split.*.bed")

    script:
        """
        # split up bed
        split -l 10000 $tr_call_regions split. --additional-suffix=.bed
        """

    stub:
        """
        touch split.aa.bed
        """

}

process longtr {

    input:
        tuple val(sample_id), val(family_id), path(haplotagged_bam), path(haplotagged_bam_index), val(data_type)
        path split_beds
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(family_id), path("tr.*.vcf.gz"), path("tr.*.vcf.gz.tbi")

    script:
        // define alignment parameters for ONT to account for higher incidence of indels in homopolymers (defaults are tailored to pacbio hifi)
        def alignment_params_optional = data_type == 'ont' ? "--alignment-params -1.0,-0.458675,-1.0,-0.458675,-0.00005800168,-1,-1" : ''
        """
        # run longtr and index vcfs in parallel for each split bed
        parallel -j ${task.cpus} '
            LongTR --bams $haplotagged_bam --bam-samps $sample_id --bam-libs $sample_id --fasta $ref --regions {} --tr-vcf tr.{/.}.vcf.gz --phased-bam --output-gls --output-pls --output-phased-gls --output-filter $alignment_params_optional --log longtr.{/.}.log
            tabix tr.{/.}.vcf.gz
        ' ::: split.*.bed
        """

    stub:
        """
        touch tr.aa.vcf.gz
        touch tr.aa.vcf.gz.tbi
        """

}

process longtr_family {

    input:
        tuple val(family_id), val(sample_ids), path(bams, stageAs: 'bam?.bam'), path(bam_indices, stageAs: 'bam?.bam.bai'), val(data_type)
        path split_beds
        val ref
        val ref_index

    output:
        tuple val(family_id), val(sample_ids), path("tr.*.vcf.gz"), path("tr.*.vcf.gz.tbi")

    script:
        def bams_csv = bams.join(',')
        def samples_csv = sample_ids.join(',')
        // define alignment parameters for ONT to account for higher incidence of indels in homopolymers (defaults are tailored to pacbio hifi)
        def alignment_params_optional = data_type == 'ont' ? "--alignment-params -1.0,-0.458675,-1.0,-0.458675,-0.00005800168,-1,-1" : ''
        """
        # run longtr and index vcfs in parallel for each split bed
        parallel -j ${task.cpus} '
            LongTR --bams $bams_csv --bam-samps $samples_csv --bam-libs $samples_csv --fasta $ref --regions {} --tr-vcf tr.{/.}.vcf.gz --phased-bam --output-gls --output-pls --output-phased-gls --output-filter $alignment_params_optional --log longtr.{/.}.log
            tabix tr.{/.}.vcf.gz
        ' ::: split.*.bed
        """

    stub:
        """
        touch tr.aa.vcf.gz
        touch tr.aa.vcf.gz.tbi
        """

}

process concat_tr_vcf {

    publishDir "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "${sample_id}.${ref_name}.longtr.$filename" }, pattern: 'tr.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(tr_vcfs), path(tr_vcf_indices)
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path("tr.vcf.gz"), path("tr.vcf.gz.tbi")

    script:
        """
        # get list of vcfs
        VCFS=(tr.*.vcf.gz)
        # concat vcfs (or use single vcf), then naturally sort variants
        if [[ \${#VCFS[@]} -eq 1 ]]; then
            bcftools sort -T ./ \${VCFS[0]} -Oz -o tr.vcf.gz
        else
            bcftools concat -a \${VCFS[@]} --threads ${task.cpus} | bcftools sort -T ./ -Oz -o tr.vcf.gz
        fi
        # index vcf
        tabix tr.vcf.gz
        """

    stub:
        """
        touch tr.vcf.gz
        touch tr.vcf.gz.tbi
        """

}

process concat_tr_vcf_family {

    publishDir "$outdir/$family_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "${family_id}.${ref_name}.longtr.$filename" }, pattern: 'tr.vcf.gz*'

    input:
        tuple val(family_id), val(sample_ids), path(tr_vcfs), path(tr_vcf_indices)
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(family_id), path("tr.vcf.gz"), path("tr.vcf.gz.tbi")

    script:
        """
        # get list of vcfs
        VCFS=(tr.*.vcf.gz)
        # concat vcfs (or use single vcf), reorder samples since longtr sorts samples alphabetically, then naturally sort variants
        if [[ \${#VCFS[@]} -eq 1 ]]; then
            bcftools view -s ${sample_ids.join(',')} \${VCFS[0]} | bcftools sort -T ./ -Oz -o tr.vcf.gz
        else
            bcftools concat -a \${VCFS[@]} --threads ${task.cpus} | bcftools view -s ${sample_ids.join(',')} | bcftools sort -T ./ -Oz -o tr.vcf.gz
        fi
        # index vcf
        tabix tr.vcf.gz
        """

    stub:
        """
        touch tr.vcf.gz
        touch tr.vcf.gz.tbi
        """

}

process sniffles {

    publishDir "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$sample_id.${ref_name}.sniffles.$filename" }, pattern: 'sv.phased.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(haplotagged_bam), path(haplotagged_bam_index), val(family_position)
        val ref
        val ref_index
        val tandem_repeat
        val outdir
        val outdir2
        val ref_name
        val sv_mapq

    output:
        tuple val(sample_id), val(family_id), path("sv.phased.vcf.gz")
        tuple val(sample_id), val(family_id), path("sv.phased.vcf.gz"), path("sv.phased.vcf.gz.tbi")
        tuple val(sample_id), val(family_id), val(family_position), path("sv.phased.vcf.gz"), path(haplotagged_bam)

    script:
        // optionally pass tandem repeat bed file
        def tandem_repeat_optional = tandem_repeat != 'NONE' ? "--tandem-repeats $tandem_repeat" : ''
        // optionally set mapq filter threshold
        def mapq_optional = sv_mapq != 'NONE' ? "--mapq ${sv_mapq}" : ''
        """
        # run sniffles
        sniffles --reference $ref --input $haplotagged_bam --threads ${task.cpus} --sample-id $sample_id --vcf sv.phased.vcf.gz --output-rnames --minsvlen 50 --phase $tandem_repeat_optional $mapq_optional
        """

    stub:
        """
        touch sv.phased.vcf.gz
        touch sv.phased.vcf.gz.tbi
        """

}

process cutesv {

    publishDir "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$sample_id.${ref_name}.cutesv.$filename" }, pattern: 'sv.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(haplotagged_bam), path(haplotagged_bam_index), val(data_type), val(family_position)
        val ref
        val ref_index
        val tandem_repeat
        val outdir
        val outdir2
        val ref_name
        val sv_mapq

    output:
        tuple val(sample_id), val(family_id), path("sv.vcf.gz")
        tuple val(sample_id), val(family_id), path("sv.vcf.gz"), path("sv.vcf.gz.tbi")
        tuple val(sample_id), val(family_id), val(family_position), path("sv.vcf.gz"), path(haplotagged_bam)

    script:
        // conditionally define platform specific settings
        def settings = data_type == 'ont'
            ? '--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3'
            : '--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5'
        // optionally set mapq filter threshold
        def mapq_optional = sv_mapq != 'NONE' ? "--min_mapq ${sv_mapq}" : ''
        """
        # run cuteSV
        cuteSV $haplotagged_bam $ref sv.vcf ./ --sample ${sample_id} -t ${task.cpus} --genotype --report_readid --min_size 50 $settings $mapq_optional
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

process jasmine {

    publishDir "$outdir/$family_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$family_id.$ref_name.${sv_caller}.jasmine.$filename" }, pattern: '*.vcf.gz*'

    input:
        tuple val(proband_sample_id), val(family_id), val(sample_ids), val(family_positions), path(sv_vcfs, stageAs: 'sv?.vcf.gz'), path(bams, stageAs: 'bam?.bam'), val(data_type), val(sv_caller)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(proband_sample_id), val(family_id), val(sv_caller), path("*.vcf.gz")
        tuple val(family_id), val(sv_caller), path("*.vcf.gz"), path("*.vcf.gz.tbi")

    script:
        def out_vcf = sv_caller == 'sniffles' ? 'sv.phased' : 'sv'
        // conditionally define iris arguments (by default, iris will pass minimap -x map-ont, the --pacbio flag passed to iris will pass minimap -x map-pb)
        def iris_args = '--run_iris iris_args=min_ins_length=20,--rerunracon,--keep_long_variants' + (data_type == 'pacbio' ? ',--pacbio' : '')
        """
        # gunzip vcfs and create file lists (proband first for --require_first_sample)
        FAMILY_POSITIONS=(${family_positions.join(' ')})
        for i in \${!FAMILY_POSITIONS[@]}; do
            [[ \${FAMILY_POSITIONS[\$i]} == "proband" ]] && { gunzip -c sv\$((i+1)).vcf.gz > \${FAMILY_POSITIONS[\$i]}.${out_vcf}.vcf; realpath \${FAMILY_POSITIONS[\$i]}.${out_vcf}.vcf >> vcfs.txt; realpath bam\$((i+1)).bam >> bams.txt; }
        done
        for i in \${!FAMILY_POSITIONS[@]}; do
            [[ \${FAMILY_POSITIONS[\$i]} != "proband" ]] && { gunzip -c sv\$((i+1)).vcf.gz > \${FAMILY_POSITIONS[\$i]}.${out_vcf}.vcf; realpath \${FAMILY_POSITIONS[\$i]}.${out_vcf}.vcf >> vcfs.txt; realpath bam\$((i+1)).bam >> bams.txt; }
        done
        # run jasmine
        jasmine threads=${task.cpus} out_dir=./ genome_file=$ref file_list=vcfs.txt bam_list=bams.txt out_file=${out_vcf}.tmp.vcf min_support=1 --mark_specific spec_reads=7 spec_len=20 --pre_normalize --output_genotypes --clique_merging --dup_to_ins --normalize_type --require_first_sample --default_zero_genotype $iris_args
        # fix vcf header (remove prefix to sample names that jasmine adds)
        grep '##' ${out_vcf}.tmp.vcf > ${out_vcf}.vcf
        grep '#CHROM' ${out_vcf}.tmp.vcf | sed -E 's/\t[0-9]+_/\t/g' >> ${out_vcf}.vcf
        grep -v '#' ${out_vcf}.tmp.vcf >> ${out_vcf}.vcf
        # sort
        bcftools sort ${out_vcf}.vcf -o ${out_vcf}.vcf
        # compress and index vcf
        bgzip -@ ${task.cpus} ${out_vcf}.vcf
        tabix ${out_vcf}.vcf.gz
        """

    stub:
        def out_vcf = sv_caller == 'sniffles' ? 'sv.phased' : 'sv'
        """
        touch ${out_vcf}.vcf.gz
        touch ${out_vcf}.vcf.gz.tbi
        """

}

process vep_sv {

    publishDir { task ->
        if (mode in ['singleton', 'NONE']) {
            return "$outdir/${family_id != 'NONE' ? family_id : sample_id}/$outdir2/$sample_id"
        } else {
            return "$outdir/$family_id/$outdir2"
        }
    }, mode: params.publish_mode, overwrite: true, saveAs: { filename ->
        if (mode in ['singleton', 'NONE']) {
            return "$sample_id.$ref_name.$sv_caller.$filename"
        } else {
            return "$family_id.$ref_name.${sv_caller}.jasmine.$filename"
        }
    }, pattern: '*.annotated.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(sv_vcf), val(sv_caller)
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
        tuple val(sample_id), val(family_id), path("*.annotated.vcf.gz"), path("*.annotated.vcf.gz.tbi")

    script:
        // conditionally define output sv caller specific filename
        def out_vcf = sv_caller == 'sniffles' ? 'sv.phased.annotated' : 'sv.annotated'
        """
        # run vep
        vep -i $sv_vcf -o ${out_vcf}.vcf.gz --format vcf --vcf --fasta $ref --dir $vep_db --assembly GRCh38 --species homo_sapiens --cache --offline --merged \
            --sift b --polyphen b --symbol --hgvs --hgvsg --uploaded_allele --check_existing --filter_common --distance 0 --nearest gene --canonical --mane --pick \
            --fork ${task.cpus} --no_stats --compress_output bgzip --dont_skip --plugin CADD,sv=$cadd_sv_db
        # index vcf
        tabix ${out_vcf}.vcf.gz
        """

    stub:
        def out_vcf = sv_caller == 'sniffles' ? 'sv.phased.annotated' : 'sv.annotated'
        """
        touch ${out_vcf}.vcf.gz
        touch ${out_vcf}.vcf.gz.tbi
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
    clair3_config = "${params.clair3_config}".trim()
    sv_caller = "${params.sv_caller}".trim()
    sv_mapq = "${params.sv_mapq}".trim()
    annotate = "${params.annotate}".trim()
    haploidaware = "${params.haploidaware}".trim()
    annotate_override = "${params.annotate_override}".trim()
    calculate_depth = "${params.calculate_depth}".trim()
    analyse_base_mods = "${params.analyse_base_mods}".trim()
    tr_calling = "${params.tr_calling}".trim()
    tr_call_regions = "${params.tr_call_regions}".trim()
    check_relatedness = "${params.check_relatedness}".trim()
    sites = "${params.sites}".trim()
    somatic_calling = "${params.somatic_calling}".trim()
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
    ref_name = file(ref).getSimpleName()

    // check user provided parameters
    // check for empty entries
    [in_data: in_data, ref: ref, ref_index: ref_index, outdir: outdir].each { param, val ->
        if (!val) {
            exit 1, "No value provided for '${param}'."
        }
    }
    [tandem_repeat: tandem_repeat, clair3_config: clair3_config, sites: sites].each { param, val ->
        if (!val) {
            exit 1, "No value provided for '${param}'. Set to 'NONE' if not required."
        }
    }
    // check file existence
    def required_files = [in_data: in_data, ref: ref, ref_index: ref_index]
    def optional_files = [tandem_repeat: tandem_repeat, tr_call_regions: tr_call_regions, parbed: parbed, clair3_config: clair3_config, sites: sites]
    required_files.each { param, val ->
        if (!file(val).exists()) {
            exit 1, "File does not exist, ${param} = '${val}' provided."
        }
    }
    optional_files.each { param, val ->
        if (val != 'NONE' && !file(val).exists()) {
            exit 1, "File does not exist, ${param} = '${val}' provided. Set to 'NONE' if not required."
        }
    }
    if (annotate == 'yes') {
        [vep_db: vep_db, revel_db: revel_db, gnomad_db: gnomad_db, clinvar_db: clinvar_db, cadd_snv_db: cadd_snv_db, cadd_indel_db: cadd_indel_db, spliceai_snv_db: spliceai_snv_db, spliceai_indel_db: spliceai_indel_db, alphamissense_db: alphamissense_db].each { param, val ->
            if (!file(val).exists()) {
                exit 1, "Annotation database file does not exist, ${param} = '${val}' provided."
            }
        }
    }
    // check parameter constraints and cross-parameter compatibility
    if (!(in_data_format in ['ubam_fastq', 'aligned_bam', 'snp_indel_vcf', 'sv_vcf'])) {
        exit 1, "In data format should be 'ubam_fastq', 'aligned_bam', 'snp_indel_vcf' or 'sv_vcf', in_data_format = '${in_data_format}' provided."
    }
    [haploidaware: haploidaware, annotate: annotate, calculate_depth: calculate_depth, analyse_base_mods: analyse_base_mods, tr_calling: tr_calling, check_relatedness: check_relatedness, somatic_calling: somatic_calling].each { param, val ->
        if (!(val in ['yes', 'no'])) {
            exit 1, "'${param}' should be either 'yes' or 'no', ${param} = '${val}' provided."
        }
    }
    if (haploidaware == "yes") {
        if (mode != 'singleton') {
            exit 1, "Haploid-aware mode is only supported in singleton mode, mode = '${mode}' provided."
        }
        if (sex != "XY") {
            exit 1, "Haploid-aware mode is only supported when sex = 'XY', sex = '${sex}' provided."
        }
        if (parbed == "NONE") {
            exit 1, "In haploid-aware mode, provide a valid PAR BED file, parbed = 'NONE' provided."
        }
        def ref_index_content = file(ref_index).text
        if (!ref_index_content.contains(chrXseq) || !ref_index_content.contains(chrYseq)) {
            exit 1, "Haploid-aware mode requires both chrX and chrY to be present in ${ref_index}."
        }
    }
    else if (haploidaware == "no") {
        // sex can be anything, including unset
        if (parbed != "NONE") {
            exit 1, "When not in haploid-aware mode, set the PAR BED file to 'NONE', haploidaware = '${haploidaware}' and parbed = '${parbed}' provided."
        }
    }
    if (in_data_format in ['ubam_fastq', 'aligned_bam']) {
        if (!(mode in ['singleton', 'duo', 'trio'])) {
            exit 1, "Mode should be 'singleton', 'duo' or 'trio', mode = '${mode}' provided."
        }
        if (mode in ['singleton', 'duo'] && !(snp_indel_caller in ['clair3', 'deepvariant'])) {
            exit 1, "When in ${mode} mode, the SNP/indel caller should be either 'clair3' or 'deepvariant', mode = '${mode}' and snp_indel_caller = '${snp_indel_caller}' provided."
        }
        if (mode == 'trio' && !(snp_indel_caller in ['clair3', 'deeptrio'])) {
            exit 1, "When in trio mode, the SNP/indel caller should be either 'clair3' or 'deeptrio', mode = '${mode}' and snp_indel_caller = '${snp_indel_caller}' provided."
        }
        if (snp_indel_caller == 'clair3' && !clair3_config) {
            exit 1, "When clair3 is selected as the SNP/indel calling software, provide a path to an appropriate clair3 config file for GLnexus (clair3_config), snp_indel_caller = '${snp_indel_caller}' provided."
        }
        if (!(sv_caller in ['cutesv', 'sniffles', 'both'])) {
            exit 1, "SV calling software should be 'sniffles', 'cutesv', or 'both', sv_caller = '${sv_caller}' provided."
        }
        if (sv_mapq != 'NONE' && !(sv_mapq.isInteger() && sv_mapq.toInteger() <= 60)) {
            exit 1, "sv_mapq should be a positive integer (maximum 60) or 'NONE', sv_mapq = '${sv_mapq}' provided."
        }
    }
    if (annotate == 'yes') {
        if (!ref.toLowerCase().contains('hg38') && !ref.toLowerCase().contains('grch38') && annotate_override != 'yes') {
            exit 1, "Only hg38/GRCh38 is supported for annotation. It looks like you may not be passing a hg38/GRCh38 reference genome based on the filename of the reference genome. ref = '${ref}' provided. Pass '--annotate_override yes' on the command line to override this error."
        }
    }
    if (tr_calling == 'yes' && tr_call_regions == 'NONE') {
        exit 1, "When calling tandem repeats, provide a valid tandem repeat call regions file, tr_calling = '${tr_calling}' and tr_call_regions = 'NONE' provided."
    }
    else if (tr_calling == 'no' && tr_call_regions != 'NONE') {
        exit 1, "When not calling tandem repeats, set tandem repeat call regions file to 'NONE', tr_calling = '${tr_calling}' and tr_call_regions = '${tr_call_regions}' provided."
    }
    if (check_relatedness == 'yes' && sites == 'NONE') {
        exit 1, "When checking relatedness, set an appropriate sites file. sites = '${sites}' provided."
    }
    else if (check_relatedness == 'no' && sites != 'NONE') {
        exit 1, "When not checking relatedness, set sites file to 'NONE', check_relatedness = '${check_relatedness}' and sites = '${sites}' provided."
    }
    if (in_data_format in ['snp_indel_vcf', 'sv_vcf']) {
        [annotate: [annotate, 'yes'], mode: [mode, 'NONE'], tandem_repeat: [tandem_repeat, 'NONE'], calculate_depth: [calculate_depth, 'no'], analyse_base_mods: [analyse_base_mods, 'no'], tr_calling: [tr_calling, 'no'], check_relatedness: [check_relatedness, 'no'], somatic_calling: [somatic_calling, 'no']].each { param, vals ->
            def (val, expected) = vals
            if (val != expected) {
                exit 1, "When the in data format is SNP/indel VCF or SV VCF, set ${param} to '${expected}', in_data_format = '${in_data_format}' and ${param} = '${val}' provided."
            }
        }
    }
    if (in_data_format == 'snp_indel_vcf') {
        if (sv_caller != 'NONE') {
            exit 1, "When the in data format is SNP/indel VCF, set the SV calling software to 'NONE', sv_caller = '${sv_caller}' provided."
        }
        if (snp_indel_caller == 'NONE') {
            exit 1, "When the in data format is SNP/indel VCF, pass the SNP/indel calling software which was used to generate the input data (not 'NONE'), snp_indel_caller = '${snp_indel_caller}' provided."
        }
        if (ref == 'NONE' || ref_index == 'NONE') {
            exit 1, "When the in data format is SNP/indel VCF, pass the reference genome and its index which was used to generate the input data (not 'NONE'), ref = '${ref}' and ref_index = '${ref_index}' provided."
        }
    }
    else if (in_data_format == 'sv_vcf') {
        if (snp_indel_caller != 'NONE') {
            exit 1, "When the in data format is SV VCF, set the SNP/indel calling software to 'NONE', snp_indel_caller = '${snp_indel_caller}' provided."
        }
        if (sv_caller == 'NONE') {
            exit 1, "When the in data format is SV VCF, pass the SV calling software which was used to generate the input data (not 'NONE'), sv_caller = '${sv_caller}' provided."
        }
        if (sv_caller != 'NONE' && !(sv_caller in ['sniffles', 'cutesv', 'both'])) {
            exit 1, "SV calling software should be 'sniffles', 'cutesv', or 'both', sv_caller = '${sv_caller}' provided."
        }
    }
    if (!(params.publish_mode in ['copy', 'copyNoFollow', 'link', 'move', 'rellink', 'symlink'])) {
        exit 1, "Choice of publishing mode should be 'copy', 'copyNoFollow', 'link', 'move', 'rellink' or 'symlink', publish_mode = '$params.publish_mode' provided."
    }

    // read in data
    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            if (row.file.isEmpty()) {
                exit 1, "There is an empty entry in the 'file' column of '${in_data}'."
            }
            row
        }
        .multiMap { row ->
            in_data:                     tuple(row.sample_id, row.family_id, row.file, row.data_type, row.regions_of_interest, row.clair3_model)
            id:                          tuple(row.sample_id, row.family_id)
            family_position:             tuple(row.sample_id, row.family_id, row.family_position)
            extension:                   tuple(row.sample_id, row.family_id, file(row.file).getExtension())
            files:                       tuple(row.sample_id, row.family_id, row.file)
            index:                       tuple(row.sample_id, row.family_id, "${row.file}.bai")
            data_type:                   tuple(row.sample_id, row.family_id, row.data_type)
            regions_of_interest:         tuple(row.sample_id, row.family_id, row.regions_of_interest)
            clair3_model:                tuple(row.sample_id, row.family_id, row.clair3_model)
            clairs_to_platform:             tuple(row.sample_id, row.family_id, row.clairs_to_platform)
            clairs_to_validation:        tuple(row.sample_id, row.family_id, row.clairs_to_platform)
            row_validation:              tuple(row.sample_id, row.family_id, row.family_position, row.file, row.data_type, row.regions_of_interest, row.clair3_model)
            family_validation:           tuple(row.sample_id, row.family_id, row.family_position, row.file, row.data_type, row.regions_of_interest, row.clair3_model)
            sample_family_validation:    tuple(row.sample_id, row.family_id)
            family_data_type_validation: tuple(row.sample_id, row.family_id, row.data_type)
        }
        .set { csv }

    // build channels
    csv.in_data
        .map { sample_id, family_id, files, data_type, regions_of_interest, clair3_model ->
            if (haploidaware == 'yes' && regions_of_interest != 'NONE' && file(regions_of_interest).exists()) {
                def content = file(regions_of_interest).text
                if (!content.contains(chrXseq) || !content.contains(chrYseq)) {
                    exit 1, "Haploid-aware mode requires both chrX and chrY to be present in ${regions_of_interest}."
                }
            }
            return tuple(sample_id, family_id, files, data_type, regions_of_interest, clair3_model)
        }
        .groupTuple(by: [0,1,3,4,5])
        .set { in_data_ch }

    csv.id
        .groupTuple(by: [0,1])
        .set { id_ch }

    csv.family_position
        .groupTuple(by: [0,1,2])
        .set { family_position_ch }

    csv.extension
        .groupTuple(by: [0,1,2])
        .set { extension_ch }

    csv.files
        .groupTuple(by: [0,1])
        .set { files_ch }

    csv.index
        .groupTuple(by: [0,1])
        .set { index_ch }

    csv.data_type
        .groupTuple(by: [0,1,2])
        .set { data_type_ch }

    csv.regions_of_interest
        .groupTuple(by: [0,1,2])
        .set { regions_of_interest_ch }

    csv.clair3_model
        .groupTuple(by: [0,1,2])
        .set { clair3_model_ch }

    csv.clairs_to_platform
        .groupTuple(by: [0,1,2])
        .set { clairs_to_platform_ch }

    // check user provided in data
    csv.row_validation
        .map { sample_id, family_id, family_position, in_file, data_type, regions_of_interest, clair3_model ->
            // check for empty entries
            def required_cols = [sample_id: sample_id, data_type: data_type]
            def optional_cols = [family_id: family_id, regions_of_interest: regions_of_interest, clair3_model: clair3_model]
            required_cols.each { col, val ->
                if (val.isEmpty()) {
                    exit 1, "There is an empty entry in the '${col}' column of '${in_data}'."
                }
            }
            optional_cols.each { col, val ->
                if (val.isEmpty()) {
                    exit 1, "There is an empty entry in the '${col}' column of '${in_data}'. Set to 'NONE' if not required."
                }
            }
            if (mode == 'singleton' && family_position.isEmpty()) {
                exit 1, "There is an empty entry in the 'family_position' column of '${in_data}'. Set to 'NONE' if not required."
            }
            // check file existence
            if (!file(in_file).exists()) {
                exit 1, "There is an entry in the 'file' column of '${in_data}' which doesn't exist. Check file '${in_file}'."
            }
            [regions_of_interest: regions_of_interest, clair3_model: clair3_model].each { col, val ->
                if (val != 'NONE' && !file(val).exists()) {
                    exit 1, "There is an entry in the '${col}' column of '${in_data}' which doesn't exist. Check file '${val}'."
                }
            }
            // check cross-column constraints
            if (mode in ['duo', 'trio'] && family_id == 'NONE') {
                exit 1, "Entries in the 'family_id' column of '${in_data}' should not be 'NONE' in duo/trio mode, 'NONE' provided for sample '${sample_id}'."
            }
            if (clair3_model != 'NONE') {
                if (snp_indel_caller != 'clair3') {
                    exit 1, "Pass 'NONE' in the 'clair3_model' column of '${in_data}' when clair3 is NOT selected as the SNP/indel calling software, '${clair3_model}' provided."
                }
                if (in_data_format in ['snp_indel_vcf', 'sv_vcf']) {
                    exit 1, "When the in data format is SNP/indel VCF or SV VCF, set the Clair3 model (clair3_model column of '${in_data}') to 'NONE'."
                }
            }
            // check value constraints and cross-parameter compatibility
            if (sample_id == 'NONE') {
                exit 1, "Entries in the 'sample_id' column of '${in_data}' should not be 'NONE'."
            }
            if (in_data_format == 'ubam_fastq') {
                if (!(file(in_file).getExtension() in ['bam', 'gz', 'fastq'])) {
                    exit 1, "There is an entry in the 'file' column of '$in_data' which doesn't have a 'bam', 'gz' or 'fastq' file extension. Check file '${in_file}'."
                }
                if (snp_indel_caller == 'clair3' && clair3_model == 'NONE') {
                    exit 1, "When clair3 is selected as the SNP/indel calling software, provide a path to an appropriate clair3 model in the 'clair3_model' column of '${in_data}' rather than setting it to 'NONE'."
                }
            }
            else if (in_data_format == 'aligned_bam') {
                if (!file(in_file + '.bai').exists()) {
                    exit 1, "You've specified that the in data format is aligned BAM, but it looks like '${in_file}' defined in the file column of '${in_data}' isn't indexed (ie. '${in_file}.bai' doesn't exist)."
                }
                if (!in_file.contains('bam') && in_data_format_override != 'yes') {
                    exit 1, "You've specified that the in data format is aligned BAM, but it looks like you may not be passing a BAM file based on the file names defined in the file column of '${in_data}'. Check file '${in_file}'. Pass '--in_data_format_override yes' on the command line to override this error."
                }
            }
            else if (in_data_format in ['snp_indel_vcf', 'sv_vcf']) {
                if (data_type != 'NONE') {
                    exit 1, "When the in data format is SNP/indel VCF or SV VCF, set the data type (data_type column of '${in_data}') to 'NONE'."
                }
                if (!in_file.contains('vcf') && in_data_format_override != 'yes') {
                    exit 1, "You've specified that the in data format is SNP/indel VCF or SV VCF, but it looks like you may not be passing a VCF file based on the file names defined in the file column of '${in_data}'. ${in_file} passed. Pass '--in_data_format_override yes' on the command line to override this error."
                }
            }
            if (in_data_format in ['ubam_fastq', 'aligned_bam']) {
                if (!(data_type in ['ont', 'pacbio'])) {
                    exit 1, "Entries in the 'data_type' column of '${in_data}' should be 'ont' or 'pacbio', '${data_type}' provided."
                }
            }
        }

    csv.row_validation
        .groupTuple(by: 0)
        .map { sample_id, family_ids, family_positions, files, data_types, regions_of_interests, clair3_models ->
            if (family_ids.unique().size() > 1) {
                exit 1, "All entries for a given 'sample_id' in '${in_data}' should have the same 'family_id', conflicting 'family_id' values '${family_ids}' provided for sample '${sample_id}'."
            }
        }

    csv.clairs_to_validation
        .map { sample_id, family_id, clairs_to_platform ->
            // check for empty entries
            if (clairs_to_platform.isEmpty()) {
                exit 1, "There is an empty entry in the 'clairs_to_platform' column of '${in_data}'. Set to 'NONE' if not required."
            }
            // check value constraints and cross-parameter compatibility
            if (somatic_calling == 'yes' && in_data_format in ['ubam_fastq', 'aligned_bam']) {
                if (clairs_to_platform == 'NONE') {
                    exit 1, "When somatic_calling is 'yes', provide a ClairS-TO platform in the 'clairs_to_platform' column of '${in_data}' rather than setting it to 'NONE'."
                }
            }
            else if (clairs_to_platform != 'NONE') {
                exit 1, "Set the 'clairs_to_platform' column of '${in_data}' to 'NONE' when somatic_calling is not 'yes', '${clairs_to_platform}' provided."
            }
        }

    csv.sample_family_validation
        .groupTuple(by: 0)
        .map { sample_id, family_ids ->
            if (family_ids.unique().size() > 1) {
                exit 1, "All entries for a given 'sample_id' in '${in_data}' should have the same 'family_id', conflicting 'family_id' values '${family_ids}' provided for sample '${sample_id}'."
            }
        }

    csv.family_data_type_validation
        .groupTuple(by: 1)
        .map { sample_ids, family_id, data_types ->
            if (family_id != 'NONE' && data_types.unique().size() > 1) {
                exit 1, "All entries for a given 'family_id' in '${in_data}' should have the same 'data_type', conflicting 'data_type' values '${data_types}' provided for family '${family_id}'."
            }
        }

    if (mode in ['duo', 'trio']) {
        csv.family_validation
            .groupTuple(by: 1)
            .map { sample_ids, family_ids, family_positions, files, data_types, regions_of_interests, clair3_models ->
                def unique_positions = family_positions.unique()
                if (mode == "duo") {
                    if (sample_ids.unique().size() != 2) {
                        exit 1, "In duo mode, entries in the 'sample_id' column of '$in_data' should contain 2 unique values for every 'family_id', '$sample_ids' provided for family '$family_ids'."
                    }
                    if (!('proband' in unique_positions && unique_positions.any { it in ['father', 'mother'] })) {
                        exit 1, "In duo mode, entries in the 'family_position' column of '$in_data' should contain 'proband' and 'father' or 'mother' for every 'family_id', '$unique_positions' provided for family '$family_ids'."
                    }
                }
                if (mode == "trio") {
                    if (sample_ids.unique().size() != 3) {
                        exit 1, "In trio mode, entries in the 'sample_id' column of '$in_data' should contain 3 unique values for every 'family_id', '$sample_ids' provided for family '$family_ids'."
                    }
                    if (!(['proband', 'father', 'mother'].every { it in unique_positions })) {
                        exit 1, "In trio mode, entries in the 'family_position' column of '$in_data' should contain 'proband', 'father' and 'mother' for every 'family_id', '$unique_positions' provided for family '$family_ids'."
                    }
                }
        }
    }

    // helpers
    // sort family members proband first (then father, then mother) and key tuple by proband sample_id
    def group_to_proband = { sample_ids, family_id, family_positions, bams, bam_indices, gvcfs ->
        def proband_sample_id = sample_ids[family_positions.indexOf('proband')]
        def position_order = ['proband', 'father', 'mother']
        def indices = family_positions.collect { position_order.indexOf(it) }
            .withIndex()
            .sort { a, b -> a[0] <=> b[0] }
            .collect { it[1] }
        tuple(proband_sample_id, family_id,
            indices.collect { sample_ids[it] },
            indices.collect { family_positions[it] },
            indices.collect { bams[it] },
            indices.collect { bam_indices[it] },
            indices.collect { gvcfs[it] })
    }
    // group sv vcfs by family, sort proband first, key by proband sample_id and attach sv_caller
    def group_family_sv = { ch, sv_caller_val ->
        ch
            .groupTuple(by: 1)
            .map { sample_ids, family_id, family_positions, sv_vcfs, bams ->
                def position_order = ['proband', 'father', 'mother']
                def indices = family_positions.collect { position_order.indexOf(it) }
                    .withIndex()
                    .sort { a, b -> a[0] <=> b[0] }
                    .collect { it[1] }
                def proband_sample_id = sample_ids[family_positions.indexOf('proband')]
                tuple(proband_sample_id, family_id,
                    indices.collect { sample_ids[it] },
                    indices.collect { family_positions[it] },
                    indices.collect { sv_vcfs[it] },
                    indices.collect { bams[it] })
            }
            .join(data_type_ch, by: [0, 1])
            .map { proband_sample_id, family_id, sample_ids, family_positions, sv_vcfs, bams, data_type ->
                tuple(proband_sample_id, family_id, sample_ids, family_positions, sv_vcfs, bams, data_type, sv_caller_val)
            }
    }
    // sort bams proband first (then father, then mother) for joint longtr
    def sort_longtr_family = { family_id, sample_ids, bams, bam_indices, data_types, family_positions ->
        def position_order = ['proband', 'father', 'mother']
        def ordered_indices = position_order.findAll { family_positions.contains(it) }.collect { family_positions.indexOf(it) }
        tuple(family_id,
            ordered_indices.collect { sample_ids[it] },
            ordered_indices.collect { bams[it] },
            ordered_indices.collect { bam_indices[it] },
            data_types[0])
    }
    // reorder joint sv vcf tuple
    def reorder_joint_sv_vcf = { proband_sample_id, family_id, sv_caller_val, vcf ->
        tuple(proband_sample_id, family_id, vcf, sv_caller_val)
    }

    // workflow
    // pre process
    scrape_settings(in_data_ch.join(family_position_ch, by: [0,1]), pipeface_version, in_data, in_data_format, ref, ref_index, tandem_repeat, mode, snp_indel_caller, sv_caller, annotate, calculate_depth, analyse_base_mods, tr_calling, tr_call_regions, check_relatedness, sites, somatic_calling, outdir, outdir2, haploidaware, sex, parbed, sv_mapq)
    // merge runs and alignment
    if (in_data_format == 'ubam_fastq') {
        merged = merge_runs(id_ch.join(extension_ch, by: [0,1]).join(files_ch, by: [0,1]))
        bam = minimap2(merged.join(extension_ch, by: [0,1]).join(data_type_ch, by: [0,1]), ref, ref_index)
    }
    else if (in_data_format == 'aligned_bam') {
        bam = id_ch.join(files_ch, by: [0,1]).join(index_ch, by: [0,1])
    }
    // qc
    if (in_data_format in ['ubam_fastq', 'aligned_bam']) {
        if (calculate_depth == 'yes') {
            mosdepth(bam.join(regions_of_interest_ch, by: [0,1]), outdir, outdir2, ref_name)
        }
        // somatic calling
        if (somatic_calling == 'yes') {
            clairs_to(bam.join(clairs_to_platform_ch, by: [0,1]).join(regions_of_interest_ch, by: [0,1]), ref, ref_index, outdir, outdir2, ref_name)
        }
        // snp/indel calling
        if (snp_indel_caller == 'clair3') {
            if (haploidaware == 'no' || sex == 'XX') {
                (snp_indel_vcf_bam, gvcf) = clair3(bam.join(data_type_ch, by: [0,1]).join(regions_of_interest_ch, by: [0,1]).join(clair3_model_ch, by: [0,1]), ref, ref_index, outdir, outdir2, ref_name, snp_indel_caller)
            }
            else if (haploidaware == 'yes' && sex == 'XY') {
                bam_diploid_haploid_bed = clair3_pre_processing(bam.join(regions_of_interest_ch, by: [0,1]), ref, ref_index, parbed)
                haploid_diploid_vcf_gvcf = clair3_haploid_aware(bam_diploid_haploid_bed.join(data_type_ch, by: [0,1]).join(clair3_model_ch, by: [0,1]), ref, ref_index, outdir, outdir2, ref_name)
                (snp_indel_vcf_bam, gvcf) = clair3_post_processing(haploid_diploid_vcf_gvcf, outdir, outdir2, ref_name, snp_indel_caller)
            }
        }
        else if (snp_indel_caller in ['deepvariant', 'deeptrio']) {
            dv_commands = deepvariant_dry_run(bam.join(data_type_ch, by: [0,1]), ref, ref_index, sex, haploidaware, parbed)
            dv_examples = deepvariant_make_examples(dv_commands.join(regions_of_interest_ch, by: [0,1]), ref, ref_index, parbed)
            dv_calls = deepvariant_call_variants(dv_examples)
            (snp_indel_raw_vcf_bam, snp_indel_gvcf_bam, gvcf) = deepvariant_post_processing(dv_calls.join(family_position_ch, by: [0,1]), ref, ref_index, outdir, outdir2, ref_name, snp_indel_caller)
            // filter refcall variants
            snp_indel_vcf_bam = filter_ref_call(snp_indel_raw_vcf_bam)
        }
        // split multiallelic variants
        snp_indel_split_vcf_bam = split_multiallele(snp_indel_vcf_bam, ref, ref_index)
        // phasing
        (snp_indel_split_phased_vcf_bam, snp_indel_split_phased_vcf, phased_read_list) = whatshap_phase(snp_indel_split_vcf_bam, ref, ref_index, outdir, outdir2, ref_name, snp_indel_caller)
        // haplotagging
        (haplotagged_bam, haplotagged_bam_fam, haplotagged_tsv) = whatshap_haplotag(snp_indel_split_phased_vcf_bam.join(family_position_ch, by: [0,1]), ref, ref_index, outdir, outdir2, ref_name)
        // base mod analysis
        if (analyse_base_mods == 'yes') {
            minimod(haplotagged_bam.join(data_type_ch, by: [0,1]), ref, ref_index, outdir, outdir2, ref_name)
        }
        // somalier extract (singleton mode)
        if (check_relatedness == 'yes' && mode == 'singleton') {
            somalier_extract(haplotagged_bam, ref, ref_index, sites, outdir, outdir2, ref_name)
        }
        if (mode == 'trio') {
            family_bam_by_position = haplotagged_bam_fam.groupTuple(by: 1).transpose()
            proband_bam = family_bam_by_position.filter { tuple -> tuple[2].contains("proband") }
            father_bam = family_bam_by_position.filter { tuple -> tuple[2].contains("father") }
            mother_bam = family_bam_by_position.filter { tuple -> tuple[2].contains("mother") }
            // joint snp/indel calling
            dt_commands = deeptrio_dry_run(proband_bam.join(data_type_ch, by: [0,1]), father_bam.join(data_type_ch, by: [0,1]), mother_bam.join(data_type_ch, by: [0,1]), ref, ref_index)
            dt_examples = deeptrio_make_examples(dt_commands, ref, ref_index)
            dt_calls = deeptrio_call_variants(dt_examples.proband.mix(dt_examples.father, dt_examples.mother))
            snp_indel_gvcf_bam = deeptrio_postprocessing(dt_calls, ref, ref_index)
        }
        if (mode in ['duo', 'trio']) {
            if (snp_indel_caller == 'clair3') {
                clair3_gvcf_bam = snp_indel_vcf_bam
                    .map { sample_id, family_id, bam, bam_index, vcf, vcf_index -> tuple(sample_id, family_id, bam, bam_index) }
                    .join(gvcf.map { sample_id, family_id, gvcf_file, gvcf_index -> tuple(sample_id, family_id, gvcf_file) }, by: [0,1])
                    .join(family_position_ch, by: [0,1])
                // pre-process gvcfs for clair3 + glnexus compatibility
                glnexus_ready_gvcfs = glnexus_pre_processing(clair3_gvcf_bam.map { sample_id, family_id, bam, bam_index, gvcf_file, family_position -> tuple(sample_id, family_id, gvcf_file) }, ref_index)
                family_gvcf_bam = clair3_gvcf_bam
                    .map { sample_id, family_id, bam, bam_index, gvcf_file, family_position -> tuple(sample_id, family_id, family_position, bam, bam_index) }
                    .join(glnexus_ready_gvcfs, by: [0,1])
                    .groupTuple(by: 1)
                    .map(group_to_proband)
            } else {
                family_gvcf_bam = snp_indel_gvcf_bam
                    .groupTuple(by: 1)
                    .map(group_to_proband)
            }
            // somalier extract and relate (duo/trio mode)
            if (check_relatedness == 'yes') {
                somalier_files = somalier_extract(haplotagged_bam, ref, ref_index, sites, outdir, outdir2, ref_name)
                somalier_relate_input = somalier_files
                    .join(family_position_ch, by: [0,1])
                    .groupTuple(by: 1)
                    .map { sample_ids, family_id, somalier_files_list, family_positions ->
                        def proband_sample_id = sample_ids[family_positions.indexOf('proband')]
                        tuple(proband_sample_id, family_id, somalier_files_list.flatten())
                    }
                somalier_relate(somalier_relate_input, outdir, outdir2, ref_name)
            }
            // gvcf merging
            joint_snp_indel_bcf_bam = glnexus(family_gvcf_bam, snp_indel_caller, clair3_config)
            joint_snp_indel_vcf_bam = glnexus_post_processing(joint_snp_indel_bcf_bam)
            // joint split multiallelic variants
            joint_snp_indel_split_vcf_bam = split_multiallele_family(joint_snp_indel_vcf_bam, ref, ref_index)
            // joint phasing
            (joint_snp_indel_phased_vcf, joint_phased_read_list) = whatshap_phase_family(joint_snp_indel_split_vcf_bam, ref, ref_index, outdir, outdir2, ref_name, snp_indel_caller)
        }
        if (mode in ['duo', 'trio']) {
            // joint snp/indel annotation
            if (annotate == 'yes') {
                vep_snp_indel(joint_snp_indel_phased_vcf, ref, ref_index, vep_db, revel_db, gnomad_db, clinvar_db, cadd_snv_db, cadd_indel_db, spliceai_snv_db, spliceai_indel_db, alphamissense_db, outdir, outdir2, ref_name, snp_indel_caller, mode)
            }
        }
        // sv calling
        if (sv_caller in ['sniffles', 'both']) {
            (sv_vcf_sniffles, sv_vcf_sniffles_indexed, sv_vcf_haplotagged_bam_fam_sniffles) = sniffles(haplotagged_bam.join(family_position_ch, by: [0,1]), ref, ref_index, tandem_repeat, outdir, outdir2, ref_name, sv_mapq)
        }
        if (sv_caller in ['cutesv', 'both']) {
            (sv_vcf_cutesv, sv_vcf_cutesv_indexed, sv_vcf_haplotagged_bam_fam_cutesv) = cutesv(haplotagged_bam.join(data_type_ch, by: [0,1]).join(family_position_ch, by: [0,1]), ref, ref_index, tandem_repeat, outdir, outdir2, ref_name, sv_mapq)
        }
        if (tr_calling == 'yes') {
            split_beds = longtr_pre_processing(tr_call_regions)
            // tr calling
            tr_vcfs = longtr(haplotagged_bam.join(data_type_ch, by: [0,1]), split_beds, ref, ref_index)
            concat_tr_vcf(tr_vcfs, outdir, outdir2, ref_name)
            // joint tr calling
            if (mode in ['duo', 'trio']) {
                family_input = haplotagged_bam
                    .join(data_type_ch, by: [0,1])
                    .join(family_position_ch, by: [0,1])
                    .map { sample_id, family_id, bam, bam_index, data_type, family_position ->
                        tuple(family_id, sample_id, bam, bam_index, data_type, family_position)
                    }
                    .groupTuple(by: 0)
                    .map(sort_longtr_family)
                tr_vcfs_family = longtr_family(family_input, split_beds, ref, ref_index)
                concat_tr_vcf_family(tr_vcfs_family, outdir, outdir2, ref_name)
            }
        }
    }
    if (in_data_format == 'snp_indel_vcf') {
        snp_indel_split_phased_vcf = id_ch.join(files_ch, by: [0,1])
    }
    if (in_data_format in ['ubam_fastq', 'aligned_bam', 'snp_indel_vcf']) {
        // annotation
        if (annotate == 'yes' && !(mode in ['duo', 'trio'])) {
            vep_snp_indel(snp_indel_split_phased_vcf, ref, ref_index, vep_db, revel_db, gnomad_db, clinvar_db, cadd_snv_db, cadd_indel_db, spliceai_snv_db, spliceai_indel_db, alphamissense_db, outdir, outdir2, ref_name, snp_indel_caller, mode)
        }
    }
    if (mode in ['duo', 'trio']) {
        // sv vcf merging
        if (sv_caller in ['sniffles', 'both']) {
            family_sv_sniffles = group_family_sv(sv_vcf_haplotagged_bam_fam_sniffles, 'sniffles')
        }
        if (sv_caller in ['cutesv', 'both']) {
            family_sv_cutesv = group_family_sv(sv_vcf_haplotagged_bam_fam_cutesv, 'cutesv')
        }
        family_sv_all = sv_caller == 'both' ? family_sv_sniffles.mix(family_sv_cutesv) :
                        sv_caller == 'sniffles' ? family_sv_sniffles : family_sv_cutesv
        (joint_sv_vcfs, joint_sv_vcfs_indexed) = jasmine(family_sv_all, ref, ref_index, outdir, outdir2, ref_name)
        if (sv_caller in ['sniffles', 'both']) {
            joint_sv_vcf_sniffles = joint_sv_vcfs
                .filter { it[2] == 'sniffles' }
                .map(reorder_joint_sv_vcf)
        }
        if (sv_caller in ['cutesv', 'both']) {
            joint_sv_vcf_cutesv = joint_sv_vcfs
                .filter { it[2] == 'cutesv' }
                .map(reorder_joint_sv_vcf)
        }
    }
    if (in_data_format == 'sv_vcf' && !(mode in ['duo', 'trio'])) {
        sv_vcf_sniffles = id_ch.join(files_ch, by: [0,1])
        sv_vcf_cutesv = id_ch.join(files_ch, by: [0,1])
    }
    // sv annotation
    sv_vcf_for_vep = Channel.empty()
    if (annotate == 'yes') {
        if (in_data_format in ['ubam_fastq', 'aligned_bam', 'sv_vcf'] && !(mode in ['duo', 'trio'])) {
            if (sv_caller in ['sniffles', 'both']) {
                sv_vcf_for_vep = sv_vcf_for_vep.mix(sv_vcf_sniffles.map { s, f, v -> tuple(s, f, v, 'sniffles') })
            }
            if (sv_caller in ['cutesv', 'both']) {
                sv_vcf_for_vep = sv_vcf_for_vep.mix(sv_vcf_cutesv.map { s, f, v -> tuple(s, f, v, 'cutesv') })
            }
        }
        if (mode in ['duo', 'trio']) {
            if (sv_caller in ['sniffles', 'both']) {
                sv_vcf_for_vep = sv_vcf_for_vep.mix(joint_sv_vcf_sniffles)
            }
            if (sv_caller in ['cutesv', 'both']) {
                sv_vcf_for_vep = sv_vcf_for_vep.mix(joint_sv_vcf_cutesv)
            }
        }
        vep_sv(sv_vcf_for_vep, ref, ref_index, vep_db, gnomad_db, cadd_sv_db, outdir, outdir2, ref_name, mode)
    }
}
