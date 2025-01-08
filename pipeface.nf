
nextflow.enable.dsl=2

// create dummy NONE file for optional pipeface inputs
def filePath = "NONE"
new File(filePath).text = "Dummy file for optional pipeface inputs. Don't delete during a pipeline run unless you want a bad time.\n"

// set defaults for optional params (set default optional input files to dummy NONE file)
// secret sauce second outdir
params.outdir2 = ""
params.tandem_repeat = "NONE"
params.annotate_override = ""
params.in_data_format_override = ""

process scrape_settings {

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$filename" }, pattern: '*pipeface_settings.txt'

    input:
        tuple val(sample_id), val(family_id), val(extension), val(files), val(data_type), val(regions_of_interest), val(clair3_model)
        val in_data
        val in_data_format
        val ref
        val ref_index
        val tandem_repeat
        val snp_indel_caller
        val sv_caller
        val annotate
        val calculate_depth
        val outdir
        val outdir2

    output:
        tuple val(sample_id), val(family_id), path('pipeface_settings.txt')

    script:
        // conditionally define reported SV caller
        if( sv_caller == 'both' ) {
            reported_sv_caller = 'cutesv,sniffles'
        }
        else if ( sv_caller == 'sniffles' ) {
            reported_sv_caller = 'sniffles'
        }
        else if ( sv_caller == 'cutesv' ) {
            reported_sv_caller = 'cutesv'
        }
        else if ( sv_caller == 'NONE' ) {
            reported_sv_caller = 'NONE'
        }
        // conditionally define reported in data format
        if ( in_data_format == 'ubam_fastq' ) {
            reported_in_data_format = 'unaligned BAM and/or FASTQ'
        }
        else if ( in_data_format == 'aligned_bam' ) {
            reported_in_data_format = 'aligned BAM'
        }
        else if ( in_data_format == 'snv_vcf' ) {
            reported_in_data_format = 'SNP/indel VCF'
        }
        else if ( in_data_format == 'sv_vcf' ) {
            reported_in_data_format = 'SV VCF'
        }
        if ( in_data_format == 'ubam_fastq' | in_data_format == 'aligned_bam' )
        """
        echo "Sample ID: $sample_id" >> pipeface_settings.txt
        echo "Family ID: $family_id" >> pipeface_settings.txt
        echo "In data format: $reported_in_data_format" >> pipeface_settings.txt
        echo "Input data file/files: $files" >> pipeface_settings.txt
        echo "Data type: $data_type" >> pipeface_settings.txt
        echo "Regions of interest file: $regions_of_interest" >> pipeface_settings.txt
        echo "Clair3 model: $clair3_model" >> pipeface_settings.txt
        echo "In data csv path: $in_data" >> pipeface_settings.txt
        echo "Reference genome: $ref" >> pipeface_settings.txt
        echo "Reference genome index: $ref_index" >> pipeface_settings.txt
        echo "Tandem repeat file: $tandem_repeat" >> pipeface_settings.txt
        echo "SNP/indel caller: $snp_indel_caller" >> pipeface_settings.txt
        echo "SV caller: $reported_sv_caller" >> pipeface_settings.txt
        echo "Annotate: $annotate" >> pipeface_settings.txt
        echo "Calculate depth: $calculate_depth" >> pipeface_settings.txt
        echo "Outdir: $outdir" >> pipeface_settings.txt
        """
        else if( in_data_format == 'snv_vcf' | in_data_format == 'sv_vcf' )
        """
        echo "Sample ID: $sample_id" >> pipeface_settings.txt
        echo "Family ID: $family_id" >>	pipeface_settings.txt
        echo "In data format: $reported_in_data_format" >> pipeface_settings.txt
        echo "Input data file/files: $files" >> pipeface_settings.txt
        echo "In data csv path: $in_data" >> pipeface_settings.txt
        echo "Annotate: $annotate" >> pipeface_settings.txt
        echo "Outdir: $outdir" >> pipeface_settings.txt
        """

    stub:
        """
        touch pipeface_settings.txt
        """

}

process scrape_bam_header {

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$file_short.$filename" }

    input:
        tuple val(sample_id), val(family_id), val(extension), val(file), val(file_short)
        val outdir
        val outdir2

    output:
        tuple val(sample_id), val(family_id), val(file_short), path('header')

    script:
    if( extension == 'gz' )
        """
        echo "none" > header
        """
    else if ( extension == 'fastq' )
        """
        echo "none" > header
        """
    else if( extension == 'bam' )
        """
        samtools head $file > header
        """

    stub:
        """
        touch header
        """

}

process merge_runs {

    input:
        tuple val(sample_id), val(family_id), val(extension), path(files)

    output:
        tuple val(sample_id), val(family_id), path('merged')

    script:
    if( extension == 'gz' )
        """
        ${
            if (files instanceof List && files.size() > 1) {
                "cat ${files} > merged"
            } else {
                "ln -s ${files} merged"
            }
        }
        """
    else if ( extension == 'fastq' )
        """
        ${
            if (files instanceof List && files.size() > 1) {
                "cat ${files} | bgzip -@ ${task.cpus} > merged"
            } else {
                "ln -s ${files} merged"
            }
        }
        """
    else if( extension == 'bam' )
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
        tuple val(sample_id), val(family_id), path('sorted.bam')

    script:
    // conditionally define preset
    if( data_type == 'ont' ) {
        preset = 'lr:hq'
    }
    else if( data_type == 'pacbio' ) {
        preset = 'map-hifi'
    }
    if( extension == 'bam' )
        """
        # run minimap
        samtools fastq \
        -@ ${task.cpus} \
        -T '*' \
        $merged | minimap2 \
        -y \
        -Y \
        --secondary=no \
        --MD \
        -a \
        -x $preset \
        -t ${task.cpus} \
        $ref - | samtools sort -@ ${task.cpus} -o sorted.bam -
        # index bam
        samtools index \
        -@ ${task.cpus} \
        sorted.bam
        """
    else if( extension == 'gz' | extension == 'fastq' )
        """
        # run minimap
        minimap2 \
        --secondary=no \
        --MD \
        -a \
        -Y \
        -x $preset \
        -t ${task.cpus} \
        $ref \
        $merged | samtools sort -@ ${task.cpus} -o sorted.bam -
        # index bam
        samtools index \
        -@ ${task.cpus} \
        sorted.bam
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
        tuple val(sample_id), val(family_id), path(bam), val(regions_of_interest)
        val mosdepth_binary
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path('depth.txt')

    script:
    // define a string to optionally pass regions of interest bed file
    def regions_of_interest_optional = file(regions_of_interest).name != 'NONE' ? "-b $regions_of_interest" : ''
        """
        # stage bam and bam index
        # do this here instead of input tuple so I can handle processing an aligned bam as an input file without requiring a bam index for ubam input
        bam_loc=\$(realpath ${bam})
        ln -sf \${bam_loc} sorted.bam
        ln -sf \${bam_loc}.bai .
        ln -sf \${bam_loc}.bai sorted.bam.bai
        # run mosdepth
        $mosdepth_binary \
        depth \
        $bam \
        $regions_of_interest_optional \
        -t ${task.cpus}
        # rename file
        ln -s depth.mosdepth.summary.txt depth.txt
        """

    stub:
        """
        touch depth.txt
        """

}

process clair3 {

    input:
        tuple val(sample_id), val(family_id), path(bam), val(data_type), val(regions_of_interest), val(clair3_model)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(family_id), path('sorted.bam'), path('sorted.bam.bai'), path('snp_indel.vcf.gz'), path('snp_indel.vcf.gz.tbi')

    script:
    // define a string to optionally pass regions of interest bed file
    def regions_of_interest_optional = file(regions_of_interest).name != 'NONE' ? "--bed_fn=$regions_of_interest" : ''
    // conditionally define platform
    if( data_type == 'ont' ) {
        platform = 'ont'
    }
    else if( data_type == 'pacbio' ) {
        platform = 'hifi'
    }
        """
        # stage bam and bam index
        # do this here instead of input tuple so I can handle processing an aligned bam as an input file without requiring a bam index for ubam input
        bam_loc=\$(realpath ${bam})
        ln -sf \${bam_loc} sorted.bam
        ln -sf \${bam_loc}.bai .
        ln -sf \${bam_loc}.bai sorted.bam.bai
        # run clair3
        run_clair3.sh \
        --bam_fn=$bam \
        --ref_fn=$ref \
        --output=./ \
        --threads=${task.cpus} \
        --platform=$platform \
        --model_path=$clair3_model \
        --sample_name=$sample_id \
        --gvcf \
        --include_all_ctgs \
        $regions_of_interest_optional
        # rename files
        ln -s merge_output.vcf.gz snp_indel.vcf.gz
        ln -s merge_output.vcf.gz.tbi snp_indel.vcf.gz.tbi
        """

    stub:
        """
        touch snp_indel.vcf.gz
        touch snp_indel.vcf.gz.tbi
        """

}

process deepvariant_dry_run {

    input:
        tuple val(sample_id), val(family_id), path(bam), val(data_type)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(family_id), path('sorted.bam'), path('sorted.bam.bai'), env(make_examples_args), env(call_variants_args)

    script:
    // conditionally define model type
    if( data_type == 'ont' ) {
        model = 'ONT_R104'
    }
    else if ( data_type == 'pacbio' ) {
        model = 'PACBIO'
    }
        """
        # stage bam and bam index
        # do this here instead of input tuple so I can handle processing an aligned bam as an input file without requiring a bam index for ubam input
        bam_loc=\$(realpath ${bam})
        ln -sf \${bam_loc} sorted.bam
        ln -sf \${bam_loc}.bai .
        ln -sf \${bam_loc}.bai sorted.bam.bai
        # do a dry-run of deepvariant
        run_deepvariant \
        --reads=$bam \
        --ref=$ref \
        --sample_name=$sample_id \
        --output_vcf=snp_indel.raw.vcf.gz \
        --model_type=$model \
        --dry_run=true > commands.txt
        # extract arguments for make_examples and call_variants stages
        make_examples_args=\$(grep "/opt/deepvariant/bin/make_examples" commands.txt | awk '{split(\$0, arr, "--add_hp_channel"); print "--add_hp_channel" arr[2]}' | sed 's/--sample_name "[^"]*"//g')
        call_variants_args=\$(grep "/opt/deepvariant/bin/call_variants" commands.txt | awk '{split(\$0, arr, "--checkpoint"); print "--checkpoint" arr[2]}')
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

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(call_variants_args), path('*.gz{,.example_info.json}')

    script:
    // define an optional string to pass regions of interest bed file
    def regions_of_interest_optional = file(regions_of_interest).name != 'NONE' ? "--regions $regions_of_interest" : ''
        """
        seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer make_examples \\
            --mode calling --ref "${ref}" --reads "${bam}" --sample_name "${sample_id}" ${regions_of_interest_optional} --examples "make_examples.tfrecord@${task.cpus}.gz" ${make_examples_args}
        """

    stub:
        """
        touch make_examples.tfrecord-00000-of-00104.gz
        touch make_examples.tfrecord-00000-of-00104.gz.example_info.json
        """

}

process deepvariant_call_variants {

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(call_variants_args), path(make_examples_out)

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(call_variants_args), path('*.gz')

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

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), val(call_variants_args), path(call_variants_out)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path('snp_indel.vcf.gz'), path('snp_indel.vcf.gz.tbi')

    script:
        """
        # postprocess_variants and vcf_stats_report stages in deepvariant
        postprocess_variants --ref "${ref}" --infile "call_variants_output.tfrecord.gz" --outfile "snp_indel.raw.vcf.gz" --cpus "${task.cpus}" --sample_name "${sample_id}"
        vcf_stats_report --input_vcf "snp_indel.raw.vcf.gz" --outfile_base "snp_indel.raw"
        # filter out refcall variants
        bcftools view -f 'PASS' snp_indel.raw.vcf.gz -o snp_indel.vcf.gz
        # index vcf
        tabix snp_indel.vcf.gz
        """

    stub:
        """
        touch snp_indel.vcf.gz
        touch snp_indel.vcf.gz.tbi
        """

}

process whatshap_phase {

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$snp_indel_caller.$filename" }, pattern: 'snp_indel.phased.*'

    input:
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path(snp_indel_vcf), path(snp_indel_vcf_index)
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
        # run whatshap phase
        whatshap phase \
        --reference $ref \
        --output snp_indel.phased.vcf.gz \
        --output-read-list snp_indel.phased.read_list.txt \
        --sample $sample_id \
        --ignore-read-groups $snp_indel_vcf $bam
        # index vcf
        tabix snp_indel.phased.vcf.gz
        # run whatshap stats
        whatshap stats \
        snp_indel.phased.vcf.gz \
        --gtf snp_indel.phased.stats.gtf \
        --sample $sample_id
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
        tuple val(sample_id), val(family_id), path(bam), path(bam_index), path(snp_indel_vcf), path(snp_indel_vcf_index)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path('sorted.haplotagged.bam'), path('sorted.haplotagged.bam.bai')
        tuple val(sample_id), val(family_id), path('sorted.haplotagged.tsv')

    script:
        """
        # run whatshap haplotag
        whatshap haplotag \
        --reference $ref \
        --output sorted.haplotagged.bam \
        --sample $sample_id \
        --tag-supplementary \
        --ignore-read-groups \
        --output-threads ${task.cpus} \
        --output-haplotag-list sorted.haplotagged.tsv \
        $snp_indel_vcf $bam
        # index bam
        samtools index \
        -@ ${task.cpus} \
        sorted.haplotagged.bam
        """

    stub:
        """
        touch sorted.haplotagged.bam
        touch sorted.haplotagged.bam.bai
        touch sorted.haplotagged.tsv
        """

}

process vep_snv {

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$snp_indel_caller.$filename" }, pattern: 'snp_indel.phased.annotated.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(snp_indel_phased_vcf)
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

    output:
        tuple val(sample_id), val(family_id), path('snp_indel.phased.annotated.vcf.gz'), path('snp_indel.phased.annotated.vcf.gz.tbi')

    script:
        """
        # run vep
        vep -i $snp_indel_phased_vcf \
        -o snp_indel.phased.annotated.vcf.gz \
        --format vcf \
        --vcf \
        --fasta $ref \
        --dir $vep_db \
        --assembly GRCh38 \
        --species homo_sapiens \
        --cache \
        --offline \
        --merged \
        --sift b \
        --polyphen b \
        --symbol \
        --hgvs \
        --hgvsg \
        --plugin REVEL,file=$revel_db \
        --custom file=$gnomad_db,short_name=gnomAD,format=vcf,type=exact,fields=AF_joint%AF_exomes%AF_genomes%nhomalt_joint%nhomalt_exomes%nhomalt_genomes \
        --custom file=$clinvar_db,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG \
        --plugin CADD,snv=$cadd_snv_db,indels=$cadd_indel_db \
        --plugin SpliceAI,snv=$spliceai_snv_db,indel=$spliceai_indel_db \
        --plugin AlphaMissense,file=$alphamissense_db \
        --uploaded_allele \
        --check_existing \
        --filter_common \
        --no_intergenic \
        --pick \
        --fork ${task.cpus} \
        --no_stats \
        --compress_output bgzip
        # index vcf
        tabix snp_indel.phased.annotated.vcf.gz
        """

    stub:
        """
        touch snp_indel.phased.annotated.vcf.gz
        touch snp_indel.phased.annotated.vcf.gz.tbi
        """

}

process pbcpgtools {

    def software = "pbcpgtools"

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$software.$filename" }, pattern: 'cpg_scores*'

    input:
        tuple val(sample_id), val(family_id), path(haplotagged_bam), path(haplotagged_bam_index), val(data_type)
        val pbcpgtools_binary
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path('cpg_scores_hap1.bw'), path('cpg_scores_hap1.bed'), path('cpg_scores_hap2.bw'), path('cpg_scores_hap2.bed'), path('cpg_scores_combined.bw'), path('cpg_scores_combined.bed'), optional: true

    script:
    if( data_type == 'pacbio' )
        """
        # run pb-cpg-tools
        $pbcpgtools_binary/bin/aligned_bam_to_cpg_scores \
        --bam $haplotagged_bam \
        --ref $ref \
        --pileup-mode model \
        --model $pbcpgtools_binary/models/pileup_calling_model.v1.tflite \
        --modsites-mode denovo \
        --hap-tag HP \
        --threads ${task.cpus}
        # rename files
        ln -s aligned_bam_to_cpg_scores.hap1.bw cpg_scores_hap1.bw
        ln -s aligned_bam_to_cpg_scores.hap1.bed cpg_scores_hap1.bed
        ln -s aligned_bam_to_cpg_scores.hap2.bw cpg_scores_hap2.bw
        ln -s aligned_bam_to_cpg_scores.hap2.bed cpg_scores_hap2.bed
        ln -s aligned_bam_to_cpg_scores.combined.bw cpg_scores_combined.bw
        ln -s aligned_bam_to_cpg_scores.combined.bed cpg_scores_combined.bed
        """
    else if( data_type == 'ont' )
        """
        echo "Data type is ONT, not running pb-CpG-tools on this data."
        """

    stub:
        """
        touch cpg_scores_hap1.bw
        touch cpg_scores_hap1.bed
        touch cpg_scores_hap2.bw
        touch cpg_scores_hap2.bed
        touch cpg_scores_combined.bw
        touch cpg_scores_combined.bed
        """

}

process sniffles {

    def sv_caller = "sniffles"

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$sv_caller.$filename" }, pattern: 'sv.phased.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(haplotagged_bam), path(haplotagged_bam_index)
        val ref
        val ref_index
        val tandem_repeat
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path('sv.phased.vcf.gz')
        tuple val(sample_id), val(family_id), path('sv.phased.vcf.gz'), path('sv.phased.vcf.gz.tbi')

    script:
    // define a string to optionally pass tandem repeat bed file
    def tandem_repeat_optional = file(tandem_repeat).name != 'NONE' ? "--tandem-repeats $tandem_repeat" : ''
        """
        # run sniffles
        sniffles \
        --reference $ref \
        --input $haplotagged_bam \
        --threads ${task.cpus} \
        --sample-id $sample_id \
        --vcf sv.phased.vcf.gz \
        --output-rnames \
        --minsvlen 50 \
        --phase $tandem_repeat_optional
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
        tuple val(sample_id), val(family_id), path(haplotagged_bam), path(haplotagged_bam_index), val(data_type)
        val ref
        val ref_index
        val tandem_repeat
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path('sv.vcf.gz')
        tuple val(sample_id), val(family_id), path('sv.vcf.gz'), path('sv.vcf.gz.tbi')

    script:
    if( data_type == 'ont' ) {
        settings = '--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3'
    }
    else if( data_type == 'pacbio' ) {
        settings = '--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5'
    }
        """
        # run cuteSV
        cuteSV \
        $haplotagged_bam \
        $ref \
        sv.vcf \
        ./ \
        --sample ${sample_id} \
        -t ${task.cpus} \
        --genotype \
        --report_readid \
        --min_size 50 \
        $settings
        # compress and index vcf
        bgzip \
        -@ ${task.cpus} \
        sv.vcf
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

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$sv_caller.$filename" }, pattern: 'sv.phased.annotated.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(sv_phased_vcf)
        val ref
        val ref_index
        val vep_db
        val gnomad_db
        val gnomad_sv_db
        val clinvar_db
        val cadd_sv_db
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path('sv.phased.annotated.vcf.gz'), path('sv.phased.annotated.vcf.gz.tbi')

    script:
        """
        # run vep
        vep -i $sv_phased_vcf \
        -o sv.phased.annotated.vcf.gz \
        --format vcf \
        --vcf \
        --fasta $ref \
        --dir $vep_db \
        --assembly GRCh38 \
        --species homo_sapiens \
        --cache \
        --offline \
        --merged \
        --sift b \
        --polyphen b \
        --symbol \
        --hgvs \
        --hgvsg \
        --custom file=$gnomad_sv_db,short_name=gnomAD_sv,format=vcf,type=overlap,reciprocal=1,overlap_cutoff=80,same_type=1,num_records=50,fields=ALGORITHMS%BOTHSIDES_SUPPORT%CHR2%CPX_INTERVALS%CPX_TYPE%END%END2%EVIDENCE%LOW_CONFIDENCE_REPETITIVE_LARGE_DUP%MEMBERS%MULTIALLELIC%NCR%OUTLIER_SAMPLE_ENRICHED_LENIENT%PAR%PCRMINUS_NCR%PCRPLUS_NCR%PESR_GT_OVERDISPERSION%POS2%PREDICTED_BREAKEND_EXONIC%PREDICTED_COPY_GAIN%PREDICTED_DUP_PARTIAL%PREDICTED_INTERGENIC%PREDICTED_INTRAGENIC_EXON_DUP%PREDICTED_INTRONIC%PREDICTED_INV_SPAN%PREDICTED_LOF%PREDICTED_MSV_EXON_OVERLAP%PREDICTED_NEAREST_TSS%PREDICTED_NONCODING_BREAKPOINT%PREDICTED_NONCODING_SPAN%PREDICTED_PARTIAL_DISPERSED_DUP%PREDICTED_PARTIAL_EXON_DUP%PREDICTED_PROMOTER%PREDICTED_TSS_DUP%PREDICTED_UTR%RESOLVED_POSTHOC%SOURCE%SVLEN%SVTYPE%UNRESOLVED_TYPE%AN%AC%AF%N_BI_GENOS%N_HOMREF%N_HET%N_HOMALT%FREQ_HOMREF%FREQ_HET%FREQ_HOMALT%CN_NUMBER%CN_COUNT%CN_STATUS%CN_FREQ%CN_NONREF_COUNT%CN_NONREF_FREQ \
        --custom file=$clinvar_db,short_name=ClinVar,format=vcf,type=overlap,reciprocal=1,overlap_cutoff=50,same_type=1,num_records=50,fields=CLNSIG \
        --plugin CADD,sv=$cadd_sv_db \
        --uploaded_allele \
        --check_existing \
        --filter_common \
        --no_intergenic \
        --pick \
        --fork ${task.cpus} \
        --no_stats \
        --compress_output bgzip
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

    publishDir "$outdir/$family_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$sv_caller.$filename" }, pattern: 'sv.annotated.vcf.gz*'

    input:
        tuple val(sample_id), val(family_id), path(sv_vcf)
        val ref
        val ref_index
        val vep_db
        val gnomad_db
        val gnomad_sv_db
        val clinvar_db
        val cadd_sv_db
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(sample_id), val(family_id), path('sv.annotated.vcf.gz'), path('sv.annotated.vcf.gz.tbi')

    script:
        """
        # run vep
        vep -i $sv_vcf \
        -o sv.annotated.vcf.gz \
        --format vcf \
        --vcf \
        --fasta $ref \
        --dir $vep_db \
        --assembly GRCh38 \
        --species homo_sapiens \
        --cache \
        --offline \
        --merged \
        --sift b \
        --polyphen b \
        --symbol \
        --hgvs \
        --hgvsg \
        --custom file=$gnomad_sv_db,short_name=gnomAD_sv,format=vcf,type=overlap,reciprocal=1,overlap_cutoff=80,same_type=1,num_records=50,fields=ALGORITHMS%BOTHSIDES_SUPPORT%CHR2%CPX_INTERVALS%CPX_TYPE%END%END2%EVIDENCE%LOW_CONFIDENCE_REPETITIVE_LARGE_DUP%MEMBERS%MULTIALLELIC%NCR%OUTLIER_SAMPLE_ENRICHED_LENIENT%PAR%PCRMINUS_NCR%PCRPLUS_NCR%PESR_GT_OVERDISPERSION%POS2%PREDICTED_BREAKEND_EXONIC%PREDICTED_COPY_GAIN%PREDICTED_DUP_PARTIAL%PREDICTED_INTERGENIC%PREDICTED_INTRAGENIC_EXON_DUP%PREDICTED_INTRONIC%PREDICTED_INV_SPAN%PREDICTED_LOF%PREDICTED_MSV_EXON_OVERLAP%PREDICTED_NEAREST_TSS%PREDICTED_NONCODING_BREAKPOINT%PREDICTED_NONCODING_SPAN%PREDICTED_PARTIAL_DISPERSED_DUP%PREDICTED_PARTIAL_EXON_DUP%PREDICTED_PROMOTER%PREDICTED_TSS_DUP%PREDICTED_UTR%RESOLVED_POSTHOC%SOURCE%SVLEN%SVTYPE%UNRESOLVED_TYPE%AN%AC%AF%N_BI_GENOS%N_HOMREF%N_HET%N_HOMALT%FREQ_HOMREF%FREQ_HET%FREQ_HOMALT%CN_NUMBER%CN_COUNT%CN_STATUS%CN_FREQ%CN_NONREF_COUNT%CN_NONREF_FREQ \
        --custom file=$clinvar_db,short_name=ClinVar,format=vcf,type=overlap,reciprocal=1,overlap_cutoff=50,same_type=1,num_records=50,fields=CLNSIG \
        --plugin CADD,sv=$cadd_sv_db \
        --uploaded_allele \
        --check_existing \
        --filter_common \
        --no_intergenic \
        --pick \
        --fork ${task.cpus} \
        --no_stats \
        --compress_output bgzip
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
    in_data = "$params.in_data"
    in_data_format = "$params.in_data_format"
    in_data_format_override = "$params.in_data_format_override"
    ref = "$params.ref"
    ref_index = "$params.ref_index"
    tandem_repeat = "$params.tandem_repeat"
    snp_indel_caller = "$params.snp_indel_caller"
    sv_caller = "$params.sv_caller"
    annotate = "$params.annotate"
    annotate_override = "$params.annotate_override"
    calculate_depth = "$params.calculate_depth"
    outdir = "$params.outdir"
    outdir2 = "$params.outdir2"
    mosdepth_binary = "$params.mosdepth_binary"
    pbcpgtools_binary = "$params.pbcpgtools_binary"
    vep_db = "$params.vep_db"
    revel_db = "$params.revel_db"
    gnomad_db = "$params.gnomad_db"
    gnomad_sv_db = "$params.gnomad_sv_db"
    clinvar_db = "$params.clinvar_db"
    cadd_snv_db = "$params.cadd_snv_db"
    cadd_indel_db = "$params.cadd_indel_db"
    cadd_sv_db = "$params.cadd_sv_db"
    spliceai_snv_db = "$params.spliceai_snv_db"
    spliceai_indel_db = "$params.spliceai_indel_db"
    alphamissense_db = "$params.alphamissense_db"

    // check user provided parameters
    if ( !in_data ) {
        exit 1, "No in data csv file specified. Either include in parameter file or pass to --in_data on the command line."
    }
    if ( !file(in_data).exists() ) {
        exit 1, "In data csv file path does not exist, '${in_data}' provided."
    }
    if ( !in_data_format ) {
        exit 1, "No in data format selected. Either include in parameter file or pass to --in_data_format on the command line. Should be 'ubam_fastq', 'aligned_bam' or 'snv_vcf'."
    }
    if ( in_data_format != 'ubam_fastq' && in_data_format != 'aligned_bam' && in_data_format != 'snv_vcf' && in_data_format != 'sv_vcf' ) {
        exit 1, "In data format should be 'ubam_fastq', 'aligned_bam', 'snv_vcf' or 'sv_vcf', '${in_data_format}' selected."
    }
    if ( in_data_format == 'snv_vcf' && tandem_repeat != 'NONE' ) {
        exit 1, "In data format is SNP/indel VCF, but you haven't set the tandem repeat file to 'NONE'. Either set tandem_repeat to 'NONE' in parameter file or pass '--tandem_repeat NONE' on the command line"
    }
    if ( in_data_format == 'sv_vcf' && tandem_repeat != 'NONE' ) {
        exit 1, "In data format is SV VCF, but you haven't set the tandem repeat file to 'NONE'. Either set tandem_repeat to 'NONE' in parameter file or pass '--tandem_repeat NONE' on the command line"
    }
    if ( in_data_format == 'snv_vcf' && sv_caller != 'NONE' ) {
        exit 1, "In data format is SNP/indel VCF, but you haven't set the SV calling software to 'NONE'. Either set sv_caller to 'NONE' in parameter file or pass '--sv_caller NONE' on the command line"
    }
    if ( in_data_format == 'sv_vcf' && snp_indel_caller != 'NONE' ) {
        exit 1, "In data format is SV VCF, but you haven't set the SNP/indel calling software to 'NONE'. Either set snp_indel_caller to 'NONE' in parameter file or pass '--snp_indel_caller NONE' on the command line."
    }
    if ( in_data_format == 'snv_vcf' && annotate == 'no' ) {
        exit 1, "In data format is SNP/indel VCF, but you've chosen not to annotate. Nothing for pipeface to do."
    }
    if ( in_data_format == 'sv_vcf' && annotate == 'no' ) {
        exit 1, "In data format is SV VCF, but you've chosen not to annotate. Nothing for pipeface to do."
    }
    if ( in_data_format == 'snv_vcf' && calculate_depth == 'yes' ) {
        exit 1, "In data format is SNP/indel VCF, but you've chosen to calculate depth (which requires a bam file). Either set calculate_depth to 'no' in parameter file or pass '--calculate_depth no' on the command line"
    }
    if ( in_data_format == 'sv_vcf' && calculate_depth == 'yes' ) {
        exit 1, "In data format is SV VCF, but you've chosen to calculate depth (which requires a bam file). Either set calculate_depth to 'no' in parameter file or pass '--calculate_depth no' on the command line"
    }
    if ( !ref ) {
        exit 1, "No reference genome provided. Either include in parameter file or pass to --ref on the command line."
    }
    if ( !file(ref).exists() ) {
        exit 1, "Reference genome file path does not exist, '${ref}' provided."
    }
    if ( !ref_index ) {
        exit 1, "No reference genome index provided. Either include in parameter file or pass to --ref_index on the command line."
    }
    if ( !file(ref_index).exists() ) {
        exit 1, "Reference genome index file path does not exist, '${ref_index}' provided."
    }
    if ( !tandem_repeat ) {
        exit 1, "No tandem repeat bed file provided. Either include in parameter file or pass to --tandem_repeat on the command line. Set to 'NONE' if you do not wish to use a tandem repeat bed file."
    }
    if ( !file(tandem_repeat).exists() ) {
        exit 1, "Tandem repeat bed file path does not exist, '${tandem_repeat}' provided."
    }
    if ( !snp_indel_caller ) {
        exit 1, "No SNP/indel calling software selected. Either include in parameter file or pass to --snp_indel_caller on the command line. Should be either 'clair3' or 'deepvariant'."
    }
    if ( in_data_format != 'snv_vcf' && in_data_format != 'sv_vcf' && snp_indel_caller != 'clair3' && snp_indel_caller != 'deepvariant' ) {
        exit 1, "SNP/indel calling software should be either 'clair3' or 'deepvariant', '${snp_indel_caller}' selected."
    }
    if ( !sv_caller ) {
        exit 1, "No SV calling software selected. Either include in parameter file or pass to --sv_caller on the command line. Should be 'sniffles', 'cutesv', or 'both'."
    }
    if ( in_data_format != 'snv_vcf' && in_data_format != 'sv_vcf' && sv_caller != 'sniffles' && sv_caller != 'cutesv' && sv_caller != 'both' ) {
        exit 1, "SV calling software should be 'sniffles', 'cutesv', or 'both', '${sv_caller}' selected."
    }
    if ( in_data_format == 'sv_vcf' && sv_caller != 'sniffles' && sv_caller != 'cutesv' ) {
        exit 1, "When the input data format is 'sv_vcf', SV calling software should be 'sniffles' or 'cutesv', '${sv_caller}' selected."
    }
    if ( in_data_format == 'snv_vcf' && ref == 'NONE' ) {
        exit 1, "When the input data format is 'snv_vcf', please pass the reference genome used to generate the input data to 'ref' instead of setting it to 'NONE'."
    }
    if ( in_data_format == 'snv_vcf' && ref_index == 'NONE' ) {
        exit 1, "When the input data format is 'snv_vcf', please pass the reference genome index used to generate the input data to 'ref_index' instead of setting it to 'NONE'."
    }
    if ( in_data_format == 'snv_vcf' && snp_indel_caller == 'NONE' ) {
        exit 1, "When the input data format is 'snv_vcf', please pass the SNP/indel calling software used to generate the input data to 'snp_indel_caller' instead of setting it to 'NONE'."
    }
    if ( in_data_format == 'sv_vcf' && sv_caller == 'NONE' ) {
        exit 1, "When the input data format is 'sv_vcf', please pass the SV calling software used to generate the input data to 'sv_caller' instead of setting it to 'NONE'."
    }
    if ( !annotate ) {
        exit 1, "Choice to annotate not made. Either include in parameter file or pass to --annotate on the command line. Should be either 'yes' or 'no'."
    }
    if ( annotate != 'yes' && annotate != 'no' ) {
        exit 1, "Choice to annotate should be either 'yes', or 'no', '${annotate}' selected."
    }
    if ( annotate == 'yes' && !ref.toLowerCase().contains('hg38') && !ref.toLowerCase().contains('grch38') && annotate_override != 'yes' ) {
        exit 1, "Annotation of only hg38/GRCh38 is supported. You've chosen to annotate, but it looks like you may not be passing a hg38/GRCh38 reference genome based on the filename of the reference genome. '${ref}' passed. Pass '--annotate_override yes' on the command line to override this error."
    }
    if ( annotate == 'yes' && !file(vep_db).exists() ) {
        exit 1, "VEP cache directory does not exist, '${vep_db}' provided."
    }
    if ( annotate == 'yes' && !file(revel_db).exists() ) {
        exit 1, "REVEL database file path does not exist, '${revel_db}' provided."
    }
    if ( annotate == 'yes' && !file(gnomad_db).exists() ) {
        exit 1, "gnomAD database file path does not exist, '${gnomad_db}' provided."
    }
    if ( annotate == 'yes' && !file(gnomad_sv_db).exists() ) {
        exit 1, "nomAD SV database file path does not exist, '${gnomad_sv_db}' provided."
    }
    if ( annotate == 'yes' && !file(clinvar_db).exists() ) {
        exit 1, "ClinVar database file path does not exist, '${clinvar_db}' provided."
    }
    if ( annotate == 'yes' && !file(cadd_snv_db).exists() ) {
        exit 1, "CADD SNV database file path does not exist, '${cadd_snv_db}' provided."
    }
    if ( annotate == 'yes' && !file(cadd_indel_db).exists() ) {
        exit 1, "CADD indel database file path does not exist, '${cadd_indel_db}' provided."
    }
    if ( annotate == 'yes' && !file(cadd_sv_db).exists() ) {
        exit 1, "CADD SV database file path does not exist, '${cadd_sv_db}' provided."
    }
    if ( annotate == 'yes' && !file(spliceai_snv_db).exists() ) {
        exit 1, "SpliceAI SNV database file path does not exist, '${spliceai_snv_db}' provided."
    }
    if ( annotate == 'yes' && !file(spliceai_indel_db).exists() ) {
        exit 1, "SpliceAI indel database file path does not exist, '${spliceai_indel_db}' provided."
    }
    if ( annotate == 'yes' && !file(alphamissense_db).exists() ) {
        exit 1, "AlphaMissense database file path does not exist, '${alphamissense_db}' provided."
    }
    if ( !calculate_depth ) {
        exit 1, "Choice to calculate depth not made. Either include in parameter file or pass to --calculate_depth on the command line. Should be either 'yes' or 'no'."
    }
    if ( calculate_depth != 'yes' && calculate_depth != 'no' ) {
        exit 1, "Choice to calculate depth should be either 'yes', or 'no', '${calculate_depth}' selected."
    }
    if ( !outdir ) {
        exit 1, "No output directory provided. Either include in parameter file or pass to --outdir on the command line."
    }
    if ( !mosdepth_binary ) {
        exit 1, "No mosdepth binary provided. Either include in parameter file or pass to --mosdepth_binary on the command line. Set to 'NONE' if not running depth calculation."
    }
    if ( mosdepth_binary != 'NONE' && calculate_depth == 'no') {
        exit 1, "Pass 'NONE' to 'mosdepth_binary' when choosing to NOT calculate depth, '${mosdepth_binary}' and '${calculate_depth}' respectively provided'."
    }
    if ( mosdepth_binary == 'NONE' && calculate_depth == 'yes') {
        exit 1, "Pass an appropriate path to 'mosdepth_binary' when choosing to calculate depth, '${mosdepth_binary}' and '${calculate_depth}' respectively provided'."
    }
    if ( !pbcpgtools_binary ) {
        exit 1, "No pb-CpG-tools binary provided. Either include in parameter file or pass to --pbcpgtools_binary on the command line. Set to 'NONE' if not analysing any pacbio data."
    }
    if ( in_data_format == 'snv_vcf' && pbcpgtools_binary != 'NONE' ) {
        exit 1, "When the input data format is 'snv_vcf', please set the pb-CpG-tools binary (pbcpgtools_binary) to 'NONE'."
    }
    if ( in_data_format == 'sv_vcf' && pbcpgtools_binary != 'NONE' ) {
        exit 1, "When the input data format is 'sv_vcf', please set the pb-CpG-tools binary (pbcpgtools_binary) to 'NONE'."
    }
    if ( !file(in_data).exists() ) {
        exit 1, "In data csv file path does not exist, '${in_data}' provided."
    }
    if ( !file(ref).exists() ) {
        exit 1, "Reference genome file path does not exist, '${ref}' provided."
    }
    if ( !file(ref_index).exists() ) {
        exit 1, "Reference genome index file path does not exist, '${ref_index}' provided."
    }
    if ( !file(tandem_repeat).exists() ) {
        exit 1, "Tandem repeat bed file path does not exist, '${tandem_repeat}' provided."
    }
    if ( !file(mosdepth_binary).exists() ) {
        exit 1, "mosdepth binary file path does not exist, '${mosdepth_binary}' provided."
    }

    // build variable
    ref_name = file(ref).getSimpleName()

    // build channel from in_data.csv file for main workflow
    // groupTuple will collapse by sample_id (as defined in the in_data.csv file), creating a list of files per sample_id
    Channel
        .fromPath( in_data )
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple( row.sample_id, row.family_id, file(row.file).getExtension(), row.file, row.data_type, row.regions_of_interest, row.clair3_model ) }
        .groupTuple(by: [0,1,2,4,5,6] )
        .set { in_data_tuple }

    // build a list of files NOT collaped by sample_id (as defined in the in_data.csv file) for reporting
    Channel
	.fromPath( in_data )
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple( row.sample_id, row.family_id, file(row.file).getExtension(), row.file,file(row.file).getName() ) }
        .set { in_data_list }

    // build channels from in_data.csv file for describing each metadata associated with samples/families
    // I can use these downstream to join to input tuples throughout the pipeline as needed to avoid needing to have all sample metadata in all process input/output tuples
    Channel
	.fromPath( in_data )
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple( row.sample_id, row.family_id ) }
        .groupTuple(by: [0,1] )
        .set { id_tuple }

    Channel
        .fromPath( in_data )
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple( row.sample_id, row.family_id, file(row.file).getExtension() ) }
        .groupTuple(by: [0,1,2] )
        .set { extension_tuple }

    Channel
	.fromPath( in_data )
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple( row.sample_id, row.family_id, row.file ) }
        .groupTuple(by: [0,1] )
        .set { files_tuple }

    Channel
	.fromPath( in_data )
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple( row.sample_id, row.family_id, row.data_type ) }
        .groupTuple(by: [0,1,2] )
        .set { data_type_tuple }

    Channel
	.fromPath( in_data )
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple( row.sample_id, row.family_id, row.regions_of_interest ) }
        .groupTuple(by: [0,1,2] )
        .set { regions_of_interest_tuple }

    Channel
        .fromPath( in_data )
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple( row.sample_id, row.family_id, row.clair3_model ) }
        .groupTuple(by: [0,1,2] )
	.set { clair3_model_tuple }

    // build channel from in_data.csv file for user input checks
    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            def sample_id = row.sample_id
            def family_id = row.family_id
            def in_file = row.file
            def data_type = row.data_type
            def regions_of_interest = row.regions_of_interest
            def clair3_model = row.clair3_model

    // check user provided parameters in in_data.csv file
    if ( sample_id.isEmpty() ) {
       exit 1, "There is an empty entry in the 'sample_id' column of '$in_data'."
    }
    if ( in_file.isEmpty() ) {
       exit 1, "There is an empty entry in the 'file' column of '$in_data'."
    }
    if ( !file(in_file).exists() ) {
       exit 1, "There is an entry in the 'file' column of '$in_data' which doesn't exist. Check file '$in_file'."
    }
    if ( file(in_file).getExtension() != 'bam' && file(in_file).getExtension() != 'gz' && file(in_file).getExtension() != 'fastq' ) {
       exit 1, "There is an entry in the 'file' column of '$in_data' which doesn't have a 'bam', 'gz' or 'fastq' file extension. '$in_file' provided."
    }
    if ( in_data_format == 'aligned_bam' && !file(in_file + '.bai').exists() ) {
       exit 1, "You've specified that the in data format is aligned BAM, but it looks '${in_file}' defined in the file column of '$in_data' isn't indexed (ie. '${in_file}.bai doesn't exist')."
    }
    if ( in_data_format == 'aligned_bam' && !in_file.contains('bam') && in_data_format_override != 'yes' ) {
       exit 1, "You've specified that the in data format is aligned BAM, but it looks like you may not be passing a BAM file based on the file names defined in the file column of '$in_data'. ${in_file} passed. Pass '--in_data_format_override yes' on the command line to override this error."
    }
    if ( in_data_format == 'snv_vcf' && !in_file.contains('vcf') && in_data_format_override != 'yes' ) {
       exit 1, "You've specified that the in data format is SNV vcf, but it looks like you may not be passing a VCF file based on the file names defined in the file column of '$in_data'. ${in_file} passed. Pass '--in_data_format_override yes' on the command line to override this error."
    }
    if ( in_data_format == 'sv_vcf' && !in_file.contains('vcf') && in_data_format_override != 'yes' ) {
       exit 1, "You've specified that the in data format is SV vcf, but it looks like you may not be passing a VCF file based on the file names defined in the file column of '$in_data'. ${in_file} passed. Pass '--in_data_format_override yes' on the command line to override this error."
    }
    if ( data_type.isEmpty() ) {
       exit 1, "There is an empty entry in the 'data_type' column of '$in_data'."
    }
    if ( in_data_format != 'snv_vcf' && in_data_format != 'sv_vcf' && data_type != 'ont' && data_type != 'pacbio' ) {
       exit 1, "There is an entry in the 'data_type' column of '$in_data' that is not 'ont' or 'pacbio', '$data_type' provided."
    }
    if ( in_data_format == 'snv_vcf' && data_type != 'NONE' ) {
       exit 1, "When the input data format is 'snv_vcf', please set the data type (data_type column of '$in_data') to 'NONE'."
    }
    if ( in_data_format == 'sv_vcf' && data_type != 'NONE' ) {
       exit 1, "When the input data format is 'sv_vcf', please set the data type (data_type column of '$in_data') to 'NONE'."
    }
    if ( regions_of_interest.isEmpty() ) {
       exit 1, "There is an empty entry in the 'regions_of_interest' column of '$in_data'."
    }
    if ( !file(regions_of_interest).exists() ) {
       exit 1, "There is an entry in the 'regions_of_interest' column of '$in_data' which doesn't exist. Check file '$regions_of_interest'."
    }
    if ( clair3_model.isEmpty() ) {
       exit 1, "There is an empty entry in the 'clair3_model' column of '$in_data'."
    }
    if ( !file(clair3_model).exists() ) {
       exit 1, "There is an entry in the 'clair3_model' column of '$in_data' which doesn't exist. Check path '$clair3_model'."
    }
    if ( snp_indel_caller != 'clair3' && clair3_model != 'NONE' ) {
       exit 1, "Pass 'NONE' in the 'clair3_model' column of '$in_data' when clair3 is NOT selected as the SNP/indel calling software, '$clair3_model' provided'."
    }
    if ( in_data_format != 'snv_vcf' && in_data_format != 'sv_vcf' && snp_indel_caller == 'clair3' && clair3_model == 'NONE' ) {
       exit 1, "When clair3 is selected as the SNP/indel calling software, provide a path to an appropriate clair3 model in the 'clair3_model' column of '$in_data' rather than setting it to 'NONE'."
    }
    if ( in_data_format == 'snv_vcf' && clair3_model != 'NONE' ) {
       exit 1, "When the input data format is 'snv_vcf', please set the Clair3 model (clair3_model column of '$in_data') to 'NONE'."
    }
    if ( in_data_format == 'sv_vcf' && clair3_model != 'NONE' ) {
       exit 1, "When the input data format is 'sv_vcf', please set the Clair3 model (clair3_model column of '$in_data') to 'NONE'."
    }
    }

    // workflow
    // pre-process, alignment and qc
    scrape_settings(in_data_tuple, in_data, in_data_format, ref, ref_index, tandem_repeat, snp_indel_caller, sv_caller, annotate, calculate_depth, outdir, outdir2)
    if ( in_data_format == 'ubam_fastq' | in_data_format == 'aligned_bam' ) {
        bam_header = scrape_bam_header(in_data_list, outdir, outdir2)
    }
    if ( in_data_format == 'ubam_fastq' ) {
        merged = merge_runs(id_tuple.join(extension_tuple, by: [0,1]).join(files_tuple, by: [0,1]))
        bam = minimap2(merged.join(extension_tuple, by: [0,1]).join(data_type_tuple, by: [0,1]), ref, ref_index)
    }
    if ( in_data_format == 'aligned_bam' ) {
        bam = id_tuple.join(files_tuple, by: [0,1])
    }
    if ( in_data_format == 'ubam_fastq' | in_data_format == 'aligned_bam' ) {
        if ( calculate_depth == 'yes' ) {
            mosdepth(bam.join(regions_of_interest_tuple, by: [0,1]), mosdepth_binary, outdir, outdir2, ref_name)
        }
        // snp/indel calling
        if ( snp_indel_caller == 'clair3' ) {
            snp_indel_vcf_bam = clair3(bam.join(data_type_tuple, by: [0,1]).join(regions_of_interest_tuple, by: [0,1]).join(clair3_model_tuple, by: [0,1]), ref, ref_index)
        }
        else if ( snp_indel_caller == 'deepvariant' ) {
            deepvariant_dry_run(bam.join(data_type_tuple, by: [0,1]), ref, ref_index)
            deepvariant_make_examples(deepvariant_dry_run.out.join(regions_of_interest_tuple, by: [0,1]), ref, ref_index)
            deepvariant_call_variants(deepvariant_make_examples.out)
            snp_indel_vcf_bam = deepvariant_post_processing(deepvariant_call_variants.out, ref, ref_index)
        }
        // phasing
        (snp_indel_phased_vcf_bam, snp_indel_phased_vcf, phased_read_list) = whatshap_phase(snp_indel_vcf_bam, ref, ref_index, outdir, outdir2, ref_name, snp_indel_caller)
        // haplotagging
        (haplotagged_bam, haplotagged_bam_fam, haplotagged_tsv) = whatshap_haplotag(snp_indel_phased_vcf_bam, ref, ref_index, outdir, outdir2, ref_name)
        // methylation analysis
        pbcpgtools(haplotagged_bam.join(data_type_tuple, by: [0,1]), pbcpgtools_binary, ref, ref_index, outdir, outdir2, ref_name)
        // sv calling
        if ( sv_caller == 'sniffles' | sv_caller == 'both' ) {
            (sv_vcf_sniffles, sv_vcf_sniffles_indexed) = sniffles(haplotagged_bam, ref, ref_index, tandem_repeat, outdir, outdir2, ref_name)
        }
        if ( sv_caller == 'cutesv' | sv_caller == 'both' ) {
            (sv_vcf_cutesv, sv_vcf_cutesv_indexed) = cutesv(haplotagged_bam.join(data_type_tuple, by: [0,1]), ref, ref_index, tandem_repeat, outdir, outdir2, ref_name)
        }
    }
    if ( in_data_format == 'snv_vcf' ) {
        snp_indel_phased_vcf = id_tuple.join(files_tuple, by: [0,1])
    }
    if ( in_data_format == 'ubam_fastq' | in_data_format == 'aligned_bam' | in_data_format == 'snv_vcf' ) {
        // annotation
        if ( annotate == 'yes' ) {
            vep_snv(snp_indel_phased_vcf, ref, ref_index, vep_db, revel_db, gnomad_db, clinvar_db, cadd_snv_db, cadd_indel_db, spliceai_snv_db, spliceai_indel_db, alphamissense_db, outdir, outdir2, ref_name, snp_indel_caller)
        }
    }
    if ( in_data_format == 'sv_vcf' ) {
        sv_vcf_sniffles = id_tuple.join(files_tuple, by: [0,1])
        sv_vcf_cutesv = id_tuple.join(files_tuple, by: [0,1])
    }
    if ( in_data_format == 'ubam_fastq' | in_data_format == 'aligned_bam' | in_data_format == 'sv_vcf' ) {
        if ( annotate == 'yes' ) {
            if ( sv_caller == 'sniffles' | sv_caller == 'both' ) {
                vep_sniffles_sv(sv_vcf_sniffles, ref, ref_index, vep_db, gnomad_db, gnomad_sv_db, clinvar_db, cadd_sv_db, outdir, outdir2, ref_name)
            }
            if ( sv_caller == 'cutesv' | sv_caller == 'both' ) {
                vep_cutesv_sv(sv_vcf_cutesv, ref, ref_index, vep_db, gnomad_db, gnomad_sv_db, clinvar_db, cadd_sv_db, outdir, outdir2, ref_name)
            }
        }
    }
}

