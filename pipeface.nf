
nextflow.enable.dsl=2

// set defaults for optional params (set default optional input files to dummy NONE file)
// secret sauce second outdir
params.outdir2 = ""
params.tandem_repeat = "NONE"
params.annotate_override = ""

process scrape_settings {

    input:
        tuple val(sample_id), val(extension), val(files), val(data_type), val(regions_of_interest), val(clair3_model)
        val in_data
        val ref
        val ref_index
        val tandem_repeat
        val snp_indel_caller
        val sv_caller
        val annotate
        val outdir

    output:
        tuple val(sample_id), path('pipeface_settings.txt')

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
        """
        echo "Sample ID: $sample_id" >> pipeface_settings.txt
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
        echo "Outdir: $outdir" >> pipeface_settings.txt
        """

    stub:
        """
        touch pipeface_settings.txt
        """

}

process publish_settings {

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$filename" }

    input:
        tuple val(sample_id), path(pipeface_settings)
        val outdir
        val outdir2
        val ref_name

    output:
        path 'pipeface_settings.txt'

    script:
        """
        OUT_PATH="$outdir/$sample_id/$outdir2"
        # if a pipeface_settings.txt file exists ...
        if [ -f \${OUT_PATH}/pipeface_settings.txt ]; then
            # ... add +1 to the suffix of each file that has a suffix
            for FILE in `ls -1vr \${OUT_PATH}/pipeface_settings.*.txt`; do
                # find the current suffix number by stripping the filepath and returning the only the number
                SUFFIX=\$(echo \${FILE} | sed 's/.*pipeface_settings//g' | tr -cd '[:digit:]')
                NEW_SUFFIX=\$((\${SUFFIX} + 1 ))
                # update suffix
                mv \${FILE} \${OUT_PATH}/pipeface_settings.\${NEW_SUFFIX}.txt
            done
            # ... and start suffix for pipeface_settings.txt
            mv \${OUT_PATH}/pipeface_settings.txt \${OUT_PATH}/pipeface_settings.1.txt
        fi
        echo "Publishing files"
        """

    stub:
        """
        touch pipeface_settings.txt
        """

}

process scrape_bam_header {

    input:
        tuple val(sample_id), val(extension), val(file)

    output:
        tuple val(sample_id), val(file), path('header')

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

process publish_bam_header {

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$file.$filename" }

    input:
        tuple val(sample_id), path(file), path(header)
        val outdir
        val outdir2

    output:
        path 'header'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch header
        """

}

process merge_runs {

    input:
        tuple val(sample_id), val(extension), path(files), val(data_type), val(regions_of_interest), val(clair3_model)

    output:
        tuple val(sample_id), val(extension), val(data_type), val(regions_of_interest), val(clair3_model), path('merged')

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
        tuple val(sample_id), val(extension), val(data_type), val(regions_of_interest), val(clair3_model), path(merged)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(data_type), val(regions_of_interest), val(clair3_model), path('sorted.bam'), path('sorted.bam.bai')
        tuple val(sample_id), path('sorted.bam')

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

process depth {

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path('depth.tsv')

    script:
        """
        # calculate average depth per chromosome
        samtools depth \
        -@ ${task.cpus} \
        -a \
        $bam | awk '{cov[\$1]+=\$3; len[\$1]++} END {for (chr in cov) print chr, "	"cov[chr]/len[chr]}' >> tmp.tsv
        # write tsv header
        echo "chromosome\taverage_depth" > depth.tsv
        # sort file
        cat tmp.tsv | sort -V >> depth.tsv
        # cleanup
        rm tmp.tsv
        """

    stub:
        """
        touch depth.tsv
        """

}

process publish_depth {

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$filename" }

    input:
        tuple val(sample_id), path(depth)
        val outdir
        val outdir2
        val ref_name

    output:
        path 'depth.tsv'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch depth.tsv
        """

}

process clair3 {

    input:
        tuple val(sample_id), val(data_type), val(regions_of_interest), val(clair3_model), path(bam), path(bam_index)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(data_type), path(bam), path(bam_index), path('snp_indel.vcf.gz'), path('snp_indel.vcf.gz.tbi')

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

process deepvariant {

    input:
        tuple val(sample_id), val(data_type), val(regions_of_interest), val(clair3_model), path(bam), path(bam_index)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(data_type), path(bam), path(bam_index), path('snp_indel.vcf.gz'), path('snp_indel.vcf.gz.tbi')

    script:
    // conditionally define model type
    if( data_type == 'ont' ) {
        model = 'ONT_R104'
    }
    else if ( data_type == 'pacbio' ) {
        model = 'PACBIO'
    }
    // define an optional string to pass regions of interest bed file
    def regions_of_interest_optional = file(regions_of_interest).name != 'NONE' ? "--regions $regions_of_interest" : ''
        """
        # run deepvariant
        singularity run /g/data/ox63/install/deepvariant_1.6.1-gpu.sif run_deepvariant \
        --reads=$bam \
        --ref=$ref \
        --sample_name=$sample_id \
        --output_vcf=snp_indel.vcf.gz \
        --model_type=$model \
        $regions_of_interest_optional \
        --num_shards=${task.cpus} \
        --postprocess_cpus=${task.cpus}
        """

    stub:
        """
        touch snp_indel.vcf.gz
        touch snp_indel.vcf.gz.tbi
        """

}

process whatshap_phase {

    input:
        tuple val(sample_id), val(data_type), path(bam), path(bam_index), path(snp_indel_vcf), path(snp_indel_vcf_index)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(data_type), path(bam), path(bam_index), path('snp_indel.phased.vcf.gz'), path('snp_indel.phased.vcf.gz.tbi')
        tuple val(sample_id), path('snp_indel.phased.vcf.gz'), path('snp_indel.phased.vcf.gz.tbi'), path('snp_indel.phased.read_list.txt')
        tuple val(sample_id), path('snp_indel.phased.vcf.gz'), path('snp_indel.phased.vcf.gz.tbi')

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
        """

    stub:
        """
        touch snp_indel.phased.vcf.gz
        touch snp_indel.phased.vcf.gz.tbi
        touch snp_indel.phased.read_list.txt
        """

}

process publish_whatshap_phase {

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$snp_indel_caller.$filename" }

    input:
        tuple val(sample_id), path(snp_indel_phased_vcf), path(snp_indel_phased_vcf_index), path(snp_indel_phased_read_list)
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller

    output:
        val sample_id
        path 'snp_indel.phased.vcf.gz'
        path 'snp_indel.phased.vcf.gz.tbi'
        path 'snp_indel.phased.read_list.txt'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch snp_indel.phased.vcf.gz
        touch snp_indel.phased.vcf.gz.tbi
        touch snp_indel.phased.read_list.txt
        """

}

process vep_snv {

    input:
        tuple val(sample_id), path(snp_indel_phased_vcf), path(snp_indel_phased_vcf_index)
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

    output:
        tuple val(sample_id), path('snp_indel.phased.annotated.vcf.gz'), path('snp_indel.phased.annotated.vcf.gz.tbi')

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

process publish_vep_snv {

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$snp_indel_caller.$filename" }

    input:
        tuple val(sample_id), path(snp_indel_phased_annotated_vcf), path(snp_indel_phased_annotated_vcf_index)
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller

    output:
        path 'snp_indel.phased.annotated.vcf.gz'
        path 'snp_indel.phased.annotated.vcf.gz.tbi'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch snp_indel.phased.annotated.vcf.gz
        touch snp_indel.phased.annotated.vcf.gz.tbi
        """

}

process whatshap_haplotag {

    input:
        tuple val(sample_id), val(data_type), path(bam), path(bam_index), path(snp_indel_vcf), path(snp_indel_vcf_index)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(data_type), path('sorted.haplotagged.bam'), path('sorted.haplotagged.bam.bai')
        tuple val(sample_id), path('sorted.haplotagged.bam'), path('sorted.haplotagged.bam.bai'), path('sorted.haplotagged.tsv')

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

process publish_whatshap_haplotag {

    def mapper_phaser = "minimap2.whatshap"

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$mapper_phaser.$filename" }

    input:
        tuple val(sample_id), path(haplotagged_bam), path(haplotagged_bam_index), path(haplotagged_bam_tsv)
        val outdir
        val outdir2
        val ref_name

    output:
        path 'sorted.haplotagged.bam'
        path 'sorted.haplotagged.bam.bai'
        path 'sorted.haplotagged.tsv'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch sorted.haplotagged.bam
        touch sorted.haplotagged.bam.bai
        touch sorted.haplotagged.tsv
        """

}

process sniffles {

    input:
        tuple val(sample_id), val(data_type), path(haplotagged_bam), path(haplotagged_bam_index)
        val ref
        val ref_index
        val tandem_repeat

    output:
        tuple val(sample_id), path('sv.phased.vcf.gz'), path('sv.phased.vcf.gz.tbi')

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
        --minsvlen 20 \
        --phase $tandem_repeat_optional
        """

    stub:
        """
        touch sv.phased.vcf.gz
        touch sv.phased.vcf.gz.tbi
        """

}

process publish_sniffles {

    def sv_caller = "sniffles"

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$sv_caller.$filename" }

    input:
        tuple val(sample_id), path(sv_vcf), path(sv_vcf_index)
        val outdir
        val outdir2
        val ref_name

    output:
        path 'sv.phased.vcf.gz'
        path 'sv.phased.vcf.gz.tbi'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch sv.phased.vcf.gz
        touch sv.phased.vcf.gz.tbi
        """

}

process cutesv {

    input:
        tuple val(sample_id), val(data_type), path(haplotagged_bam), path(haplotagged_bam_index)
        val ref
        val ref_index
        val tandem_repeat

    output:
        tuple val(sample_id), path('sv.vcf.gz'), path('sv.vcf.gz.tbi')

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

process publish_cutesv {

    def sv_caller = "cutesv"

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$sv_caller.$filename" }

    input:
        tuple val(sample_id), path(sv_vcf), path(sv_vcf_index)
        val outdir
        val outdir2
        val ref_name

    output:
        path 'sv.vcf.gz'
        path 'sv.vcf.gz.tbi'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch sv.vcf.gz
        touch sv.vcf.gz.tbi
        """

}

process vep_sniffles_sv {

    input:
        tuple val(sample_id), path(sv_phased_vcf), path(sv_phased_vcf_index)
        val ref
        val ref_index
        val sv_caller
        val vep_db
        val cadd_sv_db

    output:
        tuple val(sample_id), path('sv.phased.annotated.vcf.gz'), path('sv.phased.annotated.vcf.gz.tbi')

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

    input:
        tuple val(sample_id), path(sv_phased_vcf), path(sv_phased_vcf_index)
        val ref
        val ref_index
        val sv_caller
        val vep_db
        val cadd_sv_db

    output:
        tuple val(sample_id), path('sv.annotated.vcf.gz'), path('sv.annotated.vcf.gz.tbi')

    script:
        """
        # run vep
        vep -i $sv_phased_vcf \
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

process publish_vep_sniffles_sv {

    def sv_caller = "sniffles"

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$sv_caller.$filename" }

    input:
        tuple val(sample_id), path(sv_vcf), path(sv_vcf_index)
        val outdir
        val outdir2
        val ref_name

    output:
        path 'sv.phased.annotated.vcf.gz'
        path 'sv.phased.annotated.vcf.gz.tbi'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch sv.phased.annotated.vcf.gz
        touch sv.phased.annotated.vcf.gz.tbi
        """

}

process publish_vep_cutesv_sv {

    def sv_caller = "cutesv"

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$sv_caller.$filename" }

    input:
        tuple val(sample_id), path(sv_vcf), path(sv_vcf_index)
        val outdir
        val outdir2
        val ref_name

    output:
        path 'sv.annotated.vcf.gz'
        path 'sv.annotated.vcf.gz.tbi'

    script:
        """
        echo "Publishing files"
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
    ref = "$params.ref"
    ref_index = "$params.ref_index"
    tandem_repeat = "$params.tandem_repeat"
    snp_indel_caller = "$params.snp_indel_caller"
    sv_caller = "$params.sv_caller"
    annotate = "$params.annotate"
    annotate_override = "$params.annotate_override"
    outdir = "$params.outdir"
    outdir2 = "$params.outdir2"
    vep_db = "$params.vep_db"
    revel_db = "$params.revel_db"
    gnomad_db = "$params.gnomad_db"
    clinvar_db = "$params.clinvar_db"
    cadd_snv_db = "$params.cadd_snv_db"
    cadd_indel_db = "$params.cadd_indel_db"
    cadd_sv_db = "$params.cadd_sv_db"
    spliceai_snv_db = "$params.spliceai_snv_db"
    spliceai_indel_db = "$params.spliceai_indel_db"
    alphamissense_db = "$params.alphamissense_db"

    // check user provided parameters
    if ( !in_data ) {
        exit 1, "Error: No in data csv file specified. Either include in parameter file or pass to --in_data on the command line."
    }
    if ( !ref ) {
        exit 1, "Error: No reference genome provided. Either include in parameter file or pass to --ref on the command line."
    }
    if ( !ref_index ) {
        exit 1, "Error: No reference genome index provided. Either include in parameter file or pass to --ref_index on the command line."
    }
    if ( !tandem_repeat ) {
        exit 1, "Error: No tandem repeat bed file provided. Either include in parameter file or pass to --tandem_repeat on the command line. Set to 'NONE' if you do not wish to use a tandem repeat bed file."
    }
    if ( !snp_indel_caller ) {
        exit 1, "Error: No SNP/indel calling software selected. Either include in parameter file or pass to --snp_indel_caller on the command line. Should be either 'clair3' or 'deepvariant'."
    }
    if ( snp_indel_caller != 'clair3' && snp_indel_caller != 'deepvariant' ) {
        exit 1, "Error: SNP/indel calling software should be either 'clair3' or 'deepvariant', '${snp_indel_caller}' selected."
    }
    if ( !sv_caller ) {
        exit 1, "Error: No SV calling software selected. Either include in parameter file or pass to --sv_caller on the command line. Should be 'sniffles', 'cutesv', or 'both'."
    }
    if ( sv_caller != 'sniffles' && sv_caller != 'cutesv' && sv_caller != 'both' ) {
        exit 1, "Error: SV calling software should be 'sniffles', 'cutesv', or 'both', '${sv_caller}' selected."
    }
    if ( !annotate ) {
        exit 1, "Error: Choice to annotate not made. Either include in parameter file or pass to --annotate on the command line. Should be either 'yes' or 'no'."
    }
    if ( annotate != 'yes' && annotate != 'no' ) {
        exit 1, "Error: Choice to annotate should be either 'yes', or 'no', '${annotate}' selected."
    }
    if ( annotate == 'yes' && !ref.toLowerCase().contains('hg38') && !ref.toLowerCase().contains('grch38') && annotate_override != 'yes' ) {
        exit 1, "Warning: Annotation of only hg38/GRCh38 is supported. You've chosen to annotate, but it looks like you may not be passing a hg38/GRCh38 reference genome based on the filename of the reference genome. '${ref}' passed. Pass '--annotate_override yes' on the command line to override this error."
    }
    if ( annotate == 'yes' && !file(vep_db).exists() ) {
        exit 1, "Error VEP cache directory does not exist, '${vep_db}' provided."
    }
    if ( annotate == 'yes' && !file(revel_db).exists() ) {
        exit 1, "Error REVEL database file path does not exist, '${revel_db}' provided."
    }
    if ( annotate == 'yes' && !file(gnomad_db).exists() ) {
        exit 1, "Error gnomAD database file path does not exist, '${gnomad_db}' provided."
    }
    if ( annotate == 'yes' && !file(clinvar_db).exists() ) {
        exit 1, "Error ClinVar database file path does not exist, '${clinvar_db}' provided."
    }
    if ( annotate == 'yes' && !file(cadd_snv_db).exists() ) {
        exit 1, "Error CADD SNV database file path does not exist, '${cadd_snv_db}' provided."
    }
    if ( annotate == 'yes' && !file(cadd_indel_db).exists() ) {
        exit 1, "Error CADD indel database file path does not exist, '${cadd_indel_db}' provided."
    }
    if ( annotate == 'yes' && !file(cadd_sv_db).exists() ) {
        exit 1, "Error CADD SV database file path does not exist, '${cadd_snv_db}' provided."
    }
    if ( annotate == 'yes' && !file(spliceai_snv_db).exists() ) {
        exit 1, "Error SpliceAI SNV database file path does not exist, '${spliceai_snv_db}' provided."
    }
    if ( annotate == 'yes' && !file(spliceai_indel_db).exists() ) {
        exit 1, "Error SpliceAI indel database file path does not exist, '${spliceai_indel_db}' provided."
    }
    if ( annotate == 'yes' && !file(alphamissense_db).exists() ) {
        exit 1, "Error AlphaMissense database file path does not exist, '${alphamissense_db}' provided."
    }
    if ( !outdir ) {
        exit 1, "Error: No output directory provided. Either include in parameter file or pass to --outdir on the command line."
    }
    if ( !file(in_data).exists() ) {
        exit 1, "Error: In data csv file path does not exist, '${in_data}' provided."
    }
    if ( !file(ref).exists() ) {
        exit 1, "Error: Reference genome file path does not exist, '${ref}' provided."
    }
    if ( !file(ref_index).exists() ) {
        exit 1, "Error: Reference genome index file path does not exist, '${ref_index}' provided."
    }
    if ( !file(tandem_repeat).exists() ) {
        exit 1, "Error: Tandem repeat bed file path does not exist, '${tandem_repeat}' provided."
    }

    // build variable
    ref_name = file(ref).getSimpleName()

    // build a list of files NOT collaped by sample_id (as defined in the in_data.csv file) for reporting
    Channel
	.fromPath( in_data )
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple( row.sample_id, file(row.file).getExtension(), row.file ) }
        .set { in_data_list }

    // build channel from in_data.csv file for analysis
    // groupTuple will collapse by sample_id (as defined in the in_data.csv file), creating a list of files per sample_id
    Channel
        .fromPath( in_data )
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple( row.sample_id, file(row.file).getExtension(), row.file, row.data_type, row.regions_of_interest, row.clair3_model ) }
        .groupTuple(by: [0,1,3,4,5] )
        .set { in_data_tuple }

    // build channel from in_data.csv file for user input checks
    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            def sample_id = row.sample_id
            def in_file = row.file
            def data_type = row.data_type
            def regions_of_interest = row.regions_of_interest
            def clair3_model = row.clair3_model

    // check user provided parameters in in_data.csv file
    if ( sample_id.isEmpty() ) {
       exit 1, "Error processing '$in_data' file. There is an empty entry in the 'sample_id' column"
    }
    if ( in_file.isEmpty() ) {
       exit 1, "Error processing '$in_data' file. There is an empty entry in the 'file' column."
    }
    if ( data_type.isEmpty() ) {
       exit 1, "Error processing '$in_data' file. There is an empty entry in the 'data_type' column."
    }
    if ( regions_of_interest.isEmpty() ) {
       exit 1, "Error processing '$in_data' file. There is an empty entry in the 'regions_of_interest' column."
    }
    if ( clair3_model.isEmpty() ) {
       exit 1, "Error processing '$in_data' file. There is an empty entry in the 'clair3_model' column."
    }
    if ( data_type != 'ont' && data_type != 'pacbio' ) {
       exit 1, "Error processing '$in_data' file. There is an entry in the 'data_type' column that is not 'ont' or 'pacbio', '$data_type' provided."
    }
    if ( !file(in_file).exists() ) {
       exit 1, "Error processing '$in_data' file. There is an entry in the 'file' column which doesn't exist. Check file '$in_file'."
    }
    if ( !file(regions_of_interest).exists() ) {
       exit 1, "Error processing '$in_data' file. There is an entry in the 'regions_of_interest' column which doesn't exist. Check file '$regions_of_interest'."
    }
    if ( !file(clair3_model).exists() ) {
       exit 1, "Error processing '$in_data' file. There is an entry in the 'clair3_model' column which doesn't exist. Check path '$clair3_model'."
    }
    if ( file(in_file).getExtension() != "bam" && file(in_file).getExtension() != "gz" && file(in_file).getExtension() != "fastq" ) {
       exit 1, "Error processing '$in_data' file. There is an entry in the 'file' column which doesn't have a 'bam', 'gz' or 'fastq' file extension. '$in_file' provided."
    }
    if ( snp_indel_caller != "clair3" && clair3_model != "NONE" ) {
       exit 1, "Error processing '$in_data' file. Pass 'NONE' in the 'clair3_model' column when clair3 is NOT selected as the SNP/indel calling software, '$clair3_model' provided'."
    }
    if ( snp_indel_caller == "clair3" && clair3_model == "NONE" ) {
       exit 1, "Error processing '$in_data' file. When clair3 is selected as the SNP/indel calling software, provide a path to an appropriate clair3 model in the 'clair3_model' column rather than setting it to 'NONE'."
    }
    }

    // workflow
    scrape_settings_to_publish = scrape_settings(in_data_tuple, in_data, ref, ref_index, tandem_repeat, snp_indel_caller, sv_caller, annotate, outdir)
    publish_settings(scrape_settings_to_publish, outdir, outdir2, ref_name)
    bam_header = scrape_bam_header(in_data_list)
    publish_bam_header(bam_header, outdir, outdir2)
    merged = merge_runs(in_data_tuple)
    (bam, minimap_to_publish1) = minimap2(merged, ref, ref_index)
    minimap_to_publish2 = depth(minimap_to_publish1)
    publish_depth(minimap_to_publish2, outdir, outdir2, ref_name)
    if ( snp_indel_caller == 'clair3' ) {
        (snp_indel_vcf_bam, snp_indel_vcf) = clair3(bam, ref, ref_index)
    }
    else if ( snp_indel_caller == 'deepvariant' ) {
        (snp_indel_vcf_bam, snp_indel_vcf) = deepvariant(bam, ref, ref_index)
    }
    (snp_indel_phased_vcf_bam, whatshap_phase_to_publish, snp_indel_phased_vcf) = whatshap_phase(snp_indel_vcf_bam, ref, ref_index)
    publish_whatshap_phase(whatshap_phase_to_publish, outdir, outdir2, ref_name, snp_indel_caller)
    if ( annotate == 'yes' ) {
        vep_snv_to_publish = vep_snv(snp_indel_phased_vcf, ref, ref_index, vep_db, revel_db, gnomad_db, clinvar_db, cadd_snv_db, cadd_indel_db, spliceai_snv_db, spliceai_indel_db, alphamissense_db)
        publish_vep_snv(vep_snv_to_publish, outdir, outdir2, ref_name, snp_indel_caller)
    }
    (haplotagged_bam, whatshap_haplotag_to_publish) = whatshap_haplotag(snp_indel_phased_vcf_bam, ref, ref_index)
    publish_whatshap_haplotag(whatshap_haplotag_to_publish, outdir, outdir2, ref_name)
    if ( sv_caller == 'sniffles' ) {
        sv_vcf = sniffles(haplotagged_bam, ref, ref_index, tandem_repeat)
        publish_sniffles(sv_vcf, outdir, outdir2, ref_name)
        if ( annotate == 'yes' ) {
            annotated_sv_vcf_sniffles = vep_sniffles_sv(sv_vcf, ref, ref_index, sv_caller, vep_db, cadd_sv_db)
            publish_vep_sniffles_sv(annotated_sv_vcf_sniffles, outdir, outdir2, ref_name)
        }
    }
    else if ( sv_caller == 'cutesv' ) {
        sv_vcf = cutesv(haplotagged_bam, ref, ref_index, tandem_repeat)
        publish_cutesv(sv_vcf, outdir, outdir2, ref_name)
        if ( annotate == 'yes' ) {
            annotated_sv_vcf_cutesv = vep_cutesv_sv(sv_vcf, ref, ref_index, sv_caller, vep_db, cadd_sv_db)
            publish_vep_cutesv_sv(annotated_sv_vcf_cutesv, outdir, outdir2, ref_name)
        } 
    }
    else if ( sv_caller == 'both' ) {
        sv_vcf_sniffles = sniffles(haplotagged_bam, ref, ref_index, tandem_repeat)
        publish_sniffles(sv_vcf_sniffles, outdir, outdir2, ref_name)
        sv_vcf_cutesv = cutesv(haplotagged_bam, ref, ref_index, tandem_repeat)
        publish_cutesv(sv_vcf_cutesv, outdir, outdir2, ref_name)
        if ( annotate == 'yes' ) {
            annotated_sv_vcf_sniffles = vep_sniffles_sv(sv_vcf_sniffles, ref, ref_index, sv_caller, vep_db, cadd_sv_db)
            publish_vep_sniffles_sv(annotated_sv_vcf_sniffles, outdir, outdir2, ref_name)
            annotated_sv_vcf_cutesv = vep_cutesv_sv(sv_vcf_cutesv, ref, ref_index, sv_caller, vep_db, cadd_sv_db)
            publish_vep_cutesv_sv(annotated_sv_vcf_cutesv, outdir, outdir2, ref_name)
        }
    }
}

