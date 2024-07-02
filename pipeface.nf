nextflow.enable.dsl=2

// set defaults for optional params (set default optional input files to dummy NONE file)
// secret sauce second outdir
params.outdir2 = ""
params.tandem_repeat = "NONE"

process scrape_settings {

    input:
        tuple val(sample_id), val(extension), val(files), val(data_type), val(regions_of_interest), val(clair3_model)
        val in_data
        val ref
        val ref_index
        val tandem_repeat
        val snp_indel_caller
        val sv_caller
        val outdir

    output:
        tuple val(sample_id), path('pipeface_settings.txt')

    script:
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
        echo "SV caller: $sv_caller" >> pipeface_settings.txt
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
                "cat ${files} | bgzip > merged"
            } else {
                "ln -s ${files} merged"
            }
        }
        """
    else if( extension == 'bam' )
        """
        ${
            if (files instanceof List && files.size() > 1) { 
                "samtools merge ${files} -o merged"
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
        tuple val(sample_id), val(data_type), val(regions_of_interest), val(clair3_model), path('minimap2.sorted.bam'), path('minimap2.sorted.bam.bai')
        tuple val(sample_id), path('minimap2.version.txt')

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
        $ref - | samtools sort -@ ${task.cpus} -o minimap2.tmp.bam -
        # add sample id to bam header for downstream deepvariant
        samtools addreplacerg \
        -r ID:S1 \
        -r SM:$sample_id \
        -o minimap2.sorted.bam minimap2.tmp.bam
        # index bam
        samtools index \
        -@ ${task.cpus} \
        minimap2.sorted.bam
        # cleanup
        rm minimap2.tmp.bam
        # grab version
        minimap2 --version > minimap2.version.txt
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
        $merged | samtools sort -@ ${task.cpus} -o minimap2.tmp.bam -
        # add sample id to bam header for downstream deepvariant
        samtools addreplacerg \
        -r ID:S1 \
        -r SM:$sample_id \
        -o minimap2.sorted.bam minimap2.tmp.bam
        # index bam
        samtools index \
        -@ ${task.cpus} \
        minimap2.sorted.bam
        # cleanup
        rm minimap2.tmp.bam
        # grab version
        minimap2 --version > minimap2.version.txt
        """

    stub:
        """
        touch minimap2.sorted.bam
        touch minimap2.sorted.bam.bai
        touch minimap2.version.txt
        """

}

process publish_minimap2 {

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$filename" }

    input:
        tuple val(sample_id), path(minimap_version)
        val outdir
        val outdir2
        val ref_name

    output:
        path 'minimap2.version.txt'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch minimap2.version.txt
        """

}

process clair3 {

    input:
        tuple val(sample_id), val(data_type), val(regions_of_interest), val(clair3_model), path(bam), path(bam_index)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(data_type), path(bam), path(bam_index), path('clair3.snp_indel.vcf.gz'), path('clair3.snp_indel.vcf.gz.tbi')
        tuple val(sample_id), path('clair3.snp_indel.g.vcf.gz'), path('clair3.snp_indel.g.vcf.gz.tbi'), path('clair3.version.txt')

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
        $regions_of_interest_optional
        # rename files
        ln -s merge_output.vcf.gz clair3.snp_indel.vcf.gz
        ln -s merge_output.vcf.gz.tbi clair3.snp_indel.vcf.gz.tbi
        ln -s merge_output.gvcf.gz clair3.snp_indel.g.vcf.gz
        ln -s merge_output.gvcf.gz.tbi clair3.snp_indel.g.vcf.gz.tbi
        # grab version
        run_clair3.sh --version > clair3.version.txt
        """

    stub:
        """
        touch clair3.snp_indel.vcf.gz
        touch clair3.snp_indel.vcf.gz.tbi
        touch clair3.snp_indel.g.vcf.gz
        touch clair3.snp_indel.g.vcf.gz.tbi
        touch clair3.version.txt
        """

}

process publish_clair3 {

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$filename" }

    input:
        tuple val(sample_id), path(snp_indel_gvcf), path(snp_indel_gvcf_index), path(version)
        val outdir
        val outdir2
        val ref_name

    output:
        val sample_id
        path 'clair3.snp_indel.g.vcf.gz'
        path 'clair3.snp_indel.g.vcf.gz.tbi'
        path 'clair3.version.txt'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch clair3.snp_indel.g.vcf.gz
        touch clair3.snp_indel.g.vcf.gz.tbi
        touch clair3.version.txt
        """

}

process whatshap_phase_clair3 {

    input:
        tuple val(sample_id), val(data_type), path(bam), path(bam_index), path(snp_indel_vcf), path(snp_indel_vcf_index)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(data_type), path(bam), path(bam_index), path('clair3.snp_indel.phased.vcf.gz'), path('clair3.snp_indel.phased.vcf.gz.tbi')
        tuple val(sample_id), path('clair3.snp_indel.phased.vcf.gz'), path('clair3.snp_indel.phased.vcf.gz.tbi'), path('clair3.snp_indel.phased.read_list.txt')

    script:
        """
        # run whatshap phase
        whatshap phase \
        --reference $ref \
        --output clair3.snp_indel.phased.vcf.gz \
        --output-read-list clair3.snp_indel.phased.read_list.txt \
        --sample $sample_id \
        --ignore-read-groups $snp_indel_vcf $bam
        # index vcf
        tabix clair3.snp_indel.phased.vcf.gz
        """

    stub:
        """
        touch clair3.snp_indel.phased.vcf.gz
        touch clair3.snp_indel.phased.vcf.gz.tbi
        touch clair3.snp_indel.phased.read_list.txt
        """

}

process publish_whatshap_phase_clair3 {

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$filename" }

    input:
        tuple val(sample_id), path(snp_indel_phased_vcf), path(snp_indel_phased_vcf_index), path(snp_indel_phased_read_list)
        val outdir
        val outdir2
        val ref_name

    output:
        val sample_id
        path 'clair3.snp_indel.phased.vcf.gz'
        path 'clair3.snp_indel.phased.vcf.gz.tbi'
        path 'clair3.snp_indel.phased.read_list.txt'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch clair3.snp_indel.phased.vcf.gz
        touch clair3.snp_indel.phased.vcf.gz.tbi
        touch clair3.snp_indel.phased.read_list.txt
        """

}

process deepvariant {

    input:
        tuple val(sample_id), val(data_type), val(regions_of_interest), val(clair3_model), path(bam), path(bam_index)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(data_type), path(bam), path(bam_index), path('deepvariant.snp_indel.vcf.gz'), path('deepvariant.snp_indel.vcf.gz.tbi')
        tuple val(sample_id), path('deepvariant.snp_indel.g.vcf.gz'), path('deepvariant.snp_indel.g.vcf.gz.tbi'), path('deepvariant.version.txt')

    script:
    // define an optional string to pass regions of interest bed file
    def regions_of_interest_optional = file(regions_of_interest).name != 'NONE' ? "-L $regions_of_interest" : ''
        """
        # run deepvariant
        pbrun deepvariant \
        --mode $data_type \
        --ref $ref \
        --in-bam $bam \
        --gvcf \
        --out-variants deepvariant.snp_indel.g.vcf.gz \
        --phase-reads \
        --num-gpus ${task.gpus} \
        --num-cpu-threads-per-stream ${task.cpus} \
        --num-streams-per-gpu 1 \
        $regions_of_interest_optional
        # compress and index vcf
        bgzip \
        -@ ${task.cpus} \
        deepvariant.snp_indel.vcf
        tabix deepvariant.snp_indel.vcf.gz
        # grab version
        pbrun deepvariant --version >> deepvariant.version.txt
        """

    stub:
        """
        touch deepvariant.snp_indel.vcf.gz
        touch deepvariant.snp_indel.vcf.gz.tbi
        touch deepvariant.snp_indel.g.vcf.gz
        touch deepvariant.snp_indel.g.vcf.gz.tbi
        touch deepvariant.version.txt
        """

}

process publish_deepvariant {

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$filename" }

    input:
        tuple val(sample_id), path(snp_indel_gvcf), path(snp_indel_gvcf_index), path(deepvariant_version)
        val outdir
        val outdir2
        val ref_name

    output:
        val sample_id
        path 'deepvariant.snp_indel.g.vcf.gz'
        path 'deepvariant.snp_indel.g.vcf.gz.tbi'
        path 'deepvariant.version.txt'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch deepvariant.snp_indel.g.vcf.gz
        touch deepvariant.snp_indel.g.vcf.gz.tbi
        touch deepvariant.version.txt
        """

}

process whatshap_phase_dv {

    input:
        tuple val(sample_id), val(data_type), path(bam), path(bam_index), path(snp_indel_vcf), path(snp_indel_vcf_index)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(data_type), path(bam), path(bam_index), path('deepvariant.snp_indel.phased.vcf.gz'), path('deepvariant.snp_indel.phased.vcf.gz.tbi')
        tuple val(sample_id), path('deepvariant.snp_indel.phased.vcf.gz'), path('deepvariant.snp_indel.phased.vcf.gz.tbi'), path('deepvariant.snp_indel.phased.read_list.txt')

    script:
        """
        # run whatshap phase
        whatshap phase \
        --reference $ref \
        --output deepvariant.snp_indel.phased.vcf.gz \
        --output-read-list deepvariant.snp_indel.phased.read_list.txt \
        --sample $sample_id \
        --ignore-read-groups $snp_indel_vcf $bam
        # index vcf
        tabix deepvariant.snp_indel.phased.vcf.gz
        """

    stub:
        """
        touch deepvariant.snp_indel.phased.vcf.gz
        touch deepvariant.snp_indel.phased.vcf.gz.tbi
        touch deepvariant.snp_indel.phased.read_list.txt
        """

}

process publish_whatshap_phase_dv {

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$filename" }

    input:
        tuple val(sample_id), path(snp_indel_phased_vcf), path(snp_indel_phased_vcf_index), path(snp_indel_phased_read_list)
        val outdir
        val outdir2
        val ref_name

    output:
        val sample_id
        path 'deepvariant.snp_indel.phased.vcf.gz'
        path 'deepvariant.snp_indel.phased.vcf.gz.tbi'
        path 'deepvariant.snp_indel.phased.read_list.txt'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch deepvariant.snp_indel.phased.vcf.gz
        touch deepvariant.snp_indel.phased.vcf.gz.tbi
        touch deepvariant.snp_indel.phased.read_list.txt
        """

}

process whatshap_haplotag {

    input:
        tuple val(sample_id), val(data_type), path(bam), path(bam_index), path(snp_indel_vcf), path(snp_indel_vcf_index)
        val ref
        val ref_index

    output:
        tuple val(sample_id), val(data_type), path('minimap2.whatshap.sorted.haplotagged.bam'), path('minimap2.whatshap.sorted.haplotagged.bam.bai')
        tuple val(sample_id), path('minimap2.whatshap.sorted.haplotagged.bam'), path('minimap2.whatshap.sorted.haplotagged.bam.bai'), path('minimap2.whatshap.sorted.haplotagged.tsv'), path('whatshap.version.txt')

    script:
        """
        # run whatshap haplotag
        whatshap haplotag \
        --reference $ref \
        --output minimap2.whatshap.sorted.haplotagged.bam \
        --sample $sample_id \
        --tag-supplementary \
        --ignore-read-groups \
        --output-threads ${task.cpus} \
        --output-haplotag-list minimap2.whatshap.sorted.haplotagged.tsv \
        $snp_indel_vcf $bam
        # index bam
        samtools index \
        -@ ${task.cpus} \
        minimap2.whatshap.sorted.haplotagged.bam
        # grab version
        whatshap --version > whatshap.version.txt
        """

    stub:
        """
        touch minimap2.whatshap.sorted.haplotagged.bam
        touch minimap2.whatshap.sorted.haplotagged.bam.bai
        touch minimap2.whatshap.sorted.haplotagged.tsv
        touch whatshap.version.txt
        """

}

process publish_whatshap_haplotag {

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$filename" }

    input:
        tuple val(sample_id), path(haplotagged_bam), path(haplotagged_bam_index), path(haplotagged_bam_tsv), path(whatshap_version)
        val outdir
        val outdir2
        val ref_name

    output:
        path 'minimap2.whatshap.sorted.haplotagged.bam'
        path 'minimap2.whatshap.sorted.haplotagged.bam.bai'
        path 'minimap2.whatshap.sorted.haplotagged.tsv'
        path 'whatshap.version.txt'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch minimap2.whatshap.sorted.haplotagged.bam
        touch minimap2.whatshap.sorted.haplotagged.bam.bai
        touch minimap2.whatshap.sorted.haplotagged.tsv
        touch whatshap.version.txt
        """

}

process sniffles {

    input:
        tuple val(sample_id), val(data_type), path(haplotagged_bam), path(haplotagged_bam_index)
        val ref
        val ref_index
        val tandem_repeat

    output:
        tuple val(sample_id), path('sniffles.sv.phased.vcf.gz'), path('sniffles.sv.phased.vcf.gz.tbi'), path('sniffles.sv.phased.snf'), path('sniffles.version.txt')

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
        --vcf sniffles.sv.phased.vcf.gz \
        --snf sniffles.sv.phased.snf \
        --output-rnames \
        --minsvlen 20 \
        --phase $tandem_repeat_optional
        # grab version
        sniffles --version > sniffles.version.txt
        """

    stub:
        """
        touch sniffles.sv.phased.vcf.gz
        touch sniffles.sv.phased.vcf.gz.tbi
        touch sniffles.sv.phased.snf
        touch sniffles.version.txt
        """

}

process publish_sniffles {

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$filename" }

    input:
        tuple val(sample_id), path(sv_vcf), path(sv_vcf_index), path(sv_snf), path(sniffles_version)
        val outdir
        val outdir2
        val ref_name

    output:
        path 'sniffles.sv.phased.vcf.gz'
        path 'sniffles.sv.phased.vcf.gz.tbi'
        path 'sniffles.sv.phased.snf'
        path 'sniffles.version.txt'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch sniffles.sv.phased.vcf.gz
        touch sniffles.sv.phased.vcf.gz.tbi
        touch sniffles.sv.phased.snf
        touch sniffles.version.txt
        """

}

process cutesv {

    input:
        tuple val(sample_id), val(data_type), path(haplotagged_bam), path(haplotagged_bam_index)
        val ref
        val ref_index
        val tandem_repeat

    output:
        tuple val(sample_id), path('cutesv.sv.vcf.gz'), path('cutesv.sv.vcf.gz.tbi'), path('cutesv.version.txt')

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
        cutesv.sv.vcf \
        ./ \
        --sample ${sample_id} \
        -t ${task.cpus} \
        --genotype \
        --report_readid \
        $settings
        # compress and index vcf
        bgzip \
        -@ ${task.cpus} \
        cutesv.sv.vcf
        tabix cutesv.sv.vcf.gz
        # grab version
        cuteSV --version > cutesv.version.txt
        """

    stub:
        """
        touch cutesv.sv.vcf.gz
        touch cutesv.sv.vcf.gz.tbi
        touch cutesv.version.txt
        """

}

process publish_cutesv {

    publishDir "$outdir/$sample_id/$outdir2", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$filename" }

    input:
        tuple val(sample_id), path(sv_vcf), path(sv_vcf_index), path(cutesv_version)
        val outdir
        val outdir2
        val ref_name

    output:
        path 'cutesv.sv.vcf.gz'
        path 'cutesv.sv.vcf.gz.tbi'
        path 'cutesv.version.txt'

    script:
        """
        echo "Publishing files"
        """

    stub:
        """
        touch cutesv.sv.vcf.gz
        touch cutesv.sv.vcf.gz.tbi
        touch cutesv.version.txt
        """

}

workflow {

    // grab metrics from parameter file
    in_data = "$params.in_data"
    ref = "$params.ref"
    ref_index = "$params.ref_index"
    tandem_repeat = "$params.tandem_repeat"
    snp_indel_caller = "$params.snp_indel_caller"
    sv_caller = "$params.sv_caller"
    outdir = "$params.outdir"
    outdir2 = "$params.outdir2"

    // check user provided inputs
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
        exit 1, "Error: No SV calling software selected. Either include in parameter file or pass to --sv_caller on the command line. Should be either 'sniffles' or 'cutesv'."
    }
    if ( sv_caller != 'sniffles' && sv_caller != 'cutesv' ) {
        exit 1, "Error: SV calling software should be either 'sniffles' or 'cutesv', '${sv_caller}' selected."
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

    // build channel from input csv file
    // groupTuple will collapse by sample_id (as defined by the input csv file), creating a list of files per sample_id
    Channel
        .fromPath( in_data )
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple( row.sample_id, file(row.file).getExtension(), row.file, row.data_type, row.regions_of_interest, row.clair3_model ) }
        .groupTuple(by: [0,1,3,4,5] )
        .set { in_data_tuple }

    // also build a list of files NOT collaped by sample_id (as defined by the input csv file)
    Channel
        .fromPath( in_data )
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row-> tuple( row.sample_id, file(row.file).getExtension(), row.file ) }
        .set { in_data_list }

    // workflow
    scrape_settings_to_publish = scrape_settings(in_data_tuple, in_data, ref, ref_index, tandem_repeat, snp_indel_caller, sv_caller, outdir)
    publish_settings(scrape_settings_to_publish, outdir, outdir2, ref_name)
    bam_header = scrape_bam_header(in_data_list)
    publish_bam_header(bam_header, outdir, outdir2)
    merged = merge_runs(in_data_tuple)
    (bam, minimap_to_publish) = minimap2(merged, ref, ref_index)
    publish_minimap2(minimap_to_publish, outdir, outdir2, ref_name)
    if ( snp_indel_caller == 'clair3' ) {
        (snp_indel_vcf, clair3_to_publish) = clair3(bam, ref, ref_index)
        publish_clair3(clair3_to_publish, outdir, outdir2, ref_name)
        (snp_indel_phased_vcf, whatshap_phase_to_publish) = whatshap_phase_clair3(snp_indel_vcf, ref, ref_index)
        publish_whatshap_phase_clair3(whatshap_phase_to_publish, outdir, outdir2, ref_name)
    }
    else if ( snp_indel_caller == 'deepvariant' ) {
        (snp_indel_vcf, deepvariant_to_publish) = deepvariant(bam, ref, ref_index)
        publish_deepvariant(deepvariant_to_publish, outdir, outdir2, ref_name)
        (snp_indel_phased_vcf, whatshap_phase_to_publish) = whatshap_phase_dv(snp_indel_vcf, ref, ref_index)
        publish_whatshap_phase_dv(whatshap_phase_to_publish, outdir, outdir2, ref_name)
    }
    (haplotagged_bam, whatshap_haplotag_to_publish) = whatshap_haplotag(snp_indel_phased_vcf, ref, ref_index)
    publish_whatshap_haplotag(whatshap_haplotag_to_publish, outdir, outdir2, ref_name)
    if ( sv_caller == 'sniffles' ) {
        sniffles_to_publish = sniffles(haplotagged_bam, ref, ref_index, tandem_repeat)
        publish_sniffles(sniffles_to_publish, outdir, outdir2, ref_name)
    }
    else if ( sv_caller == 'cutesv' ) {
        cutesv_to_publish = cutesv(haplotagged_bam, ref, ref_index, tandem_repeat)
        publish_cutesv(cutesv_to_publish, outdir, outdir2, ref_name)
    }

}

