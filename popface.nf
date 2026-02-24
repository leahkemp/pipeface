nextflow.enable.dsl=2

// tag popface version
def popface_version = "dev"

// create dummy NONE file for optional pipeface inputs
new File("NONE").text = "Dummy file for optional pipeface inputs. Don't delete during a pipeline run unless you want a bad time.\n"

// set defaults for optional undocumented params
params.outdir2 = ""

process scrape_settings {

    publishDir "$outdir/$pop_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$pop_id.$filename" }, pattern: '*popface*'

    input:
        tuple val(pop_id), val(sample_ids), val(gvcfs), val(bams), val(sniffles), val(cutesv), val(somaliers), val(data_type)
        val popface_version
        val in_data
        val ref
        val ref_index
        val tr_calling
        val tr_call_regions
        val annotate
        val outdir
        val outdir2

    output:
        tuple val(pop_id), path('popface_settings.txt'), path('popface_input_files.tsv')

    script:
        """
        F1="popface_settings.txt"
        F2="popface_input_files.tsv"
        echo "Popface version: $popface_version" >> \${F1}
        echo "Population ID: $pop_id" >> \${F1}
        echo "In data csv path: $in_data" >> \${F1}
        echo "Reference genome: $ref" >> \${F1}
        echo "Reference genome index: $ref_index" >> \${F1}
        echo "Tandem repeat calling: $tr_calling" >> \${F1}
        echo "Tandem repeat call regions: $tr_call_regions" >> \${F1}
        echo "Data type: $data_type" >> \${F1}
        echo "Annotate: $annotate" >> \${F1}
        echo "Outdir: $outdir" >> \${F1}
        printf "sample_id\tgvcf\tbam\tsniffles\tcutesv\tsomalier\n" > \${F2}
        printf "%s\n" "$sample_ids" | sed "s/\\[//g; s/\\]//g; s/, /\\n/g" > sample_ids.txt
        printf "%s\n" "$gvcfs" | sed "s/\\[//g; s/\\]//g; s/, /\\n/g" > gvcfs.txt
        printf "%s\n" "$bams" | sed "s/\\[//g; s/\\]//g; s/, /\\n/g" > bams.txt
        printf "%s\n" "$sniffles" | sed "s/\\[//g; s/\\]//g; s/, /\\n/g" > sniffles.txt
        printf "%s\n" "$cutesv" | sed "s/\\[//g; s/\\]//g; s/, /\\n/g" > cutesv.txt
        printf "%s\n" "$somaliers" | sed "s/\\[//g; s/\\]//g; s/, /\\n/g" > somaliers.txt
        paste sample_ids.txt gvcfs.txt bams.txt sniffles.txt cutesv.txt somaliers.txt >> \${F2}
        """

    stub:
        """
        touch popface_settings.txt
        touch popface_input_files.tsv
        """

}

process glnexus_pre_processing {

    input:
        tuple val(pop_id), val(sample_id), path(gvcf)
        val (ref_index)

    output:
        tuple val(pop_id), val(sample_id), path("*.ammended.g.vcf.gz")

    script:
        """
        # reheader to include all contigs in gvcf header even if no variants were called in a contig to avoid this issue: https://github.com/HKU-BAL/Clair3/issues/371
        {
            bcftools view -h "$gvcf" | grep -v -E '^##contig=|^#CHROM'
            awk '{ printf "##contig=<ID=%s,length=%s>\\n", \$1, \$2 }' "$ref_index"
            bcftools view -h "$gvcf" | grep '#CHROM'
        } > ${sample_id}.ammended.g.vcf
        # convert lower cases of soft-masked sequences to upper case to avoid this issue: https://github.com/HKU-BAL/Clair3/issues/359
        zgrep -v '#' $gvcf | awk -F'\\t' -v OFS='\\t' '{ if(\$0 !~ /^#/) { \$4=toupper(\$4); \$5=toupper(\$5) } print }' >> ${sample_id}.ammended.g.vcf
        # compress vcf
        bgzip -@ ${task.cpus} ${sample_id}.ammended.g.vcf
        """

    stub:
        """
        touch sample1.ammended.g.vcf.gz
        """

}

process glnexus {

    input:
        tuple val(pop_id), val(sample_ids), path(gvcfs)
        val snp_indel_caller
        val clair3_config

    output:
        tuple val(pop_id), path('snp_indel.bcf')

    script:
        if (snp_indel_caller == 'clair3')
        """
        glnexus_cli --config $clair3_config $gvcfs > snp_indel.bcf
        """
        else if (snp_indel_caller == 'deepvariant')
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
        tuple val(pop_id), path(joint_snp_indel_bcf)

    output:
        tuple val(pop_id), path('snp_indel.raw.vcf.gz'), path('snp_indel.raw.vcf.gz.tbi')

    script:
        """
        bcftools view $joint_snp_indel_bcf | bgzip -@ ${task.cpus} -c > snp_indel.raw.vcf.gz
        tabix snp_indel.raw.vcf.gz
        """

    stub:
        """
        touch snp_indel.raw.vcf.gz
        touch snp_indel.raw.vcf.gz.tbi
        """

}

process split_multiallele {

    input:
        tuple val(pop_id), path(snp_indel_vcf), path(snp_indel_vcf_index)
        val ref
        val ref_index

    output:
        tuple val(pop_id), path('snp_indel.split.vcf.gz'), path('snp_indel.split.vcf.gz.tbi')

    script:
        """
        # run bcftools norm
        bcftools norm --threads ${task.cpus} -m -any -f $ref $snp_indel_vcf > snp_indel.split.vcf
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

process split_vcf {

    input:
        tuple val(pop_id), val(sample_id), path(joint_snp_indel_vcf), path(joint_snp_indel_vcf_index)

    output:
        tuple val(pop_id), val(sample_id), path('*snp_indel.vcf.gz'), path('*snp_indel.vcf.gz.tbi')

    script:
        """
        bcftools view -s $sample_id $joint_snp_indel_vcf -Oz -o "${sample_id}".snp_indel.vcf.gz
        tabix "${sample_id}".snp_indel.vcf.gz
        """

    stub:
        """
        touch sample1.snp_indel.vcf.gz
        touch sample1.snp_indel.vcf.gz.tbi
        """

}

process whatshap_phase {

    publishDir "$outdir/$pop_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$snp_indel_caller.$filename"}, pattern: '*.read_list.txt'
    publishDir "$outdir/$pop_id/$outdir2/$sample_id", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$snp_indel_caller.$filename"}, pattern: '*.stats.gtf'

    input:
        tuple val(pop_id), val(sample_id), path(snp_indel_vcf), path(snp_indel_vcf_index), path(bam), path(bam_index)
        val ref
        val ref_index
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller

    output:
        tuple val(pop_id), path("${sample_id}.snp_indel.phased.vcf.gz"), path("${sample_id}.snp_indel.phased.vcf.gz.tbi")
        tuple val(pop_id), path('snp_indel.phased.read_list.txt'), path('snp_indel.phased.stats.gtf')

    script:
        """
        # run whatshap phase
        whatshap phase --reference $ref --output snp_indel.phased.vcf.gz --output-read-list snp_indel.phased.read_list.txt --sample $sample_id --ignore-read-groups $snp_indel_vcf $bam
        # index vcf
        tabix snp_indel.phased.vcf.gz
        # run whatshap stats
        whatshap stats snp_indel.phased.vcf.gz --gtf snp_indel.phased.stats.gtf --sample $sample_id
        # tag vcf with sample_id for downstream vcf merge
        ln -s snp_indel.phased.vcf.gz ${sample_id}.snp_indel.phased.vcf.gz
        ln -s snp_indel.phased.vcf.gz.tbi ${sample_id}.snp_indel.phased.vcf.gz.tbi
        """

    stub:
        """
        touch sample1.snp_indel.phased.vcf.gz
        touch sample1.snp_indel.phased.vcf.gz.tbi
        """

}

process merge_vcf {

    publishDir "$outdir/$pop_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$pop_id.$ref_name.$snp_indel_caller.$filename" }, pattern: 'snp_indel.phased.*'

    input:
        tuple val(pop_id), path(snp_indel_phased_vcfs), path(snp_indel_phased_vcf_indicies)
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller

    output:
        tuple val(pop_id), path('snp_indel.phased.vcf.gz'), path('snp_indel.phased.vcf.gz.tbi')

    script:
        """
        # merge vcf
        bcftools merge -Oz -o snp_indel.phased.vcf.gz ./*.snp_indel.phased.vcf.gz --threads ${task.cpus}
        # index vcf
        tabix snp_indel.phased.vcf.gz
        """

    stub:
        """
        touch snp_indel.phased.vcf.gz
        touch snp_indel.phased.vcf.gz.tbi
        """


}

process split_sv_vcf {

    input:
        tuple val(pop_id), val(sample_id), val(sv_caller), path(sv_vcf), path(sv_vcf_index)
        path ref_index

    output:
        tuple val(pop_id), val(sv_caller), path("${sample_id}.${sv_caller}.*.vcf")

    script:
        """
        # extract BND variants (they can span multiple chromosomes)
        bcftools view -i 'TYPE="BND"' $sv_vcf -o ${sample_id}.${sv_caller}.BND.vcf
        # split rest of variants per chromosome
        cut -f1 $ref_index | while read CHR; do
            bcftools view -r \${CHR} -e 'TYPE="BND"' $sv_vcf -o ${sample_id}.${sv_caller}.\${CHR}.vcf
        done
        """

    stub:
        """
        echo chr1 | while read CHR; do
            touch ${sample_id}.${sv_caller}.\${CHR}.vcf
        done
        """

}

process jasmine {

    input:
        tuple val(pop_id), val(chr), val(sv_caller), path(split_sv_vcfs), path(bams), val(data_type)
        val ref
        val ref_index

    output:
        tuple val(pop_id), val(sv_caller), path("${chr}.sv.vcf.gz"), path("${chr}.sv.vcf.gz.tbi"), optional: true

    script:
        // conditionally define iris arguments
        // as default, iris will pass minimap -x map-ont
        // the --pacbio flag passed to iris will pass minimap -x map-pb
        // these are the only two options iris makes available for minimaps -x argument, so I can't use lr:hq and map-hifi boo
        def iris_args = ''
        if (data_type == 'ont') {
            iris_args = '--run_iris iris_args=min_ins_length=20,--rerunracon,--keep_long_variants'
        }
        else if (data_type == 'pacbio') {
            iris_args = '--run_iris iris_args=min_ins_length=20,--rerunracon,--keep_long_variants,--pacbio'
        }
        """
        # create file lists
        for VCF in ${split_sv_vcfs}; do realpath \${VCF} >> vcfs.txt; done
        for BAM in ${bams}; do realpath \${BAM} >> bams.txt; done
        # run jasmine
        jasmine threads=${task.cpus} out_dir=./ genome_file=$ref file_list=vcfs.txt bam_list=bams.txt out_file=sv.tmp.vcf min_support=1 --mark_specific spec_reads=7 spec_len=20 --pre_normalize --output_genotypes --clique_merging --dup_to_ins --normalize_type --require_first_sample --default_zero_genotype $iris_args
        # post process vcf if non-empty
        if [ \$(bcftools view -H sv.tmp.vcf | wc -l) -gt 0 ]; then
            # fix vcf header (remove prefix to sample names that jasmine adds)
            grep '##' sv.tmp.vcf > sv.vcf
            grep '#CHROM' sv.tmp.vcf | sed -E 's/\\b[0-9]+_//g' >> sv.vcf
            grep -v '#' sv.tmp.vcf >> sv.vcf
            bcftools sort sv.vcf -o sorted_sv.vcf
            # compress and index vcf
            bgzip -c sorted_sv.vcf > ${chr}.sv.vcf.gz
            tabix ${chr}.sv.vcf.gz
        fi
        """

    stub:
        """
        touch sv.vcf.gz
        touch sv.vcf.gz.tbi
        """

}

process concat_sv_vcf {

    merger = "jasmine"

    publishDir "$outdir/$pop_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$pop_id.$ref_name.$sv_caller.$merger.$filename"}, pattern: 'sv.*vcf.gz*'

    input:
        tuple val(pop_id), val(sv_caller), path(sv_vcfs), path(sv_vcf_indicies)
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(pop_id), val(sv_caller), path('sv.*vcf.gz')
        tuple val(pop_id), path('sv.*vcf.gz'), path('sv.*vcf.gz.tbi')

    script:
        """
        # get list of vcfs
        VCFS=(*sv*.vcf.gz)
        # concat vcfs, or symlink if only one vcf
        if [[ \${#VCFS[@]} -eq 1 ]]; then
            ln -s \${VCFS[0]} sv.vcf.gz
        else
            bcftools concat -a -Oz -o sv.vcf.gz \${VCFS[@]} --threads ${task.cpus}
        fi
        # index vcf
        tabix sv.vcf.gz
        # rename file for publishing purposes
        if [ $sv_caller == "sniffles" ]; then
            mv sv.vcf.gz sv.phased.vcf.gz
            mv sv.vcf.gz.tbi sv.phased.vcf.gz.tbi
        fi
        """

    stub:
        """
        touch sv.vcf.gz
        touch sv.vcf.gz.tbi
        """

}

process somalier {

    publishDir "$outdir/$pop_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$pop_id.$ref_name.$filename" }, pattern: 'somalier*'

    input:
        tuple val(pop_id), path(somaliers)
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(pop_id), path('somalier.samples.tsv'), path('somalier.pairs.tsv'), path('somalier.html')

    script:
        """
        # run somalier relate
        somalier relate $somaliers
        """

}

process longtr_pre_processing {

    input:
        val tr_call_regions

    output:
        path('split.*.bed')

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

    def software = "longtr"

    input:
        tuple val(pop_id), val(sample_ids), path(bams), path(bam_indices), val(data_type), path(split_bed)
        val ref
        val ref_index

    output:
        tuple val(pop_id), path('tr.*.vcf.gz'), path('tr.*.vcf.gz.tbi')

    script:
        // define a string to define non-default alignment parameters for ONT to account for higher incidence of indels in homopolymers (defaults are tailored to pacbio hifi)
        def alignment_params_optional = data_type == 'ont' ? "--alignment-params -1.0,-0.458675,-1.0,-0.458675,-0.00005800168,-1,-1" : ''
        """
        # comma seperate list of bams
        BAMS=\$(echo $bams | tr ' ' ',')
        SAMPLE_IDS=\$(echo $sample_ids | tr -d '[ ]')
        # run longtr
        ID=\$(echo $split_bed | sed 's/split.//;s/.bed//')
        LongTR --bams \${BAMS} --bam-samps \${SAMPLE_IDS} --bam-libs \${SAMPLE_IDS} --fasta $ref --regions $split_bed --tr-vcf tr.\${ID}.vcf.gz --phased-bam --output-gls --output-pls --output-phased-gls --output-filter $alignment_params_optional --log longtr.log
        # index vcf
        tabix tr.\${ID}.vcf.gz
        """

    stub:
        """
        touch tr.aa.vcf.gz
        touch tr.aa.vcf.gz.tbi
        """

}

process concat_tr_vcf {

    def software = "longtr"

    publishDir "$outdir/$pop_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$pop_id.$ref_name.$software.$filename" }, pattern: 'tr.vcf.gz*'

    input:
        tuple val(pop_id), path(tr_vcfs), path(tr_vcf_indicies)
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(pop_id), path('tr.vcf.gz'), path('tr.vcf.gz.tbi')

    script:
        """
        # get list of vcfs
        VCFS=(tr.*.vcf.gz)
        # concat vcfs, or symlink if only one vcf
        if [[ \${#VCFS[@]} -eq 1 ]]; then
            ln -s \${VCFS[0]} tr.vcf.gz
        else
            bcftools concat -a -Oz -o tr.vcf.gz \${VCFS[@]} --threads ${task.cpus}
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

process vep_snp_indel {

    publishDir "$outdir/$pop_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$pop_id.$ref_name.$snp_indel_caller.$filename"}, pattern: 'snp_indel.phased.annotated.vcf.gz*'

    input:
        tuple val(pop_id), path(joint_snp_indel_phased_vcf), path(joint_snp_indel_phased_vcf_index)
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
        tuple val(pop_id), path('snp_indel.phased.annotated.vcf.gz'), path('snp_indel.phased.annotated.vcf.gz.tbi')

    script:
        """
        # run vep
        vep -i $joint_snp_indel_phased_vcf -o snp_indel.phased.annotated.vcf.gz --format vcf --vcf --fasta $ref --dir $vep_db --assembly GRCh38 --species homo_sapiens --cache --offline --merged --sift b --polyphen b --symbol --hgvs --hgvsg --uploaded_allele --check_existing --filter_common --distance 0 --nearest gene --canonical --mane --pick --fork ${task.cpus} --no_stats --compress_output bgzip --dont_skip \
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

process vep_sv {

    merger = "jasmine"

    publishDir "$outdir/$pop_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$pop_id.$ref_name.$sv_caller.$merger.$filename"}, pattern: 'sv.*annotated.vcf.gz*'

    input:
        tuple val(pop_id), val(sv_caller), path(sv_vcf)
        val ref
        val ref_index
        val vep_db
        val gnomad_db
        val cadd_sv_db
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(pop_id), path('sv.*annotated.vcf.gz'), path('sv.*annotated.vcf.gz.tbi')

    script:
        """
        # run vep
        vep -i $sv_vcf -o sv.annotated.vcf.gz --format vcf --vcf --fasta $ref --dir $vep_db --assembly GRCh38 --species homo_sapiens --cache --offline --merged --sift b --polyphen b --symbol --hgvs --hgvsg --uploaded_allele --check_existing --filter_common --distance 0 --nearest gene --canonical --mane --pick --fork ${task.cpus} --no_stats --compress_output bgzip --dont_skip --plugin CADD,sv=$cadd_sv_db
        # index vcf
        tabix sv.annotated.vcf.gz
        # rename file for publishing purposes
        if [ $sv_caller == "sniffles" ]; then
            mv sv.annotated.vcf.gz sv.phased.annotated.vcf.gz
            mv sv.annotated.vcf.gz.tbi sv.phased.annotated.vcf.gz.tbi
        fi
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
    ref = "${params.ref}".trim()
    ref_index = "${params.ref_index}".trim()
    tr_calling = "${params.tr_calling}".trim()
    tr_call_regions = "${params.tr_call_regions}".trim()
    snp_indel_caller = "${params.snp_indel_caller}".trim()
    clair3_config = "${params.clair3_config}".trim()
    annotate = "${params.annotate}".trim()
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
    if (!in_data) {
        exit 1, "No in data csv file (in_data) provided."
    }
    if (!ref) {
        exit 1, "No reference genome file (ref) provided."
    }
    if (!ref_index) {
        exit 1, "No reference genome index file (ref_index) provided."
    }
    if (!(annotate in ['yes', 'no'])) {
        exit 1, "Choice to annotate should be either 'yes' or 'no', annotate = '${annotate}' provided."
    }
    if (!outdir) {
        exit 1, "No output directory (outdir) provided."
    }
    if (!file(in_data).exists()) {
        exit 1, "In data csv file does not exist, 'in_data = ${in_data}' provided."
    }
    if (!file(ref).exists()) {
        exit 1, "Reference genome file does not exist, ref = '${ref}' provided."
    }
    if (!file(ref_index).exists()) {
        exit 1, "Reference genome index file does not exist, ref_index = '${ref_index}' provided."
    }
    if (!snp_indel_caller) {
        exit 1, "No SNP/indel caller (snp_indel_caller) provided."
    }
    if (!(snp_indel_caller in ['clair3', 'deepvariant'])) {
        exit 1, "SNP/indel caller should be either 'clair3' or 'deepvariant', snp_indel_caller = '${snp_indel_caller}' provided."
    }
    if (!clair3_config) {
        exit 1, "No clair3 config file (clair3_config) provided."
    }
    if (!file(clair3_config).exists()) {
        exit 1, "Clair3 config file does not exist, 'clair3_config = ${clair3_config}' provided."
    }
    if (snp_indel_caller == 'clair3' && clair3_config == 'NONE') {
        exit 1, "When clair3 is selected as the SNP/indel calling software, provide a path to an appropriate clair3 config file for GLnexus, snp_indel_caller = '${snp_indel_caller}' and clair3_config = '${clair3_config}' provided."
    }
    if (snp_indel_caller == 'deepvariant' && clair3_config != 'NONE') {
        exit 1, "When deepvariant is selected as the SNP/indel calling software, set the clair3 config file for GLnexus to 'NONE', snp_indel_caller = '${snp_indel_caller}' and clair3_config = '${clair3_config}' provided."
    }
    if (!(tr_calling in ['yes', 'no'])) {
        exit 1, "Choice to call tandem repeats should be either 'yes', or 'no', tr_calling = '${tr_calling}' provided."
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
    if (!file(tr_call_regions).exists()) {
        exit 1, "tandem repeat call regions file does not exist, tr_call_regions = '${tr_call_regions}' provided. Set to 'NONE' if not required."
    }
    if (tr_calling == 'yes') {
        if (tr_call_regions == "NONE") {
            exit 1, "When calling tandem repeats, provide a valid tandem repeat call regions file, tr_calling = '${tr_calling}' and tr_call_regions = 'NONE' provided."
        }
    }
    else if (tr_calling == 'no') {
        if (tr_call_regions != 'NONE') {
            exit 1, "When not calling tandem repeats, set tandem repeat call regions file to 'NONE', tr_calling = '${tr_calling}' and tr_call_regions = '${tr_call_regions}' provided."
        }
    }
    if (!(params.publish_mode in ['copy', 'copyNoFollow', 'link', 'move', 'rellink', 'symlink'])) {
        exit 1, "Choice of publishing mode should be 'copy', 'copyNoFollow', 'link', 'move', 'rellink' or 'symlink', publish_mode = '$params.publish_mode' provided."
    }

    // build variable
    ref_name = file(ref).getSimpleName()

    // build input tuple
    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            def pop_id = row.pop_id
            def sample_ids = row.sample_id
            def gvcfs = row.gvcf
            def bams = row.bam
            def sniffless = row.sniffles
            def cutesvs = row.cutesv
            def somaliers = row.somalier
            def data_type = row.data_type
            return tuple(pop_id, sample_ids, gvcfs, bams, sniffless, cutesvs, somaliers, data_type)
        }
        .groupTuple(by: [0,7])
        .set { in_data_tuple }

    // build channels from in_data.csv file for describing each metadata associated with populations
    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .filter { row -> row.gvcf != 'NONE' }
        .map { row -> tuple(row.pop_id, row.sample_id) }
        .set { id_tuple }

    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .filter { row -> row.gvcf != 'NONE' }
        .map { row -> tuple(row.pop_id, row.sample_id, row.gvcf) }
        .set { gvcfs_tuple }

    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .filter { row -> row.gvcf != 'NONE' }
        .map { row -> tuple(row.pop_id, row.sample_id, row.bam, "${row.bam}.bai") }
        .set { bams_tuple }

    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .flatMap { row ->
            def tuples = []
            if (row.sniffles != 'NONE') tuples << tuple(row.pop_id, row.sample_id, 'sniffles', row.sniffles, "${row.sniffles}.tbi")
            if (row.cutesv != 'NONE') tuples << tuple(row.pop_id, row.sample_id, 'cutesv', row.cutesv, "${row.cutesv}.tbi")
            return tuples
        }
        .set { svs_tuple }

    Channel
    .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .filter { row -> row.somalier != 'NONE' }
        .map { row -> tuple(row.pop_id, row.somalier) }
        .groupTuple(by: 0)
        .set { somaliers_tuple }

    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .filter { row -> row.bam != 'NONE' }
        .map { row -> tuple(row.pop_id, row.sample_id, file(row.bam), file("${row.bam}.bai"), row.data_type) }
        .groupTuple(by: [0,4])
        .set { bams_data_type_tuple }

    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .filter { row -> row.bam != 'NONE' }
        .map { row -> tuple(row.pop_id, file(row.bam), row.data_type) }
        .groupTuple(by: [0,2])
        .set { bams_data_type_tuple2 }

    // check user provided parameters in in_data.csv file
    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            def pop_id = row.pop_id
            def sample_id = row.sample_id
            def gvcf = row.gvcf
            def bam = row.bam
            def sniffles = row.sniffles
            def cutesv = row.cutesv
            def somalier = row.somalier
            def data_type = row.data_type

            if (pop_id.isEmpty()) {
                exit 1, "There is an empty entry in the 'pop_id' column of '${in_data}'."
            }
            if (sample_id.isEmpty()) {
                exit 1, "There is an empty entry in the 'sample_id' column of '${in_data}'."
            }
            if (gvcf.isEmpty()) {
                exit 1, "There is an empty entry in the 'gvcf' column of '${in_data}'. Set to 'NONE' if not required."
            }
            if (!file(gvcf).exists()) {
                exit 1, "There is an entry in the 'gvcf' column of '${in_data}' which doesn't exist. Check file '${gvcf}'."
            }
            if (bam.isEmpty()) {
                exit 1, "There is an empty entry in the 'bam' column of '${in_data}'. Set to 'NONE' if not required."
            }
            if (!file(bam).exists()) {
                exit 1, "There is an entry in the 'bam' column of '${in_data}' which doesn't exist. Check file '${bam}'."
            }
            if (gvcf != 'NONE' && bam == 'NONE') {
                exit 1, "When a GVCF file is provided in the 'gvcf' column of '${in_data}', the associated BAM file must be provided in the 'bam' column. gvcf = '${gvcf}' and bam = '${bam}' provided."
            }
            if (bam != 'NONE') {
                def bai_exists = file("${bam}.bai").exists()
                def crai_exists = file("${bam}.crai").exists()
                if (!bai_exists && !crai_exists) {
                    exit 1, "There is an entry in the 'bam' column of '${in_data}' which doesn't look to have an associated index. Expecting '${bam}.bai' or '${bam}.crai'."
                }
            }
            if (sniffles.isEmpty()) {
                exit 1, "There is an empty entry in the 'sniffles' column of '${in_data}'. Set to 'NONE' if not required."
            }
            if (!file(sniffles).exists()) {
                exit 1, "There is an entry in the 'sniffles' column of '${in_data}' which doesn't exist. Check file '${sniffles}'."
            }
            if (sniffles != 'NONE' && bam == 'NONE') {
                exit 1, "When a Sniffles SV VCF is provided in the 'sniffles' column of '${in_data}', the associated BAM file must be provided in the 'bam' column. sniffles = '${sniffles}' and bam = '${bam}' provided."
            }
            if (sniffles != 'NONE') {
                if (!file("${sniffles}.tbi").exists()) {
                    exit 1, "There is an entry in the 'sniffles' column of '${in_data}' which doesn't look to have an associated index. Expecting '${sniffles}.tbi'."
                }
            }
            if (cutesv.isEmpty()) {
                exit 1, "There is an empty entry in the 'cutesv' column of '${in_data}'. Set to 'NONE' if not required."
            }
            if (!file(cutesv).exists()) {
                exit 1, "There is an entry in the 'cutesv' column of '${in_data}' which doesn't exist. Check file '${cutesv}'."
            }
            if (cutesv != 'NONE' && bam == 'NONE') {
                exit 1, "When a cuteSV VCF is provided in the 'cutesv' column of '${in_data}', the associated BAM file must be provided in the 'bam' column. cutesv = '${cutesv}' and bam = '${bam}' provided."
            }
            if (cutesv != 'NONE') {
                if (!file("${cutesv}.tbi").exists()) {
                    exit 1, "There is an entry in the 'cutesv' column of '${in_data}' which doesn't look to have an associated index. Expecting '${cutesv}.tbi'."
                }
            }
            if (somalier.isEmpty()) {
                exit 1, "There is an empty entry in the 'somalier' column of '${in_data}'. Set to 'NONE' if not required."
            }
            if (!file(somalier).exists()) {
                exit 1, "There is an entry in the 'somalier' column of '${in_data}' which doesn't exist. Check file '${somalier}'."
            }
            if (data_type.isEmpty()) {
                exit 1, "There is an empty entry in the 'data_type' column of '${in_data}'."
            }
            if (!(data_type in ['ont', 'pacbio'])) {
                exit 1, "Entries in the 'data_type' column of '${in_data}' should be 'ont' or 'pacbio', '${data_type}' provided."
            }
        }

    // check user provided parameters relating in in_data.csv file relating to populations
    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row -> tuple(row.pop_id, row.sample_id, row.gvcf, row.bam, row.sniffles, row.cutesv, row.somalier, row.data_type) }
        .groupTuple(by: 0)
        .map { pop_id, sample_ids, gvcfs, bams, sniffless, cutesvs, somaliers, data_types ->
            def gvcfs_any_none = gvcfs.any { it == 'NONE' }
            def gvcfs_all_none = gvcfs.every { it == 'NONE' }
            def bams_any_none = bams.any { it == 'NONE' }
            def bams_all_none = bams.every { it == 'NONE' }
            def sniffless_any_none = sniffless.any { it == 'NONE' }
            def sniffless_all_none = sniffless.every { it == 'NONE' }
            def cutesvs_any_none = cutesvs.any { it == 'NONE' }
            def cutesvs_all_none = cutesvs.every { it == 'NONE' }
            def somaliers_any_none = somaliers.any { it == 'NONE' }
            def somaliers_all_none = somaliers.every { it == 'NONE' }
            if (sample_ids.unique().size() < 2) {
                exit 1, "Entries in the 'sample_id' column of '$in_data' should contain 2 or more unique values for every 'pop_id', '$sample_ids' provided for population '$pop_id'."
            }
            if (gvcfs_any_none && !gvcfs_all_none) {
                exit 1, "Entries in the 'gvcf' column of '$in_data' must be either all 'NONE' or all real files for a given 'pop_id', '$gvcfs' provided for population '$pop_id'."
            }
            if (bams_any_none && !bams_all_none) {
                exit 1, "Entries in the 'bam' column of '$in_data' must be either all 'NONE' or all real files for a given 'pop_id', '$bams' provided for population '$pop_id'."
            }
            if (sniffless_any_none && !sniffless_all_none) {
                exit 1, "Entries in the 'sniffles' column of '$in_data' must be either all 'NONE' or all real files for a given 'pop_id', '$sniffless' provided for population '$pop_id'."
            }
            if (cutesvs_any_none && !cutesvs_all_none) {
                exit 1, "Entries in the 'cutesv' column of '$in_data' must be either all 'NONE' or all real files for a given 'pop_id', '$cutesvs' provided for population '$pop_id'."
            }
            if (somaliers_any_none && !somaliers_all_none) {
                exit 1, "Entries in the 'somalier' column of '$in_data' must be either all 'NONE' or all real files for a given 'pop_id', '$somaliers' provided for population '$pop_id'."
            }
            if (data_types.unique().size() != 1) {
                exit 1, "Entries in the 'data_type' column of '$in_data' must be identical for a given 'pop_id'. '$data_types' provided for population '$pop_id'."
            }
        }

    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row -> row.sample_id }
        .collect()
        .map { sample_ids ->
            def duplicates = sample_ids.findAll { id -> sample_ids.count(id) > 1 }.toSet()
            if (duplicates) {
                exit 1, "Entries in the 'sample_id' column of '$in_data' should all be unique values, duplicates: '$duplicates'."
            }
        }

    if (tr_calling == 'yes') {
        Channel
            .fromPath(in_data)
            .splitCsv(header: true, sep: ',', strip: true)
            .map { row -> row.bam }
            .collect()
            .map { bams ->
                def any_real_bam = bams.any { it != 'NONE' }
                if (!any_real_bam) {
                    exit 1, "When tandem repeat calling is turned on, relevant BAM files need to be provided in 'bam' column of '${in_data}', tr_calling = '${tr_calling}' and '${bams}' provided."
                }
            }
    }

        // workflow
        // pre-process
        scrape_settings(in_data_tuple, popface_version, in_data, ref, ref_index, tr_calling, tr_call_regions, annotate, outdir, outdir2)
        // gvcf merging
        if (snp_indel_caller == 'clair3') {
            gvcfs = glnexus_pre_processing(gvcfs_tuple, ref_index).groupTuple(by: 0)
        }
        if (snp_indel_caller == 'deepvariant') {
            gvcfs = gvcfs_tuple.groupTuple(by: 0)
        }
        joint_snp_indel_bcf = glnexus(gvcfs, snp_indel_caller, clair3_config)
        joint_snp_indel_vcf = glnexus_post_processing(joint_snp_indel_bcf)
        // split multiallelic variants
        joint_snp_indel_split_vcf = split_multiallele(joint_snp_indel_vcf, ref, ref_index)
        // joint phasing
        joint_snp_indel_vcf_id = joint_snp_indel_split_vcf
            .combine(id_tuple)
            .map { pop_id, joint_vcf, joint_vcf_index, pop_id2, sample_id ->
                if (pop_id != pop_id2) {
                    return null
                }
                tuple(pop_id, sample_id, joint_vcf, joint_vcf_index)
            }
        snp_indel_vcf = split_vcf(joint_snp_indel_vcf_id)
        (snp_indel_phased_vcfs, stats) = whatshap_phase(snp_indel_vcf.join(bams_tuple, by: [0,1]), ref, ref_index, outdir, outdir2, ref_name, snp_indel_caller)
        vcfs = snp_indel_phased_vcfs.groupTuple(by: 0)
        joint_snp_indel_phased_vcf = merge_vcf(vcfs, outdir, outdir2, ref_name, snp_indel_caller)
        // sv vcf merging
        split_sv_vcfs = split_sv_vcf(svs_tuple, ref_index)
        jasmine_input = split_sv_vcfs
            .flatMap { pop_id, sv_caller, vcfs ->
                vcfs.collectMany { vcf ->
                    def chr = vcf.baseName.tokenize('.')[-1]
                    [ tuple(pop_id, chr, sv_caller, vcf) ]
                }
            }
            .groupTuple(by: [0,1,2])
            .combine(bams_data_type_tuple2, by: 0)
        merged_sv_vcfs = jasmine(jasmine_input, ref, ref_index)
        (joint_sv_vcf, joint_sv_vcf_indexed) = concat_sv_vcf(merged_sv_vcfs.groupTuple(by: [0,1]), outdir, outdir2, ref_name)
        // joint somalier
        somalier(somaliers_tuple, outdir, outdir2, ref_name)
        // joint tr calling
        if (tr_calling == 'yes') {
            split_bed = longtr_pre_processing(tr_call_regions).flatten()
            longtr_input = bams_data_type_tuple
                .combine(split_bed)
                .map { pop_id, sample_ids, bams, bam_indices, data_type, split_bed ->
                    tuple(pop_id, sample_ids, bams, bam_indices, data_type, split_bed)
                }
            tr_vcfs = longtr(longtr_input, ref, ref_index)
            concat_tr_vcf(tr_vcfs.groupTuple(by: 0), outdir, outdir2, ref_name)
        }
        // joint annotation
        if (annotate == 'yes') {
            vep_snp_indel(joint_snp_indel_phased_vcf, ref, ref_index, vep_db, revel_db, gnomad_db, clinvar_db, cadd_snv_db, cadd_indel_db, spliceai_snv_db, spliceai_indel_db, alphamissense_db, outdir, outdir2, ref_name, snp_indel_caller)
            vep_sv(joint_sv_vcf, ref, ref_index, vep_db, gnomad_db, cadd_sv_db, outdir, outdir2, ref_name)
        }

}
