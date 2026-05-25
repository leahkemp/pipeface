nextflow.enable.dsl=2

// tag popface version
def popface_version = "0.11.0"

// set defaults for undocumented params
params.outdir2 = ""
params.annotate_override = ""
params.min_gap = 100000
params.chunks = 5

process scrape_settings {

    publishDir "$outdir/$pop_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$pop_id.$filename" }, pattern: '*popface*'

    input:
        tuple val(pop_id), val(sample_ids), val(related), val(family_positions), val(gvcfs), val(bams), val(sniffles), val(cutesv), val(somaliers), val(data_type)
        val popface_version
        val in_data
        val ref
        val ref_index
        val snp_indel_caller
        val tr_calling
        val tr_call_regions
        val annotate
        val outdir
        val outdir2

    output:
        tuple val(pop_id), path("popface_settings.txt"), path("popface_input_files.csv")

    script:
        """
        {
            echo "Popface version: $popface_version"
            echo "Population ID: $pop_id"
            echo "In data csv path: $in_data"
            echo "Reference genome: $ref"
            echo "Reference genome index: $ref_index"
            echo "SNP/indel caller: $snp_indel_caller"
            echo "Tandem repeat calling: $tr_calling"
            echo "Tandem repeat call regions: $tr_call_regions"
            echo "Data type: $data_type"
            echo "Related: $related"
            echo "Annotate: $annotate"
            echo "Outdir: $outdir"
        } > popface_settings.txt
        printf "sample_id,family_position,gvcf,bam,sniffles,cutesv,somalier\\n" > popface_input_files.csv
        NAMES=(sample_ids family_positions gvcfs bams sniffles cutesv somaliers)
        VALUES=("$sample_ids" "$family_positions" "$gvcfs" "$bams" "$sniffles" "$cutesv" "$somaliers")
        for i in \${!NAMES[@]}; do
            printf "%s\n" "\${VALUES[\$i]}" | sed "s/\\[//g; s/\\]//g; s/, /\\n/g" > "\${NAMES[\$i]}.txt"
        done
        paste -d',' sample_ids.txt family_positions.txt gvcfs.txt bams.txt sniffles.txt cutesv.txt somaliers.txt >> popface_input_files.csv
        """

    stub:
        """
        touch popface_settings.txt
        touch popface_input_files.csv
        """

}

process glnexus_pre_processing {

    input:
        tuple val(pop_id), val(sample_id), path(gvcf)
        val (ref_index)

    output:
        tuple val(pop_id), val(sample_id), path("*.amended.g.vcf.gz")

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
        touch sample1.amended.g.vcf.gz
        """

}

process glnexus {

    input:
        tuple val(pop_id), val(sample_ids), path(gvcfs)
        val snp_indel_caller
        val clair3_config

    output:
        tuple val(pop_id), val(sample_ids), path("snp_indel.bcf")

    script:
        if (snp_indel_caller == 'clair3')
        """
        printf '%s\\n' $gvcfs > gvcf_list.txt
        glnexus_cli --config $clair3_config --list gvcf_list.txt > snp_indel.bcf
        """
        else if (snp_indel_caller == 'deepvariant')
        """
        printf '%s\\n' $gvcfs > gvcf_list.txt
        glnexus_cli --config DeepVariant --list gvcf_list.txt > snp_indel.bcf
        """

    stub:
        """
        touch snp_indel.bcf
        """

}

process glnexus_post_processing {

    input:
        tuple val(pop_id), val(sample_ids), path(joint_snp_indel_bcf)

    output:
        tuple val(pop_id), path("snp_indel.raw.vcf.gz"), path("snp_indel.raw.vcf.gz.tbi")

    script:
        """
        # convert bcf to vcf, reorder since glnexus sorts samples alphabetically
        bcftools view $joint_snp_indel_bcf | bcftools view -s ${sample_ids.join(',')} | bgzip -@ ${task.cpus} -c > snp_indel.raw.vcf.gz
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
        tuple val(pop_id), path("snp_indel.split.vcf.gz"), path("snp_indel.split.vcf.gz.tbi")

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

process split_vcf {

    input:
        tuple val(pop_id), val(sample_id), path(joint_snp_indel_vcf), path(joint_snp_indel_vcf_index)

    output:
        tuple val(pop_id), val(sample_id), path("*snp_indel.vcf.gz"), path("*snp_indel.vcf.gz.tbi")

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

    publishDir "$outdir/$pop_id/$outdir2/phasing", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$snp_indel_caller.$filename"}, pattern: '*.read_list.txt'
    publishDir "$outdir/$pop_id/$outdir2/phasing", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$ref_name.$snp_indel_caller.$filename"}, pattern: '*.stats.gtf'

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
        tuple val(pop_id), path("snp_indel.phased.read_list.txt"), path("snp_indel.phased.stats.gtf")

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
        touch ${sample_id}.snp_indel.phased.vcf.gz
        touch ${sample_id}.snp_indel.phased.vcf.gz.tbi
        """

}

process merge_vcf {

    publishDir "$outdir/$pop_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$pop_id.$ref_name.$snp_indel_caller.$filename" }, pattern: 'snp_indel.phased.*'

    input:
        tuple val(pop_id), path(snp_indel_phased_vcfs), path(snp_indel_phased_vcf_indices)
        val outdir
        val outdir2
        val ref_name
        val snp_indel_caller

    output:
        tuple val(pop_id), path("snp_indel.phased.vcf.gz"), path("snp_indel.phased.vcf.gz.tbi")

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

process split_sv_vcfs {

    input:
        tuple val(pop_id), val(sv_caller), val(sample_ids), path(sv_vcfs), path(sv_vcf_indices)
        val ref_index
        val min_gap
        val chunks

    output:
        tuple val(pop_id), val(sv_caller), path("*.vcf")

    script:
        """
        # write sample/vcf mapping (sample_ids paired with staged vcf paths, in input order)
        paste <(printf '%s\\n' ${sample_ids.join(' ')}) <(printf '%s\\n' ${sv_vcfs.join(' ')}) > sample_vcfs.tsv

        # get chromosome lengths
        cut -f1,2 $ref_index > chrom_lengths.txt

        # extract DEL, DUP and INS positions
        while IFS=\$'\t' read -r sample vcf; do
            bcftools view -i 'SVTYPE="DEL"' \$vcf | bcftools query -f '%CHROM\\t%POS\\t%INFO/END\\n' >> all_del_pos.bed
            bcftools view -i 'SVTYPE="DUP"' \$vcf | bcftools query -f '%CHROM\\t%POS\\t%INFO/END\\n' >> all_dup_pos.bed
            bcftools view -i 'SVTYPE="INS"' \$vcf | bcftools query -f '%CHROM\\t%POS\\t%POS\\n' >> all_ins_pos.bed
        done < sample_vcfs.tsv

        # extract DEL and DUP gaps (complement of union)
        for i in del dup; do
            bedtools sort -i all_\${i}_pos.bed -g chrom_lengths.txt | bedtools merge -i - > \${i}_union.bed
            bedtools complement -i \${i}_union.bed -g chrom_lengths.txt | awk -v min=$min_gap '(\$3 - \$2) >= min' | bedtools sort -i - -g chrom_lengths.txt | bgzip > \${i}_gaps.bed.gz
            tabix -p bed \${i}_gaps.bed.gz
        done

        # extract INS gaps (regions between consecutive INS positions exceeding min gap)
        bedtools sort -i all_ins_pos.bed -g chrom_lengths.txt \
            | awk -v min=$min_gap '
                BEGIN { prev_chrom = ""; prev_end = 0 }
                {
                    chrom = \$1; start = \$2; end = \$3
                    if (chrom == prev_chrom && (start - prev_end) >= min)
                        print chrom "\t" prev_end "\t" start
                    prev_chrom = chrom
                    prev_end   = end
                }
            ' \
            | bgzip > ins_gaps.bed.gz
        tabix -p bed ins_gaps.bed.gz

        # extract INS+DUP gaps
        bedtools intersect -a <(zcat dup_gaps.bed.gz) -b <(zcat ins_gaps.bed.gz) | awk -v min=$min_gap '(\$3 - \$2) >= min' | bedtools sort -i - -g chrom_lengths.txt | bgzip > ins_dup_gaps.bed.gz
        tabix -p bed ins_dup_gaps.bed.gz

        # combine INS+DUP positions
        cat all_ins_pos.bed all_dup_pos.bed > all_ins_dup_pos.bed

        # extract chunk boundaries (pools small chunks)
        for i in del ins_dup; do
            compute_chunks.py --gaps \${i}_gaps.bed.gz --genome chrom_lengths.txt --positions all_\${i}_pos.bed --n-chunks $chunks --min-gap $min_gap --output chunk_boundaries_\${i}.bed
        done

        # extract BND and INV variants (they can span multiple chromosomes)
        while IFS=\$'\t' read -r sample vcf; do
            bcftools view -i 'INFO/SVTYPE="BND" || INFO/SVTYPE="INV"' \$vcf -o \$sample.${sv_caller}.BND_INV.vcf
        done < sample_vcfs.tsv

        # group chunk_boundaries rows by chunk name into per-chunk regions beds
        for i in del ins_dup; do
            mkdir -p chunk_regions_\${i}
            awk -v outdir=chunk_regions_\${i} '{print \$1"\t"\$2"\t"\$3 > outdir"/"\$4".bed"}' chunk_boundaries_\${i}.bed
        done

        # split DEL variants by chunk
        while IFS=\$'\t' read -r sample vcf; do
            for chunk_bed in chunk_regions_del/*.bed; do
                chunk_name=\$(basename "\$chunk_bed" .bed)
                bcftools view -R "\$chunk_bed" -i 'INFO/SVTYPE="DEL"' \$vcf -o \$sample.${sv_caller}.DEL_\${chunk_name}.vcf
            done
        done < sample_vcfs.tsv

        # split INS and DUP variants by chunk
        while IFS=\$'\t' read -r sample vcf; do
            for chunk_bed in chunk_regions_ins_dup/*.bed; do
                chunk_name=\$(basename "\$chunk_bed" .bed)
                bcftools view -R "\$chunk_bed" -i 'INFO/SVTYPE="INS" || INFO/SVTYPE="DUP"' \$vcf -o \$sample.${sv_caller}.INS_DUP_\${chunk_name}.vcf
            done
        done < sample_vcfs.tsv

        """

    stub:
        """
        touch chunk_boundaries_del.bed
        touch chunk_boundaries_ins_dup.bed
        for sample in ${sample_ids.join(' ')}; do
            touch \$sample.${sv_caller}.BND_INV.vcf
            touch \$sample.${sv_caller}.DEL_chunk_chr1_1.vcf
            touch \$sample.${sv_caller}.INS_DUP_chunk_chr1_1.vcf
        done
        """

}

process jasmine {

    input:
        tuple val(pop_id), val(partition), val(sv_caller), path(split_sv_vcfs), val(sample_ids), path(bams), path(bams_indices), val(data_type), val(related)
        val ref
        val ref_index

    output:
        tuple val(pop_id), val(sv_caller), path("${partition}.*.vcf.gz"), path("${partition}.*.vcf.gz.tbi"), optional: true

    script:
        def out_vcf = sv_caller == 'sniffles' ? 'sv.phased' : 'sv'
        // conditionally define iris arguments (by default, iris will pass minimap -x map-ont, the --pacbio flag passed to iris will pass minimap -x map-pb)
        def iris_args = '--run_iris iris_args=min_ins_length=20,--rerunracon,--keep_long_variants' + (data_type == 'pacbio' ? ',--pacbio' : '')
        // conditionally define require first flag
        def require_first_sample_optional = related == 'yes' ? '--require_first_sample' : ''
        """
        # create file lists in in_data_pipeface.csv row order (sample_ids and bams are pre-sorted by csv index)
        SAMPLES=(${sample_ids.join(' ')})
        BAMS=(${bams.join(' ')})
        for i in \${!SAMPLES[@]}; do
            realpath \${SAMPLES[\$i]}.${sv_caller}.${partition}.vcf >> vcfs.txt
            realpath \${BAMS[\$i]} >> bams.txt
        done
        # run jasmine
        # note. jasmine threads is specfically set to 1 due this issue: https://github.com/mkirsche/Jasmine/issues/49
        jasmine threads=1 out_dir=./ genome_file=$ref file_list=vcfs.txt bam_list=bams.txt out_file=${partition}.${out_vcf}.tmp.vcf min_support=1 --mark_specific spec_reads=7 spec_len=20 --pre_normalize --output_genotypes --clique_merging --dup_to_ins --normalize_type $require_first_sample_optional --default_zero_genotype $iris_args
        # fix vcf header (remove prefix to sample names that jasmine adds)
        grep '##' ${partition}.${out_vcf}.tmp.vcf > ${partition}.${out_vcf}.vcf
        grep '#CHROM' ${partition}.${out_vcf}.tmp.vcf | sed -E 's/\t[0-9]+_/\t/g' >> ${partition}.${out_vcf}.vcf
        grep -v '#' ${partition}.${out_vcf}.tmp.vcf >> ${partition}.${out_vcf}.vcf
        # add iris info tags if iris didn't process or refine any variants in current chunk
        if ! grep -q '##INFO=<ID=IRIS_PROCESSED' ${partition}.${out_vcf}.vcf; then
            sed -i '/^#CHROM/i ##INFO=<ID=IRIS_PROCESSED,Number=1,Type=String,Description="Whether or not a variant has been considered by Iris for refinement">' ${partition}.${out_vcf}.vcf
        fi
        if ! grep -q '##INFO=<ID=IRIS_REFINED' ${partition}.${out_vcf}.vcf; then
            sed -i '/^#CHROM/i ##INFO=<ID=IRIS_REFINED,Number=1,Type=String,Description="Whether or not a variant has been refined by Iris">' ${partition}.${out_vcf}.vcf
        fi
        # sort
        bcftools sort ${partition}.${out_vcf}.vcf -o ${partition}.${out_vcf}.vcf
        # compress and index vcf
        bgzip -@ ${task.cpus} ${partition}.${out_vcf}.vcf
        tabix ${partition}.${out_vcf}.vcf.gz
        """

    stub:
        def out_vcf = sv_caller == 'sniffles' ? 'sv.phased' : 'sv'
        """
        touch ${partition}.${out_vcf}.vcf.gz
        touch ${partition}.${out_vcf}.vcf.gz.tbi
        """

}

process concat_sv_vcf {

    publishDir "$outdir/$pop_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$pop_id.$ref_name.${sv_caller}.jasmine.$filename"}, pattern: 'sv.*vcf.gz*'

    input:
        tuple val(pop_id), val(sv_caller), path(sv_vcfs), path(sv_vcf_indices)
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(pop_id), val(sv_caller), path("sv.*vcf.gz")
        tuple val(pop_id), path("sv.*vcf.gz"), path("sv.*vcf.gz.tbi")

    script:
        // conditionally define output sv caller specific filename
        def out_vcf = sv_caller == 'sniffles' ? 'sv.phased' : 'sv'
        """
        # get list of vcfs
        VCFS=(*sv*.vcf.gz)
        # concat vcfs (or pass through if only one) and sort
        if [[ \${#VCFS[@]} -eq 1 ]]; then
            bcftools sort -Oz -o ${out_vcf}.vcf.gz \${VCFS[0]}
        else
            bcftools concat -a \${VCFS[@]} --threads ${task.cpus} | bcftools sort -Oz -o ${out_vcf}.vcf.gz -
        fi
        # index vcf
        tabix ${out_vcf}.vcf.gz
        """

    stub:
        def out_vcf = sv_caller == 'sniffles' ? 'sv.phased' : 'sv'
        """
        touch ${out_vcf}.vcf.gz
        touch ${out_vcf}.vcf.gz.tbi
        """

}

process somalier {

    publishDir "$outdir/$pop_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$pop_id.$ref_name.$filename" }, pattern: 'somalier*'

    input:
        tuple val(pop_id), path(somalier_files)
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(pop_id), path("somalier.samples.tsv"), path("somalier.pairs.tsv"), path("somalier.html")

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
        tuple val(pop_id), val(sample_ids), path(bams), path(bam_indices), val(data_type), path(split_bed)
        val ref
        val ref_index

    output:
        tuple val(pop_id), val(sample_ids), path("tr.*.vcf.gz"), path("tr.*.vcf.gz.tbi")

    script:
        def bams_csv = bams.join(',')
        def samples_csv = sample_ids.join(',')
        // define alignment parameters for ONT to account for higher incidence of indels in homopolymers (defaults are tailored to pacbio hifi)
        def alignment_params_optional = data_type == 'ont' ? "--alignment-params -1.0,-0.458675,-1.0,-0.458675,-0.00005800168,-1,-1" : ''
        """
        # run longtr
        ID=\$(echo $split_bed | sed 's/split.//;s/.bed//')
        LongTR --bams $bams_csv --bam-samps $samples_csv --bam-libs $samples_csv --fasta $ref --regions $split_bed --tr-vcf tr.\${ID}.vcf.gz --phased-bam --output-gls --output-pls --output-phased-gls --output-filter $alignment_params_optional --log longtr.log
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

    publishDir "$outdir/$pop_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$pop_id.${ref_name}.longtr.$filename" }, pattern: 'tr.vcf.gz*'

    input:
        tuple val(pop_id), val(sample_ids), path(tr_vcfs), path(tr_vcf_indices)
        val outdir
        val outdir2
        val ref_name

    output:
        tuple val(pop_id), path("tr.vcf.gz"), path("tr.vcf.gz.tbi")

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
        tuple val(pop_id), path("snp_indel.phased.annotated.vcf.gz"), path("snp_indel.phased.annotated.vcf.gz.tbi")

    script:
        """
        # run vep
        vep -i $joint_snp_indel_phased_vcf -o snp_indel.phased.annotated.vcf.gz --format vcf --vcf --fasta $ref --dir $vep_db --assembly GRCh38 --species homo_sapiens --cache --offline --merged \
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

process vep_sv {

    publishDir "$outdir/$pop_id/$outdir2", mode: params.publish_mode, overwrite: true, saveAs: { filename -> "$pop_id.$ref_name.${sv_caller}.jasmine.$filename"}, pattern: '*.annotated.vcf.gz*'

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
        tuple val(pop_id), path("*.annotated.vcf.gz"), path("*.annotated.vcf.gz.tbi")

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
    ref = "${params.ref}".trim()
    ref_index = "${params.ref_index}".trim()
    tr_calling = "${params.tr_calling}".trim()
    tr_call_regions = "${params.tr_call_regions}".trim()
    snp_indel_caller = "${params.snp_indel_caller}".trim()
    clair3_config = "${params.clair3_config}".trim()
    annotate = "${params.annotate}".trim()
    annotate_override = "${params.annotate_override}".trim()
    outdir = "${params.outdir}".trim()
    outdir2 = "${params.outdir2}".trim()
    min_gap = "${params.min_gap}".trim()
    chunks = "${params.chunks}".trim()
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
    [snp_indel_caller: snp_indel_caller, clair3_config: clair3_config].each { param, val ->
        if (!val) {
            exit 1, "No value provided for '${param}'. Set to 'NONE' if not required."
        }
    }
    // check file existence
    def required_files = [in_data: in_data, ref: ref, ref_index: ref_index]
    def optional_files = [tr_call_regions: tr_call_regions]
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
    // check parameter constraints
    if (!(snp_indel_caller in ['clair3', 'deepvariant', 'NONE'])) {
        exit 1, "SNP/indel caller should be 'clair3', 'deepvariant', or 'NONE', snp_indel_caller = '${snp_indel_caller}' provided."
    }
    if (!(annotate in ['yes', 'no'])) {
        exit 1, "'annotate' should be either 'yes' or 'no', annotate = '${annotate}' provided."
    }
    if (annotate == 'yes') {
        if (!ref.toLowerCase().contains('hg38') && !ref.toLowerCase().contains('grch38') && annotate_override != 'yes') {
            exit 1, "Only hg38/GRCh38 is supported for annotation. It looks like you may not be passing a hg38/GRCh38 reference genome based on the filename of the reference genome. ref = '${ref}' provided. Pass '--annotate_override yes' on the command line to override this error."
        }
    }
    if (!(tr_calling in ['yes', 'no'])) {
        exit 1, "'tr_calling' should be either 'yes' or 'no', tr_calling = '${tr_calling}' provided."
    }
    if (!(params.publish_mode in ['copy', 'copyNoFollow', 'link', 'move', 'rellink', 'symlink'])) {
        exit 1, "Choice of publishing mode should be 'copy', 'copyNoFollow', 'link', 'move', 'rellink' or 'symlink', publish_mode = '$params.publish_mode' provided."
    }
    // check cross-parameter compatibility
    if (snp_indel_caller == 'clair3') {
        if (clair3_config == 'NONE') {
            exit 1, "When clair3 is selected as the SNP/indel calling software, provide a path to an appropriate clair3 config file for GLnexus (clair3_config), clair3_config = 'NONE' provided."
        }
        if (!file(clair3_config).exists()) {
            exit 1, "File does not exist, clair3_config = '${clair3_config}' provided."
        }
    }
    if (tr_calling == 'yes' && tr_call_regions == 'NONE') {
        exit 1, "When calling tandem repeats, provide a valid tandem repeat call regions file, tr_calling = '${tr_calling}' and tr_call_regions = 'NONE' provided."
    }
    else if (tr_calling == 'no' && tr_call_regions != 'NONE') {
        exit 1, "When not calling tandem repeats, set tandem repeat call regions file to 'NONE', tr_calling = '${tr_calling}' and tr_call_regions = '${tr_call_regions}' provided."
    }

    // read in data
    Channel
        .fromPath(in_data)
        .splitCsv(header: true, sep: ',', strip: true)
        .toList()
        .flatMap { rows -> rows.withIndex().collect { row, idx -> [row, idx] } }
        .multiMap { row, index ->
            in_data:            tuple(row.pop_id, row.sample_id, row.related, row.family_position, row.gvcf, row.bam, row.sniffles, row.cutesv, row.somalier, row.data_type)
            ids:                tuple(row.pop_id, row.sample_id)
            gvcfs:              tuple(row.pop_id, row.sample_id, row.gvcf, index)
            gvcfs_bams:         tuple(row.pop_id, row.sample_id, row.gvcf, row.bam)
            svs:                tuple(row.pop_id, row.sample_id, row.sniffles, row.cutesv)
            somaliers:          tuple(row.pop_id, row.somalier)
            bams_data_type:     tuple(row.pop_id, row.sample_id, row.bam, row.data_type, index)
            related:            tuple(row.pop_id, row.related)
            row_validation:     tuple(row.pop_id, row.sample_id, row.gvcf, row.bam, row.sniffles, row.cutesv, row.somalier, row.data_type, row.related, row.family_position)
            cohort_validation:  tuple(row.pop_id, row.sample_id, row.related, row.family_position, row.gvcf, row.bam, row.sniffles, row.cutesv, row.somalier, row.data_type)
            dup_validation:     row.sample_id
            tr_validation:      row.bam
        }
        .set { csv }

    // build channels
    csv.in_data
        .groupTuple(by: [0,2,9])
        .set { in_data_ch }

    csv.ids
        .set { id_ch }

    csv.gvcfs
        .filter { pop_id, sample_id, gvcf, index -> gvcf != 'NONE' }
        .set { gvcfs_ch }

    csv.gvcfs_bams
        .filter { pop_id, sample_id, gvcf, bam -> gvcf != 'NONE' }
        .map { pop_id, sample_id, gvcf, bam -> tuple(pop_id, sample_id, bam, "${bam}.bai") }
        .set { gvcfs_bams_ch }

    csv.svs
        .flatMap { pop_id, sample_id, sniffles, cutesv ->
            def tuples = []
            if (sniffles != 'NONE') tuples << tuple(pop_id, sample_id, 'sniffles', sniffles, "${sniffles}.tbi")
            if (cutesv != 'NONE') tuples << tuple(pop_id, sample_id, 'cutesv', cutesv, "${cutesv}.tbi")
            return tuples
        }
        .set { svs_ch }

    csv.somaliers
        .filter { pop_id, somalier -> somalier != 'NONE' }
        .groupTuple(by: 0)
        .set { somalier_files_ch }

    csv.bams_data_type
        .filter { pop_id, sample_id, bam, data_type, index -> bam != 'NONE' }
        .map { pop_id, sample_id, bam, data_type, index -> tuple(pop_id, sample_id, bam, "${bam}.bai", data_type, index) }
        .groupTuple(by: [0,4])
        .map { pop_id, sample_ids, bams, bais, data_type, indices ->
            def sorted = [indices, sample_ids, bams, bais].transpose().sort { a, b -> a[0] <=> b[0] }
            tuple(pop_id, sorted.collect { it[1] }, sorted.collect { it[2] }, sorted.collect { it[3] }, data_type)
        }
        .set { bams_data_type_ch }

    csv.related
        .unique()
        .set { related_ch }

    // check user provided in data
    csv.row_validation
        .map { pop_id, sample_id, gvcf, bam, sniffles, cutesv, somalier, data_type, related, family_position ->
            // check for empty entries
            def required_cols = [pop_id: pop_id, sample_id: sample_id, data_type: data_type, related: related]
            def optional_cols = [gvcf: gvcf, bam: bam, sniffles: sniffles, cutesv: cutesv, somalier: somalier, family_position: family_position]
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
            // check file existence
            [gvcf: gvcf, bam: bam, sniffles: sniffles, cutesv: cutesv, somalier: somalier].each { col, val ->
                if (val != 'NONE' && !file(val).exists()) {
                    exit 1, "There is an entry in the '${col}' column of '${in_data}' which doesn't exist. Check file '${val}'. Set to 'NONE' if not required."
                }
            }
            [bam: [bam, "${bam}.bai"], sniffles: [sniffles, "${sniffles}.tbi"], cutesv: [cutesv, "${cutesv}.tbi"]].each { col, paths ->
                def (val, idx) = paths
                if (val != 'NONE' && !file(idx).exists()) {
                    exit 1, "There is an entry in the '${col}' column of '${in_data}' which doesn't look to have an associated index. Expecting '${idx}'."
                }
            }
            // check value constraints
            [pop_id: pop_id, sample_id: sample_id].each { col, val ->
                if (val == 'NONE') {
                    exit 1, "Entries in the '${col}' column of '${in_data}' should not be 'NONE'."
                }
            }
            if (!(related in ['yes', 'no'])) {
                exit 1, "Entries in the 'related' column of '${in_data}' should be 'yes' or 'no', '${related}' provided."
            }
            if (!(data_type in ['ont', 'pacbio'])) {
                exit 1, "Entries in the 'data_type' column of '${in_data}' should be 'ont' or 'pacbio', '${data_type}' provided."
            }
            // check cross-column compatibility
            if (related == 'no' && family_position != 'NONE') {
                exit 1, "When 'related' is 'no', 'family_position' should be 'NONE', '${family_position}' provided."
            }
            if (related == 'yes' && family_position == 'NONE') {
                exit 1, "When 'related' is 'yes', 'family_position' should not be 'NONE'."
            }
            if (gvcf != 'NONE' && bam == 'NONE') {
                exit 1, "When a GVCF file is provided in the 'gvcf' column of '${in_data}', the associated BAM file must be provided in the 'bam' column. gvcf = '${gvcf}' and bam = '${bam}' provided."
            }
            if (sniffles != 'NONE' && bam == 'NONE') {
                exit 1, "When a Sniffles SV VCF is provided in the 'sniffles' column of '${in_data}', the associated BAM file must be provided in the 'bam' column. sniffles = '${sniffles}' and bam = '${bam}' provided."
            }
            if (cutesv != 'NONE' && bam == 'NONE') {
                exit 1, "When a cuteSV VCF is provided in the 'cutesv' column of '${in_data}', the associated BAM file must be provided in the 'bam' column. cutesv = '${cutesv}' and bam = '${bam}' provided."
            }
        }

    csv.cohort_validation
        .groupTuple(by: 0)
        .map { pop_id, sample_ids, relateds, family_positions, gvcfs, bams, sniffless, cutesvs, somaliers, data_types ->
            if (sample_ids.unique().size() < 2) {
                exit 1, "Entries in the 'sample_id' column of '$in_data' should contain 2 or more unique values for every 'pop_id', '$sample_ids' provided for population '$pop_id'."
            }
            [gvcf: gvcfs, bam: bams, sniffles: sniffless, cutesv: cutesvs, somalier: somaliers].each { col, vals ->
                if (vals.any { it == 'NONE' } && !vals.every { it == 'NONE' }) {
                    exit 1, "Entries in the '${col}' column of '$in_data' should be either all 'NONE' or all real files for a given 'pop_id', '$vals' provided for population '$pop_id'."
                }
            }
            if (gvcfs.every { it == 'NONE' } && snp_indel_caller != 'NONE') {
                exit 1, "All entries in the 'gvcf' column of '$in_data' are 'NONE' for population '$pop_id', but snp_indel_caller = '${snp_indel_caller}' is set. Set snp_indel_caller to 'NONE' if SNP/indel merging is not required."
            }
            if (!gvcfs.every { it == 'NONE' } && snp_indel_caller == 'NONE') {
                exit 1, "gVCF files are provided in the 'gvcf' column of '$in_data' for population '$pop_id', but snp_indel_caller = 'NONE'. Set snp_indel_caller to 'clair3' or 'deepvariant'."
            }
            [data_type: data_types, related: relateds].each { col, vals ->
                if (vals.unique().size() != 1) {
                    exit 1, "Entries in the '${col}' column of '$in_data' should be identical for a given 'pop_id'. '$vals' provided for population '$pop_id'."
                }
            }
            if (relateds[0] == 'yes') {
                def proband_count = family_positions.count { it == 'proband' }
                if (proband_count != 1) {
                    exit 1, "When 'related' is 'yes', one sample should have 'family_position' = 'proband' for each 'pop_id'. ${proband_count} proband(s) provided for population '$pop_id'."
                }
                if (family_positions[0] != 'proband') {
                    exit 1, "When 'related' is 'yes', the proband must be listed first in '${in_data}' for each 'pop_id'. The first sample for population '$pop_id' has 'family_position' = '${family_positions[0]}'."
                }
            }
        }

    csv.dup_validation
        .collect()
        .map { sample_ids ->
            def duplicates = sample_ids.findAll { id -> sample_ids.count(id) > 1 }.toSet()
            if (duplicates) {
                exit 1, "Entries in the 'sample_id' column of '$in_data' should all be unique values, duplicates: '$duplicates'."
            }
        }

    if (tr_calling == 'yes') {
        csv.tr_validation
            .collect()
            .map { bams ->
                def any_real_bam = bams.any { it != 'NONE' }
                if (!any_real_bam) {
                    exit 1, "When tandem repeat calling is turned on, relevant BAM files should be provided in the 'bam' column of '${in_data}', tr_calling = '${tr_calling}' and '${bams}' provided."
                }
            }
    }

    // helpers
    // sort gvcfs by csv row index to preserve in data sample order
    def sort_by_index = { pop_id, indices, sample_ids, gvcfs_list ->
        def sorted = [indices, sample_ids, gvcfs_list].transpose().sort { a, b -> a[0] <=> b[0] }
        tuple(pop_id, sorted.collect { it[1] }, sorted.collect { it[2] })
    }

    // workflow
    // pre process
    scrape_settings(in_data_ch, popface_version, in_data, ref, ref_index, snp_indel_caller, tr_calling, tr_call_regions, annotate, outdir, outdir2)
    // gvcf merging
    if (snp_indel_caller != 'NONE') {
        if (snp_indel_caller == 'clair3') {
            gvcfs_for_pre = gvcfs_ch.map { pop_id, sample_id, gvcf, index -> tuple(pop_id, sample_id, gvcf) }
            gvcfs = glnexus_pre_processing(gvcfs_for_pre, ref_index)
                .join(gvcfs_ch.map { pop_id, sample_id, gvcf, index -> tuple(pop_id, sample_id, index) }, by: [0,1])
                .map { pop_id, sample_id, amended_gvcf, index -> tuple(pop_id, index as Integer, sample_id, amended_gvcf) }
                .groupTuple(by: 0)
                .map(sort_by_index)
        }
        else if (snp_indel_caller == 'deepvariant') {
            gvcfs = gvcfs_ch
                .map { pop_id, sample_id, gvcf, index -> tuple(pop_id, index as Integer, sample_id, gvcf) }
                .groupTuple(by: 0)
                .map(sort_by_index)
        }
        joint_snp_indel_bcf = glnexus(gvcfs, snp_indel_caller, clair3_config)
        joint_snp_indel_vcf = glnexus_post_processing(joint_snp_indel_bcf)
        // split multiallelic variants
        joint_snp_indel_split_vcf = split_multiallele(joint_snp_indel_vcf, ref, ref_index)
        // phasing
        joint_snp_indel_vcf_id = joint_snp_indel_split_vcf
            .combine(id_ch)
            .map { pop_id, joint_vcf, joint_vcf_index, pop_id2, sample_id ->
                if (pop_id != pop_id2) {
                    return null
                }
                tuple(pop_id, sample_id, joint_vcf, joint_vcf_index)
            }
        snp_indel_vcf = split_vcf(joint_snp_indel_vcf_id)
        (snp_indel_phased_vcfs, stats) = whatshap_phase(snp_indel_vcf.join(gvcfs_bams_ch, by: [0,1]), ref, ref_index, outdir, outdir2, ref_name, snp_indel_caller)
        vcfs = snp_indel_phased_vcfs.groupTuple(by: 0)
        joint_snp_indel_phased_vcf = merge_vcf(vcfs, outdir, outdir2, ref_name, snp_indel_caller)
    }
    // sv vcf merging
    sv_vcfs_grouped = svs_ch
        .map { pop_id, sample_id, sv_caller, vcf, tbi -> tuple(pop_id, sv_caller, sample_id, vcf, tbi) }
        .groupTuple(by: [0,1])
    split_sv_vcfs_ch = split_sv_vcfs(sv_vcfs_grouped, ref_index, min_gap, chunks)
    jasmine_input = split_sv_vcfs_ch
        .flatMap { pop_id, sv_caller, vcfs ->
            vcfs.collectMany { vcf ->
                def partition = vcf.baseName.tokenize('.')[-1]
                [ tuple(pop_id, partition, sv_caller, vcf) ]
            }
        }
        .groupTuple(by: [0,1,2])
        .combine(bams_data_type_ch, by: 0)
        .combine(related_ch, by: 0)
    merged_sv_vcfs = jasmine(jasmine_input, ref, ref_index)
    (joint_sv_vcf, joint_sv_vcf_indexed) = concat_sv_vcf(merged_sv_vcfs.groupTuple(by: [0,1]), outdir, outdir2, ref_name)
    // somalier
    somalier(somalier_files_ch, outdir, outdir2, ref_name)
    // tr calling
    if (tr_calling == 'yes') {
        split_bed = longtr_pre_processing(tr_call_regions).flatten()
        longtr_input = bams_data_type_ch
            .combine(split_bed)
        tr_vcfs = longtr(longtr_input, ref, ref_index)
        concat_tr_vcf(tr_vcfs.groupTuple(by: [0,1]), outdir, outdir2, ref_name)
    }
    // annotation
    if (annotate == 'yes') {
        if (snp_indel_caller != 'NONE') {
            vep_snp_indel(joint_snp_indel_phased_vcf, ref, ref_index, vep_db, revel_db, gnomad_db, clinvar_db, cadd_snv_db, cadd_indel_db, spliceai_snv_db, spliceai_indel_db, alphamissense_db, outdir, outdir2, ref_name, snp_indel_caller)
        }
        vep_sv(joint_sv_vcf, ref, ref_index, vep_db, gnomad_db, cadd_sv_db, outdir, outdir2, ref_name)
    }
}
