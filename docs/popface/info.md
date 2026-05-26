# Popface information

## Workflow

```mermaid
%%{init:{'theme':'dark','themeVariables':{'fontSize':'11px'}}}%%
flowchart TB

input_gvcf("DeepVariant gVCFs")
input_bam("Aligned BAMs")
input_svs("Structural Variant VCF's")
gvcf_merging("gVCF merging")
joint_split_multiallele("Split multiallelic variants into biallelic variants")
split_vcf("Split joint VCF")
snp_indel_phasing("SNP/indel phasing")
merge_vcf("Merge VCF")
joint_snp_indel_annotation("Joint SNP/indel annotation (hg38 only)")
split_sv_vcf("Split SV VCF")
sv_vcf_merging("Structural variant VCF merging")
concat_sv_vcf("Concatenate SV VCFs")
joint_sv_annotation("Joint structural variant annotation (hg38 only)")
joint_tr_calling("Joint TR calling")
concat_tr_vcf("Concatenate TR VCFs")
input_somalier("Somalier extracted files")
joint_somalier("Joint somalier relatedness/quality control check")

input_gvcf-.->gvcf_merging-.->joint_split_multiallele-.->split_vcf-.->snp_indel_phasing-.->merge_vcf-.->concat_sv_vcf-.->joint_snp_indel_annotation
input_bam-.->snp_indel_phasing
input_svs-.->split_sv_vcf-.->sv_vcf_merging-.->joint_sv_annotation
input_bam-.->sv_vcf_merging
input_bam-.->joint_tr_calling-.->concat_tr_vcf
input_somalier-.->joint_somalier

```

## Main analyses

- Quads to populations (up to 1000 individuals)

## Main tools

- [GLnexus](https://github.com/dnanexus-rnd/GLnexus)
- [WhatsHap](https://github.com/whatshap/whatshap)
- [Jasmine (customised)](https://github.com/bioinfomethods/Jasmine)
- [somalier](https://github.com/brentp/somalier)
- [LongTR](https://github.com/gymrek-lab/LongTR)
- [ensembl-vep](https://github.com/Ensembl/ensembl-vep)

*[See the list of software and their versions used by this version of popface](../software_versions.txt) as well as the [list of variant databases and their versions](../database_versions.txt) if variant annotation is carried out (assuming the default [nextflow_popface.config](../../config/nextflow_popface.config) file is used).*

## Main input files

### Required

- Indexed reference genome

### Optional

- Clair3 or DeepVariant gVCF files
- Aligned BAM files
- Sniffles and/or cuteSV SV VCF files
- Somalier extracted files

## Main output files

- Joint phased Clair3 or DeepVariant SNP/indel VCF file
- Joint phased and annotated Clair3 or DeepVariant SNP/indel VCF file (hg38 only)
- Joint phased Sniffles2 and/or un-phased cuteSV SV VCF file
- Joint phased and annotated Sniffles2 and/or un-phased and annotated cuteSV SV VCF file (hg38 only)
- Joint phased tandem repeat VCF file
- Joint relatedness and quality control somalier TSV and HTML files

