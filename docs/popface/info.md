# Popface information

## Workflow

```mermaid
%%{init: {'theme':'dark'}}%%
flowchart LR

input_gvcf("Input data: <br><br> DeepVariant gVCF's")
input_bam("Input data: <br><br> Aligned BAM's")
gvcf_merging{{"gVCF merging"}}
joint_split_multiallele{{"Split multiallelic variants into biallelic variants"}}
split_vcf{{"Split joint VCF"}}
snp_indel_phasing{{"SNP/indel phasing"}}
merge_vcf{{"Merge VCF"}}
joint_snp_indel_annotation{{"Joint SNP/indel annotation (hg38 only)"}}
tr_calling{{"TR calling"}}
concat_vcf{{"Concatenate TR VCF's"}}

input_somalier("Input data: <br><br> Somalier extracted files")
joint_somalier{{"Joint somalier relatedness/quality control check"}}

input_gvcf-.->gvcf_merging-.->joint_split_multiallele-.->split_vcf-.->snp_indel_phasing-.->merge_vcf-.->joint_snp_indel_annotation
input_bam-.->snp_indel_phasing
input_bam-.->tr_calling-.->concat_vcf
input_somalier-.->joint_somalier

```

## Main analyses

- Quads to populations

## Main tools

- [GLnexus](https://github.com/dnanexus-rnd/GLnexus)
- [WhatsHap](https://github.com/whatshap/whatshap)
- [somalier](https://github.com/brentp/somalier)
- [Samtools](https://github.com/samtools/samtools)
- [ensembl-vep](https://github.com/Ensembl/ensembl-vep)

*[See the list of software and their versions used by this version of popface](../software_versions.txt) as well as the [list of variant databases and their versions](../database_versions.txt) if variant annotation is carried out (assuming the default [nextflow_popface.config](../../config/nextflow_popface.config) file is used).*

## Main input files

### Required

- Indexed reference genome

### Optional

- DeepVariant gVCF files and aligned BAM files
- Somalier extracted files

## Main output files

- Joint phased DeepVariant SNP/indel VCF file
- Joint phased and annotated DeepVariant SNP/indel VCF file (hg38 only)
- Joint relatedness and quality control somalier TSV and HTML files

