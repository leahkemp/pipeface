# Pipeface information

## Workflow

### Singleton

```mermaid
%%{init:{'theme':'dark','themeVariables':{'fontSize':'11'},'flowchart':{'nodeSpacing':20}}}%%
flowchart TB

input_data("ONT fastq.gz <br> and/or <br> ONT fastq <br> and/or <br> ONT uBAM <br> and/or <br> pacbio HiFi uBAM")
alignment("Merge runs (if needed), bam to fastq conversion (if needed), alignment, sorting")
depth("Calculate alignment depth")
snp_indel_calling("SNP/indel variant calling")
somatic_calling("Somatic SNV/indel variant calling")
split_multiallele("Split multiallelic variants into biallelic variants")
snp_indel_phasing("SNP/indel phasing")
snp_indel_annotation("SNP/indel annotation (hg38 only)")
haplotagging("Haplotagging bams")
calculate_base_mod_freqs("Calculate base modification frequencies (uBAMs containing base modifications only)")
somalier("Somalier extract")
tr_calling("TR calling")
sv_calling("Structural variant calling")
sv_annotation("Structural variant annotation (hg38 only)")
sv_repeat_annotation("Structural variant repeat annotation (hg38 only)")
puzzleapp_preprocessing("Puzzleapp preprocessing (hg38 only)")

input_data-.->alignment-.->snp_indel_calling-.->split_multiallele-.->snp_indel_phasing-.->haplotagging-.->sv_calling
alignment-.->depth
alignment-.->somatic_calling
alignment-.->haplotagging
haplotagging-.->calculate_base_mod_freqs
haplotagging-.->somalier
haplotagging-.->tr_calling
snp_indel_phasing-.->snp_indel_annotation-.->puzzleapp_preprocessing
sv_calling-.->sv_annotation-.->sv_repeat_annotation-.->puzzleapp_preprocessing
depth-.->puzzleapp_preprocessing

```

### Duo

```mermaid
%%{init:{'theme':'dark','themeVariables':{'fontSize':'11'},'flowchart':{'nodeSpacing':20}}}%%
flowchart TB

input_data("ONT fastq.gz <br> and/or <br> ONT fastq <br> and/or <br> ONT uBAM <br> and/or <br> pacbio HiFi uBAM")
alignment("Merge runs (if needed), bam to fastq conversion (if needed), alignment, sorting")
depth("Calculate alignment depth")
snp_indel_calling("SNP/indel variant calling")
somatic_calling("Somatic SNV/indel variant calling")
split_multiallele("Split multiallelic variants into biallelic variants")
snp_indel_phasing("SNP/indel phasing")
joint_somalier("Joint somalier relatedness/quality control check")
gvcf_merging("gVCF merging")
joint_split_multiallele("Split multiallelic variants into biallelic variants")
joint_snp_indel_phasing("Joint SNP/indel phasing")
joint_snp_indel_annotation("Joint SNP/indel annotation (hg38 only)")
haplotagging("Haplotagging bams")
calculate_base_mod_freqs("Calculate base modification frequencies (uBAMs containing base modifications only)")
tr_calling("TR calling")
joint_tr_calling("Joint TR calling")
sv_calling("Structural variant calling")
sv_vcf_merging("Structural variant VCF merging")
joint_sv_annotation("Joint structural variant annotation (hg38 only)")
sv_repeat_annotation("Structural variant repeat annotation (hg38 only)")
puzzleapp_preprocessing("Puzzleapp preprocessing (hg38 only)")

input_data-.->alignment-.->snp_indel_calling-.->split_multiallele-.->snp_indel_phasing-.->haplotagging-.->sv_calling
alignment-.->depth
alignment-.->somatic_calling
alignment-.->haplotagging
haplotagging-.->calculate_base_mod_freqs
haplotagging-.->tr_calling
haplotagging-.->joint_tr_calling
haplotagging-.->joint_somalier
snp_indel_calling-.->gvcf_merging-.->joint_split_multiallele-.->joint_snp_indel_phasing-.->joint_snp_indel_annotation-.->puzzleapp_preprocessing
sv_calling-.->sv_vcf_merging-.->joint_sv_annotation-.->sv_repeat_annotation-.->puzzleapp_preprocessing
depth-.->puzzleapp_preprocessing

```

### Trio

```mermaid
%%{init:{'theme':'dark','themeVariables':{'fontSize':'11'},'flowchart':{'nodeSpacing':20}}}%%
flowchart TB

input_data("ONT fastq.gz <br> and/or <br> ONT fastq <br> and/or <br> ONT uBAM <br> and/or <br> pacbio HiFi uBAM")
alignment("Merge runs (if needed), bam to fastq conversion (if needed), alignment, sorting")
depth("Calculate alignment depth")
snp_indel_calling("SNP/indel variant calling")
somatic_calling("Somatic SNV/indel variant calling")
split_multiallele("Split multiallelic variants into biallelic variants")
snp_indel_phasing("SNP/indel phasing")
joint_snp_indel_calling("Joint SNP/indel variant calling")
joint_somalier("Joint somalier relatedness/quality control check")
gvcf_merging("gVCF merging")
joint_split_multiallele("Split multiallelic variants into biallelic variants")
joint_snp_indel_phasing("Joint SNP/indel phasing")
joint_snp_indel_annotation("Joint SNP/indel annotation (hg38 only)")
haplotagging("Haplotagging bams")
calculate_base_mod_freqs("Calculate base modification frequencies (uBAMs containing base modifications only)")
tr_calling("TR calling")
joint_tr_calling("Joint TR calling")
sv_calling("Structural variant calling")
sv_vcf_merging("Structural variant VCF merging")
joint_sv_annotation("Joint structural variant annotation (hg38 only)")
sv_repeat_annotation("Structural variant repeat annotation (hg38 only)")
puzzleapp_preprocessing("Puzzleapp preprocessing (hg38 only)")

input_data-.->alignment-.->snp_indel_calling-.->split_multiallele-.->snp_indel_phasing-.->haplotagging-.->sv_calling
alignment-.->depth
alignment-.->somatic_calling
alignment-.->haplotagging
haplotagging-.->calculate_base_mod_freqs
haplotagging-.->tr_calling
haplotagging-.->joint_tr_calling
haplotagging-.->joint_somalier
snp_indel_phasing-.->joint_snp_indel_calling-.->gvcf_merging-.->joint_split_multiallele-.->joint_snp_indel_phasing-.->joint_snp_indel_annotation-.->puzzleapp_preprocessing
sv_calling-.->sv_vcf_merging-.->joint_sv_annotation-.->sv_repeat_annotation-.->puzzleapp_preprocessing
depth-.->puzzleapp_preprocessing

```

## Main analyses

- ONT and/or pacbio HiFi data
- Singletons, duos or trios
- WGS and/or targeted
- hg38 or chm13 reference genome

## Main tools

- [Minimap2](https://github.com/lh3/minimap2)
- [Clair3](https://github.com/HKU-BAL/Clair3) or [DeepVariant](https://github.com/google/deepvariant)/[DeepTrio](https://github.com/google/deepvariant/blob/r1.10/docs/deeptrio-details.md)
- [WhatsHap](https://github.com/whatshap/whatshap)
- [GLnexus](https://github.com/dnanexus-rnd/GLnexus)
- [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) and/or [cuteSV](https://github.com/tjiangHIT/cuteSV)
- [Jasmine (customised)](https://github.com/bioinfomethods/Jasmine)
- [somalier](https://github.com/brentp/somalier)
- [mosdepth](https://github.com/brentp/mosdepth)
- [minimod](https://github.com/warp9seq/minimod?tab=readme-ov-file)
- [LongTR](https://github.com/gymrek-lab/LongTR)
- [ensembl-vep](https://github.com/Ensembl/ensembl-vep)
- [SVscanner](https://github.com/GenTechGp/SVscanner)
- [puzzleapp](https://github.com/GenTechGp/puzzleapp)
- [ClairS-TO](https://github.com/HKU-BAL/ClairS-TO)

## Main annotation databases

Used when variant annotation is turned on (hg38 only):

- [VEP cache](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html) (merged Ensembl/RefSeq)
- [REVEL](https://sites.google.com/site/revelgenomics/)
- [gnomAD](https://gnomad.broadinstitute.org/)
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
- [CADD](https://cadd.gs.washington.edu/) (SNVs and indels) and [CADD-SV](https://cadd-sv.bihealth.org/)
- [SpliceAI](https://github.com/Illumina/SpliceAI)
- [AlphaMissense](https://github.com/google-deepmind/alphamissense)
- [Dfam](https://www.dfam.org/) (repeat annotation of SVs by SVscanner)
- [STRchive](https://strchive.org/) (bundled with SVscanner)

*[See the list of software and their versions used by this version of pipeface](../software_versions.txt) as well as the [list of variant databases and their versions](../database_versions.txt) if variant annotation is carried out (assuming the default [nextflow_pipeface.config](../../config/nextflow_pipeface.config) file is used).*

## Main input files

### Required

- ONT/pacbio HiFi FASTQ (gzipped or uncompressed) or unaligned BAM
- Indexed reference genome
- Clair3 models (if running Clair3)

### Optional

- Regions of interest BED file
- Tandem repeat BED file
- PAR regions BED file (if running in haploid aware mode)
- Tandem repeat calling regions (if running tandem repeat calling)
- Somalier sites file (if running relatedness check)

## Main output files

### Singleton

- Aligned, sorted and haplotagged bam
- Alignment depth per chromosome (and per region in the case of targeted sequencing)
- Phased Clair3 or DeepVariant SNP/indel VCF file
- Phased and annotated Clair3 or DeepVariant SNP/indel VCF file (hg38 only)
- Clair3 or DeepVariant SNP/indel gVCF file
- Bed and bigwig base modification frequencies for complete read set and separate haplotypes (uBAMs containing base modifications only)
- Phased tandem repeat VCF file
- Phased Sniffles2 and/or un-phased cuteSV SV VCF file
- Phased and annotated (VEP + SVscanner) Sniffles2 and/or un-phased and annotated cuteSV SV VCF file (hg38 only)
- SVscanner text diagrams of the repeat elements annotated in each SV (hg38 only)
- Puzzleapp SNP/indel and SV TSV files and coverage/VAF quality control HTML file (hg38 only)
- Somalier extracted files
- ClairS-TO somatic SNV/indel VCF files

### Duo

- Aligned, sorted and haplotagged bam
- Alignment depth per chromosome (and per region in the case of targeted sequencing)
- DeepVariant SNP/indel gVCF file
- Joint phased DeepVariant SNP/indel VCF file
- Joint phased and annotated DeepVariant SNP/indel VCF file (hg38 only)
- Bed and bigwig base modification frequencies for complete read set and separate haplotypes (uBAMs containing base modifications only)
- Phased tandem repeat VCF file
- Joint phased Sniffles2 and/or un-phased cuteSV SV VCF file
- Joint phased and annotated (VEP + SVscanner) Sniffles2 and/or un-phased and annotated cuteSV SV VCF file (hg38 only)
- SVscanner text diagrams of the repeat elements annotated in each joint SV (hg38 only)
- Puzzleapp joint SNP/indel and SV TSV files and coverage/VAF quality control HTML file (hg38 only)
- Joint phased tandem repeat VCF file
- Somalier extracted files
- Joint relatedness and quality control somalier TSV and HTML files
- ClairS-TO somatic SNV/indel VCF files

### Trio

- Aligned, sorted and haplotagged bam
- Alignment depth per chromosome (and per region in the case of targeted sequencing)
- DeepVariant SNP/indel gVCF file
- Joint phased DeepTrio SNP/indel VCF file
- Joint phased and annotated DeepTrio SNP/indel VCF file (hg38 only)
- Bed and bigwig base modification frequencies for complete read set and separate haplotypes (uBAMs containing base modifications only)
- Phased tandem repeat VCF file
- Joint phased Sniffles2 and/or un-phased cuteSV SV VCF file
- Joint phased and annotated (VEP + SVscanner) Sniffles2 and/or un-phased and annotated cuteSV SV VCF file (hg38 only)
- SVscanner text diagrams of the repeat elements annotated in each joint SV (hg38 only)
- Puzzleapp joint SNP/indel and SV TSV files and coverage/VAF quality control HTML file (hg38 only)
- Joint phased tandem repeat VCF file
- Somalier extracted files
- Joint relatedness and quality control somalier TSV and HTML files
- ClairS-TO somatic SNV/indel VCF files

> [!NOTE]
> - Running DeepVariant/DeepTrio on ONT data assumes r10 data
> - Running base modification analyses assumes the input data is in uBAM format and base modifications are present in these data

## Haploid Aware Mode

- Enables correct handling of the haploid nature of chrX and chrY for XY samples, along with PAR regions
- Only supported for singletons at the moment

