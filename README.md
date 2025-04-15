
# Pipeface

## Overview

Pipefaceee.

Nextflow pipeline to process long read [ONT](https://nanoporetech.com/) and/or [pacbio](https://www.pacb.com/) HiFi data.

There currently exists tools and workflows that undertake comparable analyses, but pipeface serves as a central workflow to process long read data (both ONT and pacbio HiFi data). Pipeface's future hold's STR, CNV and tandem repeat calling.

<p align="center">
    <img src="./images/pipeface.png">

## Workflow

### Overview

#### Singleton

```mermaid
%%{init: {'theme':'dark'}}%%
flowchart LR

input_data("Input data: <br><br> ONT fastq.gz <br> and/or <br> ONT fastq <br> and/or <br> ONT uBAM <br> and/or <br> pacbio HiFi uBAM")
merging{{"Merge runs (if needed)"}}
alignment{{"bam to fastq conversion (if needed), alignment, sorting"}}
depth{{"Calculate alignment depth"}}
snp_indel_calling{{"SNP/indel variant calling"}}
split_multiallele{{"Split multiallelic variants into biallelic variants"}}
snp_indel_phasing{{"SNP/indel phasing"}}
snp_indel_annotation{{"SNP/indel annotation (hg38 only)"}}
haplotagging{{"Haplotagging bams"}}
calculate_base_mod_freqs{{"Calculate base modificiation frequencies (ONT uBAM's containing base modifications only)"}}
generate_meth_probs{{"Generate site methylation probabilities (pacbio uBAM's containing base modifications only)"}}
sv_calling{{"Structural variant calling"}}
sv_annotation{{"Structural variant annotation (hg38 only)"}}

input_data-.->merging-.->alignment-.->snp_indel_calling-.->split_multiallele-.->snp_indel_phasing-.->haplotagging-.->sv_calling
alignment-.->depth
alignment-.->haplotagging
haplotagging-.->calculate_base_mod_freqs
haplotagging-.->generate_meth_probs
snp_indel_phasing-.->snp_indel_annotation
sv_calling-.->sv_annotation
```

#### Cohort

```mermaid
%%{init: {'theme':'dark'}}%%
flowchart LR

input_data("Input data: <br><br> ONT fastq.gz <br> and/or <br> ONT fastq <br> and/or <br> ONT uBAM <br> and/or <br> pacbio HiFi uBAM")
merging{{"Merge runs (if needed)"}}
alignment{{"bam to fastq conversion (if needed), alignment, sorting"}}
depth{{"Calculate alignment depth"}}
snp_indel_calling{{"SNP/indel variant calling"}}
split_multiallele{{"Split multiallelic variants into biallelic variants"}}
snp_indel_phasing{{"SNP/indel phasing"}}
joint_snp_indel_calling{{"Joint SNP/indel variant calling"}}
gvcf_merging{{"gVCF merging"}}
joint_split_multiallele{{"Split multiallelic variants into biallelic variants"}}
joint_snp_indel_phasing{{"Joint SNP/indel phasing"}}
joint_snp_indel_annotation{{"Joint SNP/indel annotation (hg38 only)"}}
haplotagging{{"Haplotagging bams"}}
calculate_base_mod_freqs{{"Calculate base modificiation frequencies (ONT uBAM's containing base modifications only)"}}
generate_meth_probs{{"Generate site methylation probabilities (pacbio uBAM's containing base modifications only)"}}
sv_calling{{"Structural variant calling"}}
sv_vcf_merging{{"Structural variant VCF merging"}}
joint_sv_annotation{{"Joint structural variant annotation (hg38 only)"}}

input_data-.->merging-.->alignment-.->snp_indel_calling-.->split_multiallele-.->snp_indel_phasing-.->haplotagging-.->sv_calling
alignment-.->depth
alignment-.->haplotagging
haplotagging-.->calculate_base_mod_freqs
haplotagging-.->generate_meth_probs
snp_indel_phasing-.->joint_snp_indel_calling-.->gvcf_merging-.->joint_split_multiallele-.->joint_snp_indel_phasing-.->joint_snp_indel_annotation
sv_calling-.->sv_vcf_merging-.->joint_sv_annotation

```

### Detailed

#### Singleton

```mermaid
%%{init: {'theme':'dark'}}%%
flowchart LR

ont_data_f1("Sample 1 <br><br> Input data: <br><br> ONT fastq.gz")
ont_data_f2("Sample 1 <br><br> Input data: <br><br> ONT fastq.gz")
pacbio_data_f3("Sample 2 <br><br> Input data: <br><br> Pacbio HiFi uBAM")
pacbio_data_f4("Sample 2 <br><br> Input data: <br><br> Pacbio HiFi uBAM")
ont_data_f5("Sample 3 <br><br> Input data: <br><br> ONT fastq")
ont_data_f6("Sample 4 <br><br> Input data: <br><br> ONT uBAM")

merging_m1{{"Description: merge runs <br><br> Main tools: GNU coreutils <br><br> Commands: cat"}}
merging_m2{{"Description: merge runs <br><br> Main tools: Samtools <br><br> Commands: samtools merge"}}

alignment_s1{{"Description: alignment, sorting <br><br> Main tools: Minimap2 and Samtools <br><br> Commands: minimap2 and samtools sort"}}
alignment_s2{{"Description: bam to fastq conversion, alignment, sorting <br><br> Main tools: Minimap2 and Samtools <br><br> Commands: minimap2 and samtools sort"}}
alignment_s3{{"Description: alignment, sorting <br><br> Main tools: Minimap2 and Samtools <br><br> Commands: minimap2 and samtools sort"}}
alignment_s4{{"Description: bam to fastq conversion, alignment, sorting <br><br> Main tools: Minimap2 and Samtools <br><br> Commands: minimap2 and samtools sort"}}

depth_s1{{"Description: calculate alignment depth <br><br> Main tools: mosdepth <br><br> Commands: mosdepth depth"}}
depth_s2{{"Description: calculate alignment depth <br><br> Main tools: mosdepth <br><br> Commands: mosdepth depth"}}
depth_s3{{"Description: calculate alignment depth <br><br> Main tools: mosdepth <br><br> Commands: mosdepth depth"}}
depth_s4{{"Description: calculate alignment depth <br><br> Main tools: mosdepth <br><br> Commands: mosdepth depth"}}

snp_indel_calling_s1{{"Description: SNP/indel variant calling <br><br> Main tools: Clair3 or DeepVariant <br><br> Commands: run_clair3.sh or run_deepvariant"}}
snp_indel_calling_s2{{"Description: SNP/indel variant calling <br><br> Main tools: Clair3 or DeepVariant <br><br> Commands: run_clair3.sh or run_deepvariant"}}
snp_indel_calling_s3{{"Description: SNP/indel variant calling <br><br> Main tools: Clair3 or DeepVariant <br><br> Commands: run_clair3.sh or run_deepvariant"}}
snp_indel_calling_s4{{"Description: SNP/indel variant calling <br><br> Main tools: Clair3 or DeepVariant <br><br> Commands: run_clair3.sh or run_deepvariant"}}

split_multiallele_s1{{"Description: Split multiallelic variants into biallelic variants <br><br> Main tools: BCFtools <br><br> Commands: bcftools norm"}}
split_multiallele_s2{{"Description: Split multiallelic variants into biallelic variants <br><br> Main tools: BCFtools <br><br> Commands: bcftools norm"}}
split_multiallele_s3{{"Description: Split multiallelic variants into biallelic variants <br><br> Main tools: BCFtools <br><br> Commands: bcftools norm"}}
split_multiallele_s4{{"Description: Split multiallelic variants into biallelic variants <br><br> Main tools: BCFtools <br><br> Commands: bcftools norm"}}

snp_indel_phasing_s1{{"Description: SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase and whatshap stats"}}
snp_indel_phasing_s2{{"Description: SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase and whatshap stats"}}
snp_indel_phasing_s3{{"Description: SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase and whatshap stats"}}
snp_indel_phasing_s4{{"Description: SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase and whatshap stats"}}

snp_indel_annotation_s1{{"Description: SNP/indel annotation (hg38 only)" <br><br> Main tools: ensembl-vep <br><br> Commands: vep}}
snp_indel_annotation_s2{{"Description: SNP/indel annotation (hg38 only)" <br><br> Main tools: ensembl-vep <br><br> Commands: vep}}
snp_indel_annotation_s3{{"Description: SNP/indel annotation (hg38 only)" <br><br> Main tools: ensembl-vep <br><br> Commands: vep}}
snp_indel_annotation_s4{{"Description: SNP/indel annotation (hg38 only)" <br><br> Main tools: ensembl-vep <br><br> Commands: vep}}

haplotagging_s1{{"Description: haplotagging bams <br><br> Main tools: WhatsHap <br><br> Commands: whatshap haplotag"}}
haplotagging_s2{{"Description: haplotagging bams <br><br> Main tools: WhatsHap <br><br> Commands: whatshap haplotag"}}
haplotagging_s3{{"Description: haplotagging bams <br><br> Main tools: WhatsHap <br><br> Commands: whatshap haplotag"}}
haplotagging_s4{{"Description: haplotagging bams <br><br> Main tools: WhatsHap <br><br> Commands: whatshap haplotag"}}

generate_meth_probs_s2{{"Description: generate site methylation probabilities <br><br> Main tools: pb-CpG-tools <br><br> Commands: aligned_bam_to_cpg_scores"}}
calculate_base_mod_freqs_s4{{"Description: calculate base modificiation frequencies <br><br> Main tools: minimod <br><br> Commands: minimod mod-freq"}}

sv_calling_s1{{"Description: structural variant calling <br><br> Main tools: Sniffles2 and/or cuteSV <br><br> Commands: sniffles and/or cuteSV"}}
sv_calling_s2{{"Description: structural variant calling <br><br> Main tools: Sniffles2 and/or cuteSV <br><br> Commands: sniffles and/or cuteSV"}}
sv_calling_s3{{"Description: structural variant calling <br><br> Main tools: Sniffles2 and/or cuteSV <br><br> Commands: sniffles and/or cuteSV"}}
sv_calling_s4{{"Description: structural variant calling <br><br> Main tools: Sniffles2 and/or cuteSV <br><br> Commands: sniffles and/or cuteSV"}}

sv_annotation_s1{{"Description: structural variant annotation (hg38 only)" <br><br> Main tools: ensembl-vep <br><br> Commands: vep}}
sv_annotation_s2{{"Description: structural variant annotation (hg38 only)" <br><br> Main tools: ensembl-vep <br><br> Commands: vep}}
sv_annotation_s3{{"Description: structural variant annotation (hg38 only)" <br><br> Main tools: ensembl-vep <br><br> Commands: vep}}
sv_annotation_s4{{"Description: structural variant annotation (hg38 only)" <br><br> Main tools: ensembl-vep <br><br> Commands: vep}}

ont_data_f1-.->merging_m1-.->alignment_s1-.->snp_indel_calling_s1-.->split_multiallele_s1-.->snp_indel_phasing_s1-.->haplotagging_s1-.->sv_calling_s1
ont_data_f2-.->merging_m1
pacbio_data_f3-.->merging_m2-.->alignment_s2-.->snp_indel_calling_s2-.->split_multiallele_s2-.->snp_indel_phasing_s2-.->haplotagging_s2-.->sv_calling_s2
pacbio_data_f4-.->merging_m2
ont_data_f5-.->alignment_s3-.->snp_indel_calling_s3-.->split_multiallele_s3-.->snp_indel_phasing_s3-.->haplotagging_s3-.->sv_calling_s3
ont_data_f6-.->alignment_s4-.->snp_indel_calling_s4-.->split_multiallele_s4-.->snp_indel_phasing_s4-.->haplotagging_s4-.->sv_calling_s4

alignment_s1-.->depth_s1
alignment_s2-.->depth_s2
alignment_s3-.->depth_s3
alignment_s4-.->depth_s4

alignment_s1-.->haplotagging_s1
alignment_s2-.->haplotagging_s2
alignment_s3-.->haplotagging_s3
alignment_s4-.->haplotagging_s4

haplotagging_s2-.->generate_meth_probs_s2
haplotagging_s4-.->calculate_base_mod_freqs_s4

snp_indel_phasing_s1-.->snp_indel_annotation_s1
snp_indel_phasing_s2-.->snp_indel_annotation_s2
snp_indel_phasing_s3-.->snp_indel_annotation_s3
snp_indel_phasing_s4-.->snp_indel_annotation_s4

sv_calling_s1-.->sv_annotation_s1
sv_calling_s2-.->sv_annotation_s2
sv_calling_s3-.->sv_annotation_s3
sv_calling_s4-.->sv_annotation_s4

```

#### Cohort

```mermaid
%%{init: {'theme':'dark'}}%%
flowchart LR

ont_data_f1("Sample 1 <br><br> Input data: <br><br> ONT uBAM")
ont_data_f2("Sample 1 <br><br> Input data: <br><br> ONT uBAM")
ont_data_f3("Sample 2 <br><br> Input data: <br><br> ONT uBAM")
ont_data_f4("Sample 3 <br><br> Input data: <br><br> ONT uBAM")
ont_data_f5("Sample 4 <br><br> Input data: <br><br> ONT uBAM")
ont_data_f6("Sample 5 <br><br> Input data: <br><br> ONT uBAM")
ont_data_f7("Sample 6 <br><br> Input data: <br><br> ONT uBAM")

merging_m1{{"Description: merge runs <br><br> Main tools: Samtools <br><br> Commands: samtools merge"}}

alignment_s1{{"Description: bam to fastq conversion, alignment, sorting <br><br> Main tools: Minimap2 and Samtools <br><br> Commands: minimap2 and samtools sort"}}
alignment_s2{{"Description: bam to fastq conversion, alignment, sorting <br><br> Main tools: Minimap2 and Samtools <br><br> Commands: minimap2 and samtools sort"}}
alignment_s3{{"Description: bam to fastq conversion, alignment, sorting <br><br> Main tools: Minimap2 and Samtools <br><br> Commands: minimap2 and samtools sort"}}
alignment_s4{{"Description: bam to fastq conversion, alignment, sorting <br><br> Main tools: Minimap2 and Samtools <br><br> Commands: minimap2 and samtools sort"}}
alignment_s5{{"Description: bam to fastq conversion, alignment, sorting <br><br> Main tools: Minimap2 and Samtools <br><br> Commands: minimap2 and samtools sort"}}
alignment_s6{{"Description: bam to fastq conversion, alignment, sorting <br><br> Main tools: Minimap2 and Samtools <br><br> Commands: minimap2 and samtools sort"}}

depth_s1{{"Description: calculate alignment depth <br><br> Main tools: mosdepth <br><br> Commands: mosdepth depth"}}
depth_s2{{"Description: calculate alignment depth <br><br> Main tools: mosdepth <br><br> Commands: mosdepth depth"}}
depth_s3{{"Description: calculate alignment depth <br><br> Main tools: mosdepth <br><br> Commands: mosdepth depth"}}
depth_s4{{"Description: calculate alignment depth <br><br> Main tools: mosdepth <br><br> Commands: mosdepth depth"}}
depth_s5{{"Description: calculate alignment depth <br><br> Main tools: mosdepth <br><br> Commands: mosdepth depth"}}
depth_s6{{"Description: calculate alignment depth <br><br> Main tools: mosdepth <br><br> Commands: mosdepth depth"}}

snp_indel_calling_s1{{"Description: SNP/indel variant calling <br><br> Main tools: DeepVariant <br><br> Commands: run_deepvariant"}}
snp_indel_calling_s2{{"Description: SNP/indel variant calling <br><br> Main tools: DeepVariant <br><br> Commands: run_deepvariant"}}
snp_indel_calling_s3{{"Description: SNP/indel variant calling <br><br> Main tools: DeepVariant <br><br> Commands: run_deepvariant"}}
snp_indel_calling_s4{{"Description: SNP/indel variant calling <br><br> Main tools: DeepVariant <br><br> Commands: run_deepvariant"}}
snp_indel_calling_s5{{"Description: SNP/indel variant calling <br><br> Main tools: DeepVariant <br><br> Commands: run_deepvariant"}}
snp_indel_calling_s6{{"Description: SNP/indel variant calling <br><br> Main tools: DeepVariant <br><br> Commands: run_deepvariant"}}

split_multiallele_s1{{"Description: Split multiallelic variants into biallelic variants <br><br> Main tools: BCFtools <br><br> Commands: bcftools norm"}}
split_multiallele_s2{{"Description: Split multiallelic variants into biallelic variants <br><br> Main tools: BCFtools <br><br> Commands: bcftools norm"}}
split_multiallele_s3{{"Description: Split multiallelic variants into biallelic variants <br><br> Main tools: BCFtools <br><br> Commands: bcftools norm"}}
split_multiallele_s4{{"Description: Split multiallelic variants into biallelic variants <br><br> Main tools: BCFtools <br><br> Commands: bcftools norm"}}
split_multiallele_s5{{"Description: Split multiallelic variants into biallelic variants <br><br> Main tools: BCFtools <br><br> Commands: bcftools norm"}}
split_multiallele_s6{{"Description: Split multiallelic variants into biallelic variants <br><br> Main tools: BCFtools <br><br> Commands: bcftools norm"}}

snp_indel_phasing_s1{{"Description: SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase"}}
snp_indel_phasing_s2{{"Description: SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase"}}
snp_indel_phasing_s3{{"Description: SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase"}}
snp_indel_phasing_s4{{"Description: SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase"}}
snp_indel_phasing_s5{{"Description: SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase"}}
snp_indel_phasing_s6{{"Description: SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase"}}

haplotagging_s1{{"Description: haplotagging bams <br><br> Main tools: WhatsHap <br><br> Commands: whatshap haplotag"}}
haplotagging_s2{{"Description: haplotagging bams <br><br> Main tools: WhatsHap <br><br> Commands: whatshap haplotag"}}
haplotagging_s3{{"Description: haplotagging bams <br><br> Main tools: WhatsHap <br><br> Commands: whatshap haplotag"}}
haplotagging_s4{{"Description: haplotagging bams <br><br> Main tools: WhatsHap <br><br> Commands: whatshap haplotag"}}
haplotagging_s5{{"Description: haplotagging bams <br><br> Main tools: WhatsHap <br><br> Commands: whatshap haplotag"}}
haplotagging_s6{{"Description: haplotagging bams <br><br> Main tools: WhatsHap <br><br> Commands: whatshap haplotag"}}

calculate_base_mod_freqs_s1{{"Description: calculate base modificiation frequencies <br><br> Main tools: minimod <br><br> Commands: minimod mod-freq"}}
calculate_base_mod_freqs_s2{{"Description: calculate base modificiation frequencies <br><br> Main tools: minimod <br><br> Commands: minimod mod-freq"}}
calculate_base_mod_freqs_s3{{"Description: calculate base modificiation frequencies <br><br> Main tools: minimod <br><br> Commands: minimod mod-freq"}}
calculate_base_mod_freqs_s4{{"Description: calculate base modificiation frequencies <br><br> Main tools: minimod <br><br> Commands: minimod mod-freq"}}
calculate_base_mod_freqs_s5{{"Description: calculate base modificiation frequencies <br><br> Main tools: minimod <br><br> Commands: minimod mod-freq"}}
calculate_base_mod_freqs_s6{{"Description: calculate base modificiation frequencies <br><br> Main tools: minimod <br><br> Commands: minimod mod-freq"}}

joint_snp_indel_calling_f1{{"Description: Joint SNP/indel variant calling <br><br> Main tools: DeepTrio <br><br> Commands: run_deeptrio"}}
joint_snp_indel_calling_f2{{"Description: Joint SNP/indel variant calling <br><br> Main tools: DeepTrio <br><br> Commands: run_deeptrio"}}

gvcf_merging_f1{{"Description: gVCF merging <br><br> Main tools: GLnexus <br><br> Commands: glnexus_cli"}}
gvcf_merging_f2{{"Description: gVCF merging <br><br> Main tools: GLnexus <br><br> Commands: glnexus_cli"}}

joint_split_multiallele_f1{{"Description: Split multiallelic variants into biallelic variants <br><br> Main tools: BCFtools <br><br> Commands: bcftools norm"}}
joint_split_multiallele_f2{{"Description: Split multiallelic variants into biallelic variants <br><br> Main tools: BCFtools <br><br> Commands: bcftools norm"}}

joint_snp_indel_phasing_f1{{"Description: Joint SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase and and whatshap stats"}}
joint_snp_indel_phasing_f2{{"Description: Joint SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase and and whatshap stats"}}

joint_snp_indel_annotation_f1{{"Description: Joint SNP/indel annotation (hg38 only)" <br><br> Main tools: ensembl-vep <br><br> Commands: vep}}
joint_snp_indel_annotation_f2{{"Description: Joint SNP/indel annotation (hg38 only)" <br><br> Main tools: ensembl-vep <br><br> Commands: vep}}

sv_calling_s1{{"Description: structural variant calling <br><br> Main tools: Sniffles2 and/or cuteSV <br><br> Commands: sniffles and/or cuteSV"}}
sv_calling_s2{{"Description: structural variant calling <br><br> Main tools: Sniffles2 and/or cuteSV <br><br> Commands: sniffles and/or cuteSV"}}
sv_calling_s3{{"Description: structural variant calling <br><br> Main tools: Sniffles2 and/or cuteSV <br><br> Commands: sniffles and/or cuteSV"}}
sv_calling_s4{{"Description: structural variant calling <br><br> Main tools: Sniffles2 and/or cuteSV <br><br> Commands: sniffles and/or cuteSV"}}
sv_calling_s5{{"Description: structural variant calling <br><br> Main tools: Sniffles2 and/or cuteSV <br><br> Commands: sniffles and/or cuteSV"}}
sv_calling_s6{{"Description: structural variant calling <br><br> Main tools: Sniffles2 and/or cuteSV <br><br> Commands: sniffles and/or cuteSV"}}

sv_vcf_merging_f1{{"Description: Structural variant VCF merging <br><br> Main tools: Jasmine <br><br> Commands: jasmine"}}
sv_vcf_merging_f2{{"Description: Structural variant VCF merging <br><br> Main tools: Jasmine <br><br> Commands: jasmine"}}

joint_sv_annotation_s1{{"Description: Joint structural variant annotation (hg38 only)" <br><br> Main tools: ensembl-vep <br><br> Commands: vep}}
joint_sv_annotation_s2{{"Description: Joint structural variant annotation (hg38 only)" <br><br> Main tools: ensembl-vep <br><br> Commands: vep}}

ont_data_f1-.->merging_m1-.->alignment_s1-.->snp_indel_calling_s1-.->split_multiallele_s1-.->snp_indel_phasing_s1-.->haplotagging_s1-.->sv_calling_s1
ont_data_f2-.->merging_m1
ont_data_f3-.->alignment_s2-.->snp_indel_calling_s2-.->split_multiallele_s2-.->snp_indel_phasing_s2-.->haplotagging_s2-.->sv_calling_s2
ont_data_f4-.->alignment_s3-.->snp_indel_calling_s3-.->split_multiallele_s3-.->snp_indel_phasing_s3-.->haplotagging_s3-.->sv_calling_s3
ont_data_f5-.->alignment_s4-.->snp_indel_calling_s4-.->split_multiallele_s4-.->snp_indel_phasing_s4-.->haplotagging_s4-.->sv_calling_s4
ont_data_f6-.->alignment_s5-.->snp_indel_calling_s5-.->split_multiallele_s5-.->snp_indel_phasing_s5-.->haplotagging_s5-.->sv_calling_s5
ont_data_f7-.->alignment_s6-.->snp_indel_calling_s6-.->split_multiallele_s6-.->snp_indel_phasing_s6-.->haplotagging_s6-.->sv_calling_s6

alignment_s1-.->depth_s1
alignment_s2-.->depth_s2
alignment_s3-.->depth_s3
alignment_s4-.->depth_s4
alignment_s5-.->depth_s5
alignment_s6-.->depth_s6

alignment_s1-.->haplotagging_s1
alignment_s2-.->haplotagging_s2
alignment_s3-.->haplotagging_s3
alignment_s4-.->haplotagging_s4
alignment_s5-.->haplotagging_s5
alignment_s6-.->haplotagging_s6

haplotagging_s1-.->calculate_base_mod_freqs_s1
haplotagging_s2-.->calculate_base_mod_freqs_s2
haplotagging_s3-.->calculate_base_mod_freqs_s3
haplotagging_s4-.->calculate_base_mod_freqs_s4
haplotagging_s5-.->calculate_base_mod_freqs_s5
haplotagging_s6-.->calculate_base_mod_freqs_s6

haplotagging_s1-.->joint_snp_indel_calling_f1
haplotagging_s2-.->joint_snp_indel_calling_f1
haplotagging_s3-.->joint_snp_indel_calling_f1
haplotagging_s4-.->joint_snp_indel_calling_f2
haplotagging_s5-.->joint_snp_indel_calling_f2
haplotagging_s6-.->joint_snp_indel_calling_f2

joint_snp_indel_calling_f1-.->gvcf_merging_f1-.->joint_split_multiallele_f1-.->joint_snp_indel_phasing_f1-.->joint_snp_indel_annotation_f1
joint_snp_indel_calling_f2-.->gvcf_merging_f2-.->joint_split_multiallele_f2-.->joint_snp_indel_phasing_f2-.->joint_snp_indel_annotation_f2

sv_calling_s1-.->sv_vcf_merging_f1
sv_calling_s2-.->sv_vcf_merging_f1
sv_calling_s3-.->sv_vcf_merging_f1
sv_calling_s4-.->sv_vcf_merging_f2
sv_calling_s5-.->sv_vcf_merging_f2
sv_calling_s6-.->sv_vcf_merging_f2

sv_vcf_merging_f1-.->joint_sv_annotation_s1
sv_vcf_merging_f2-.->joint_sv_annotation_s2

```


## Main analyses

- ONT and/or pacbio HiFi data
- Individuals or cohorts
- WGS and/or targeted
- hg38 or hs1 reference genome

## Main tools

- [Minimap2](https://github.com/lh3/minimap2)
- [Clair3](https://github.com/HKU-BAL/Clair3) or [DeepVariant](https://github.com/google/deepvariant)/[DeepTrio](https://github.com/google/deepvariant/blob/r1.8/docs/deeptrio-details.md)
- [WhatsHap](https://github.com/whatshap/whatshap)
- [GLnexus](https://github.com/dnanexus-rnd/GLnexus)
- [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) and/or [cuteSV](https://github.com/tjiangHIT/cuteSV)
- [Jasmine (customised)](https://github.com/bioinfomethods/Jasmine)
- [Samtools](https://github.com/samtools/samtools)
- [mosdepth](https://github.com/brentp/mosdepth)
- [minimod](https://github.com/warp9seq/minimod?tab=readme-ov-file)
- [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools)
- [ensembl-vep](https://github.com/Ensembl/ensembl-vep)

## Main input files

### Required

- ONT/pacbio HiFi FASTQ (gzipped or uncompressed) or unaligned BAM
- Indexed reference genome
- Clair3 models (if running Clair3)

### Optional

- Regions of interest BED file
- Tandem repeat BED file

## Main output files

### Singleton

- Aligned, sorted and haplotagged bam
- Alignment depth per chromosome (and per region in the case of targeted sequencing)
- Phased Clair3 or DeepVariant SNP/indel VCF file
- Phased and annotated Clair3 or DeepVariant SNP/indel VCF file (hg38 only)
- Clair3 or DeepVariant SNP/indel gVCF file
- Bed and bigwig base modification frequencies for complete read set and separate haplotypes (ONT uBAM's containing base modifications only)
- Bed and bigwig site methylation probabilities for complete read set and separate haplotypes (pacbio uBAM's containing base modifications only)
- Phased Sniffles2 and/or un-phased cuteSV SV VCF file
- Phased and annotated Sniffles2 and/or un-phased and annotated cuteSV SV VCF file (hg38 only)

### Cohort

- Aligned, sorted and haplotagged bam
- Alignment depth per chromosome (and per region in the case of targeted sequencing)
- Clair3 or DeepVariant SNP/indel gVCF file
- Joint phased DeepTrio SNP/indel VCF file
- Joint phased and annotated DeepTrio SNP/indel VCF file (hg38 only)
- Bed and bigwig base modification frequencies for complete read set and separate haplotypes (ONT uBAM's containing base modifications only)
- Bed and bigwig site methylation probabilities for complete read set and separate haplotypes (pacbio uBAM's containing base modifications only)
- Joint phased Sniffles2 and/or un-phased cuteSV SV VCF file
- Joint phased and annotated Sniffles2 and/or un-phased and annotated cuteSV SV VCF file (hg38 only)

> **_Note:_** Running DeepVariant/DeepTrio on ONT data assumes r10 data

> **_Note:_** Running base modification analyses assume the input data is un uBAM format and base modifications are present in these data

## Haploid Aware Mode 

- Enables correct handling of the haploid nature of chrX and chrY for XY samples, along with PAR regions
- Only supported for singletons at the moment. 

## Assumptions

- Running pipeline on Australia's [National Computational Infrastructure (NCI)](https://nci.org.au/)
- Access to if89 project (to access software installs used by pipeface)
- Access to xy86 project (to access variant databases used by pipeface, only required if running variant annotation)

*[See the list of software and their versions used by this version of pipeface](./docs/software_versions.txt) as well as the [list of variant databases and their versions](./docs/database_versions.txt) if variant annotation is carried out (assuming the default [nextflow_pipeface.config](./config/nextflow_pipeface.config) file is used).*

## Run it!

See a walkthrough for how to [run pipeface on NCI](./docs/run_on_nci.md).

## Credit

This is a highly collaborative project, with many contributions from the [Genomic Technologies Lab](https://www.garvan.org.au/research/labs-groups/genomic-technologies-lab). Notably, Dr Andre Reis and Dr Ira Deveson are closely involved in the development of this pipeline. Optimisations involving DeepVariant and DeepTrio have been contributed by Dr Kisaru Liyanage and Dr Matthew Downton from the [National Computational Infrastructure](https://nci.org.au), with support from Australian BioCommons as part of the Workflow Commons project. Haploid-aware mode has been contributed by Dr Hardip Patel & Kirat Alreja from the [National Centre for Indigenous Genomics](https://ncig.anu.edu.au). The installation and hosting of software used in this pipeline has and continues to be supported by the [Australian BioCommons Tools and Workflows project (if89)](https://australianbiocommons.github.io/ables/if89/).

