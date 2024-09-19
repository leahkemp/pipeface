# Pipeface

## Overview

Pipefaceee.

Nextflow pipeline to align, variant call (SNP's, indels's, SV's) and phase long read [ONT](https://nanoporetech.com/) and/or [pacbio](https://www.pacb.com/) HiFi data.

There currently exists tools and workflows which align, variant call and phase ONT/pacbio HiFi data, but pipeface serves as a central workflow to process long read data. Pipeface's future hold's STR, CNV and tandem repeat calling, as well as the analysis of cohorts.

<p align="center">
    <img src="./images/pipeface.png">

## Workflow

### Overview

```mermaid
%%{init: {'theme':'dark'}}%%
flowchart LR

input_data("Input data: <br><br> ONT fastq.gz <br> and/or <br> ONT fastq <br> and/or <br> ONT uBAM <br> and/or <br> pacbio HiFi uBAM")
merging{{"Merge runs (if needed)"}}
alignment{{"bam to fastq conversion (if needed), alignment, sorting"}}
depth{{"Calculate alignment depth"}}
snp_indel_calling{{"SNP/indel variant calling"}}
snp_indel_phasing{{"SNP/indel phasing"}}
haplotagging{{"Haplotagging bams"}}
sv_calling{{"Structural variant calling"}}

input_data-.->merging-.->alignment-.->snp_indel_calling-.->snp_indel_phasing-.->haplotagging-.->sv_calling
alignment-.->depth
alignment-.->haplotagging

```

### Detailed

```mermaid
%%{init: {'theme':'dark'}}%%
flowchart LR

ont_data_f1("Sample 1 <br><br> Input data: <br> ONT fastq.gz")
ont_data_f2("Sample 1 <br><br> Input data: <br> ONT fastq.gz")
pacbio_data_f3("Sample 2 <br><br> Input data: <br> Pacbio HiFi uBAM")
pacbio_data_f4("Sample 2 <br><br> Input data: <br> Pacbio HiFi uBAM")
ont_data_f5("Sample 3 <br><br> Input data: <br> ONT fastq")
ont_data_f6("Sample 4 <br><br> Input data: <br> ONT uBAM")

merging_m1{{"Description: merge runs <br><br> Main tools: GNU coreutils <br><br> Commands: cat"}}
merging_m2{{"Description: merge runs <br><br> Main tools: Samtools <br><br> Commands: samtools merge"}}

alignment_s1{{"Description: alignment, sorting <br><br> Main tools: Minimap2 and Samtools <br><br> Commands: minimap2 and samtools sort"}}
alignment_s2{{"Description: alignment, sorting <br><br> Main tools: Minimap2 and Samtools <br><br> Commands: minimap2 and samtools sort"}}
alignment_s3{{"Description: bam to fastq conversion, alignment, sorting <br><br> Main tools: Minimap2 and Samtools <br><br> Commands: minimap2 and samtools sort"}}
alignment_s4{{"Description: bam to fastq conversion, alignment, sorting <br><br> Main tools: Minimap2 and Samtools <br><br> Commands: minimap2 and samtools sort"}}

depth_s1{{"Description: calculate alignment depth <br><br> Main tools: Samtools <br><br> Commands: samtools depth"}}
depth_s2{{"Description: calculate alignment depth <br><br> Main tools: Samtools <br><br> Commands: samtools depth"}}
depth_s3{{"Description: calculate alignment depth <br><br> Main tools: Samtools <br><br> Commands: samtools depth"}}
depth_s4{{"Description: calculate alignment depth <br><br> Main tools: Samtools <br><br> Commands: samtools depth"}}

snp_indel_calling_s1{{"Description: SNP/indel variant calling <br><br> Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) <br><br> Commands: run_clair3.sh or pbrun deepvariant"}}
snp_indel_calling_s2{{"Description: SNP/indel variant calling <br><br> Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) <br><br> Commands: run_clair3.sh or pbrun deepvariant"}}
snp_indel_calling_s3{{"Description: SNP/indel variant calling <br><br> Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) <br><br> Commands: run_clair3.sh or pbrun deepvariant"}}
snp_indel_calling_s4{{"Description: SNP/indel variant calling <br><br> Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) <br><br> Commands: run_clair3.sh or pbrun deepvariant"}}

snp_indel_phasing_s1{{"Description: SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase"}}
snp_indel_phasing_s2{{"Description: SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase"}}
snp_indel_phasing_s3{{"Description: SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase"}}
snp_indel_phasing_s4{{"Description: SNP/indel phasing <br><br> Main tools: WhatsHap <br><br> Commands: whatshap phase"}}

haplotagging_s1{{"Description: haplotagging bams <br><br> Main tools: WhatsHap <br><br> Commands: whatshap haplotag"}}
haplotagging_s2{{"Description: haplotagging bams <br><br> Main tools: WhatsHap <br><br> Commands: whatshap haplotag"}}
haplotagging_s3{{"Description: haplotagging bams <br><br> Main tools: WhatsHap <br><br> Commands: whatshap haplotag"}}
haplotagging_s4{{"Description: haplotagging bams <br><br> Main tools: WhatsHap <br><br> Commands: whatshap haplotag"}}

sv_calling_s1{{"Description: structural variant calling <br><br> Main tools: Sniffles2 and/or cuteSV <br><br> Commands: sniffles and/or cuteSV"}}
sv_calling_s2{{"Description: structural variant calling <br><br> Main tools: Sniffles2 and/or cuteSV <br><br> Commands: sniffles and/or cuteSV"}}
sv_calling_s3{{"Description: structural variant calling <br><br> Main tools: Sniffles2 and/or cuteSV <br><br> Commands: sniffles and/or cuteSV"}}
sv_calling_s4{{"Description: structural variant calling <br><br> Main tools: Sniffles2 and/or cuteSV <br><br> Commands: sniffles and/or cuteSV"}}

ont_data_f1-.->merging_m1-.->alignment_s1-.->snp_indel_calling_s1-.->snp_indel_phasing_s1-.->haplotagging_s1-.->sv_calling_s1
ont_data_f2-.->merging_m1
ont_data_f5-.->alignment_s2-.->snp_indel_calling_s2-.->snp_indel_phasing_s2-.->haplotagging_s2-.->sv_calling_s2
ont_data_f6-.->alignment_s3-.->snp_indel_calling_s3-.->snp_indel_phasing_s3-.->haplotagging_s3-.->sv_calling_s3

pacbio_data_f3-.->merging_m2-.->alignment_s4-.->snp_indel_calling_s4-.->snp_indel_phasing_s4-.->haplotagging_s4-.->sv_calling_s4
pacbio_data_f4-.->merging_m2

alignment_s1-.->depth_s1
alignment_s2-.->depth_s2
alignment_s3-.->depth_s3
alignment_s4-.->depth_s4

alignment_s1-.->haplotagging_s1
alignment_s2-.->haplotagging_s2
alignment_s3-.->haplotagging_s3
alignment_s4-.->haplotagging_s4

```

## Main analyses

- ONT and pacbio HiFi data
- WGS and targeted
- hg38 or hs1 reference genome

## Main tools

- [Minimap2](https://github.com/lh3/minimap2)
- [Clair3](https://github.com/HKU-BAL/Clair3) or [DeepVariant](https://github.com/google/deepvariant) (wrapped in [NVIDIA Parabricks](https://docs.nvidia.com/clara/parabricks/latest/))
- [WhatsHap](https://github.com/whatshap/whatshap)
- [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) and/or [cuteSV](https://github.com/tjiangHIT/cuteSV)
- [Samtools](https://github.com/samtools/samtools)

## Main input files

### Required

- ONT/pacbio HiFi FASTQ (gzipped or uncompressed) or unaligned BAM
- Indexed reference genome
- Clair3 models (if running Clair3)

### Optional

- Regions of interest BED file
- Tandem repeat BED file

## Main output files

- Aligned, sorted and haplotagged bam
- Clair3 or DeepVariant phased SNP/indel VCF file
- Clair3 or DeepVariant SNP/indel gVCF file
- Phased Sniffles2 and/or un-phased cuteSV SV VCF file

## Assumptions

- Running pipeline on Australia's [National Computational Infrastructure (NCI)](https://nci.org.au/)
- Access to if89 project on [National Computational Infrastructure (NCI)](https://nci.org.au/)
- Access to pipeline dependencies:
    - [Nextflow and it's java dependency](https://nf-co.re/docs/usage/installation). Validated to run on:
        - Nextflow 24.04.1
        - Java 17.0.2

## Run it!

See a walkthrough for how to [run pipeface on NCI](./docs/run_on_nci.md).

## Credit

This is a highly collaborative project, with many contributions from the [genomic technologies lab](https://www.garvan.org.au/research/labs-groups/genomic-technologies-lab). Notably, Dr Andre Reis and Dr Ira Deveson are closely involved in the development of this pipeline. The installation and hosting of software used in this pipeline has been supported by [project if89 on NCI](https://australianbiocommons.github.io/ables/if89/).
