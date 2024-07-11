# Pipeface

## Overview

Pipefaceee.

Nextflow pipeline to align, variant call (SNP's, indels's, SV's) and phase long read [ONT](https://nanoporetech.com/) and/or [pacbio](https://www.pacb.com/) HiFi data.

<p align="center">
    <img src="./images/pipeface.png">

## Workflow

### Overview

```mermaid
%%{init: {'theme':'dark'}}%%
flowchart LR

input_data("Input data: \n\n - ONT fastq.gz \n - ONT fastq \n - ONT uBAM \n - pacbio HiFi uBAM")
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

ont_data_f1("Input data: \n\n ONT fastq.gz \n\n (sample 1)")
ont_data_f2("Input data: \n\n ONT fastq.gz \n\n (sample 1)")
pacbio_data_f3("Pacbio HiFi uBAM \n (sample 2)")
pacbio_data_f4("Pacbio HiFi uBAM \n (sample 2)")
ont_data_f5("Input data: \n\n ONT fastq \n\n (sample 3)")
ont_data_f6("Input data: \n\n ONT uBAM \n\n (sample 4)")

merging_m1{{"Description: Merge runs \n\n Main tools: Samtools or GNU coreutils \n\n Commands: cat"}}
merging_m2{{"Description: Merge runs \n\n Main tools: Samtools or GNU coreutils \n\n Commands: samtools merge"}}

alignment_s1{{"Description: alignment, sorting \n\n Main tools: Minimap2 and Samtools \n\n Commands: minimap2 and samtools sort"}}
alignment_s2{{"Description: alignment, sorting \n\n Main tools: Minimap2 and Samtools \n\n Commands: minimap2 and samtools sort"}}
alignment_s3{{"Description: bam to fastq conversion, alignment, sorting \n\n Main tools: Minimap2 and Samtools \n\n Commands: minimap2 and samtools sort"}}
alignment_s4{{"Description: bam to fastq conversion, alignment, sorting \n\n Main tools: Minimap2 and Samtools \n\n Commands: minimap2 and samtools sort"}}

depth_s1{{"Description: Calculate alignment depth \n\n Main tools: Samtools \n\n Commands: samtools depth"}}
depth_s2{{"Description: Calculate alignment depth \n\n Main tools: Samtools \n\n Commands: samtools depth"}}
depth_s3{{"Description: Calculate alignment depth \n\n Main tools: Samtools \n\n Commands: samtools depth"}}
depth_s4{{"Description: Calculate alignment depth \n\n Main tools: Samtools \n\n Commands: samtools depth"}}

snp_indel_calling_s1{{"Description: SNP/indel variant calling \n\n Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) \n\n Commands: run_clair3.sh or pbrun deepvariant"}}
snp_indel_calling_s2{{"Description: SNP/indel variant calling \n\n Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) \n\n Commands: run_clair3.sh or pbrun deepvariant"}}
snp_indel_calling_s3{{"Description: SNP/indel variant calling \n\n Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) \n\n Commands: run_clair3.sh or pbrun deepvariant"}}
snp_indel_calling_s4{{"Description: SNP/indel variant calling \n\n Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) \n\n Commands: run_clair3.sh or pbrun deepvariant"}}

snp_indel_phasing_s1{{"Description: SNP/indel phasing \n\n Main tools: WhatsHap \n\n Commands: whatshap phase"}}
snp_indel_phasing_s2{{"Description: SNP/indel phasing \n\n Main tools: WhatsHap \n\n Commands: whatshap phase"}}
snp_indel_phasing_s3{{"Description: SNP/indel phasing \n\n Main tools: WhatsHap \n\n Commands: whatshap phase"}}
snp_indel_phasing_s4{{"Description: SNP/indel phasing \n\n Main tools: WhatsHap \n\n Commands: whatshap phase"}}

haplotagging_s1{{"Description: Haplotagging bams \n\n Main tools: WhatsHap \n\n Commands: whatshap haplotag"}}
haplotagging_s2{{"Description: Haplotagging bams \n\n Main tools: WhatsHap \n\n Commands: whatshap haplotag"}}
haplotagging_s3{{"Description: Haplotagging bams \n\n Main tools: WhatsHap \n\n Commands: whatshap haplotag"}}
haplotagging_s4{{"Description: Haplotagging bams \n\n Main tools: WhatsHap \n\n Commands: whatshap haplotag"}}

sv_calling_s1{{"Description: Structural variant calling \n\n Main tools: Sniffles2 and/or cuteSV \n\n Commands: sniffles and/or cuteSV"}}
sv_calling_s2{{"Description: Structural variant calling \n\n Main tools: Sniffles2 and/or cuteSV \n\n Commands: sniffles and/or cuteSV"}}
sv_calling_s3{{"Description: Structural variant calling \n\n Main tools: Sniffles2 and/or cuteSV \n\n Commands: sniffles and/or cuteSV"}}
sv_calling_s4{{"Description: Structural variant calling \n\n Main tools: Sniffles2 and/or cuteSV \n\n Commands: sniffles and/or cuteSV"}}

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
- [Clair3](https://github.com/HKU-BAL/Clair3) OR [DeepVariant](https://github.com/google/deepvariant) (wrapped in [NVIDIA Parabricks](https://docs.nvidia.com/clara/parabricks/latest/))
- [WhatsHap](https://github.com/whatshap/whatshap)
- [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) AND/OR [cuteSV](https://github.com/tjiangHIT/cuteSV)
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
        - Java version 17.0.2

## Run it!

See a walkthrough for how to [run pipeface on NCI](./docs/run_on_nci.md).
