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

input_data("Input data: \n\n ONT fastq.gz \n and/or \n ONT fastq \n and/or \n ONT uBAM \n and/or \n pacbio HiFi uBAM")
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

ont_data_f1("Sample 1 \n\n Input data: \n\n ONT fastq.gz")
ont_data_f2("Sample 1 \n\n Input data: \n\n ONT fastq.gz")
pacbio_data_f3("Sample 2 \n\n Input data: \n\n Pacbio HiFi uBAM")
pacbio_data_f4("Sample 2 \n\n Input data: \n\n Pacbio HiFi uBAM")
ont_data_f5("Sample 3 \n\n Input data: \n\n ONT fastq")
ont_data_f6("Sample 4 \n\n Input data: \n\n ONT uBAM")
ont_data_f7("Sample 5 \n\n Input data: \n\n ONT fastq")
ont_data_f8("Sample 5 \n\n Input data: \n\n ONT fastq")

merging_m1{{"Description: merge runs \n\n Main tools: GNU coreutils \n\n Commands: cat \n\n Software versions: GNU coreutils 8.30"}}
merging_m2{{"Description: merge runs \n\n Main tools: Samtools \n\n Commands: samtools merge \n\n Software versions: Samtools 1.19"}}
merging_m3{{"Description: merge runs \n\n Main tools: GNU coreutils \n\n Commands: cat \n\n Software versions: GNU coreutils 8.30, HTSlib 1.16"}}

alignment_s1{{"Description: alignment, sorting \n\n Main tools: Minimap2 and Samtools \n\n Commands: minimap2 and samtools sort \n\n Software verions: Samtools 1.19, Minimap2 2.28-r1209"}}
alignment_s2{{"Description: alignment, sorting \n\n Main tools: Minimap2 and Samtools \n\n Commands: minimap2 and samtools sort \n\n Software verions: Samtools 1.19, Minimap2 2.28-r1209"}}
alignment_s3{{"Description: bam to fastq conversion, alignment, sorting \n\n Main tools: Minimap2 and Samtools \n\n Commands: minimap2 and samtools sort \n\n Software verions: Samtools 1.19, Minimap2 2.28-r1209"}}
alignment_s4{{"Description: bam to fastq conversion, alignment, sorting \n\n Main tools: Minimap2 and Samtools \n\n Commands: minimap2 and samtools sort \n\n Software verions: Samtools 1.19, Minimap2 2.28-r1209"}}
alignment_s5{{"Description: bam to fastq conversion, alignment, sorting \n\n Main tools: Minimap2 and Samtools \n\n Commands: minimap2 and samtools sort \n\n Software verions: Samtools 1.19, Minimap2 2.28-r1209"}}

depth_s1{{"Description: calculate alignment depth \n\n Main tools: Samtools \n\n Commands: samtools depth \n\n Software versions: Samtools 1.19, GNU Coreutils 8.30, GNU Awk 4.2.1"}}
depth_s2{{"Description: calculate alignment depth \n\n Main tools: Samtools \n\n Commands: samtools depth \n\n Software versions: Samtools 1.19, GNU Coreutils 8.30, GNU Awk 4.2.1"}}
depth_s3{{"Description: calculate alignment depth \n\n Main tools: Samtools \n\n Commands: samtools depth \n\n Software versions: Samtools 1.19, GNU Coreutils 8.30, GNU Awk 4.2.1"}}
depth_s4{{"Description: calculate alignment depth \n\n Main tools: Samtools \n\n Commands: samtools depth \n\n Software versions: Samtools 1.19, GNU Coreutils 8.30, GNU Awk 4.2.1"}}
depth_s5{{"Description: calculate alignment depth \n\n Main tools: Samtools \n\n Commands: samtools depth \n\n Software versions: Samtools 1.19, GNU Coreutils 8.30, GNU Awk 4.2.1"}}

snp_indel_calling_s1{{"Description: SNP/indel variant calling \n\n Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) \n\n Commands: run_clair3.sh or pbrun deepvariant \n\n Software versions: Clair3 v1.0.9, GNU Coreutils 8.30 or Clara Parabricks 4.2.1-1.beta4 (equivilant to deepvariant 1.5.0), HTSlib 1.16"}}
snp_indel_calling_s2{{"Description: SNP/indel variant calling \n\n Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) \n\n Commands: run_clair3.sh or pbrun deepvariant \n\n Software versions: Clair3 v1.0.9, GNU Coreutils 8.30 or Clara Parabricks 4.2.1-1.beta4 (equivilant to deepvariant 1.5.0), HTSlib 1.16"}}
snp_indel_calling_s3{{"Description: SNP/indel variant calling \n\n Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) \n\n Commands: run_clair3.sh or pbrun deepvariant \n\n Software versions: Clair3 v1.0.9, GNU Coreutils 8.30 or Clara Parabricks 4.2.1-1.beta4 (equivilant to deepvariant 1.5.0), HTSlib 1.16"}}
snp_indel_calling_s4{{"Description: SNP/indel variant calling \n\n Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) \n\n Commands: run_clair3.sh or pbrun deepvariant \n\n Software versions: Clair3 v1.0.9, GNU Coreutils 8.30 or Clara Parabricks 4.2.1-1.beta4 (equivilant to deepvariant 1.5.0), HTSlib 1.16"}}
snp_indel_calling_s5{{"Description: SNP/indel variant calling \n\n Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) \n\n Commands: run_clair3.sh or pbrun deepvariant \n\n Software versions: Clair3 v1.0.9, GNU Coreutils 8.30 or Clara Parabricks 4.2.1-1.beta4 (equivilant to deepvariant 1.5.0), HTSlib 1.16"}}

snp_indel_phasing_s1{{"Description: SNP/indel phasing \n\n Main tools: WhatsHap \n\n Commands: whatshap phase \n\n Software versions: WhatsHap 2.3, Samtools 1.19, HTSlib 1.16"}}
snp_indel_phasing_s2{{"Description: SNP/indel phasing \n\n Main tools: WhatsHap \n\n Commands: whatshap phase \n\n Software versions: WhatsHap 2.3, Samtools 1.19, HTSlib 1.16"}}
snp_indel_phasing_s3{{"Description: SNP/indel phasing \n\n Main tools: WhatsHap \n\n Commands: whatshap phase \n\n Software versions: WhatsHap 2.3, Samtools 1.19, HTSlib 1.16"}}
snp_indel_phasing_s4{{"Description: SNP/indel phasing \n\n Main tools: WhatsHap \n\n Commands: whatshap phase \n\n Software versions: WhatsHap 2.3, Samtools 1.19, HTSlib 1.16"}}
snp_indel_phasing_s5{{"Description: SNP/indel phasing \n\n Main tools: WhatsHap \n\n Commands: whatshap phase \n\n Software versions: WhatsHap 2.3, Samtools 1.19, HTSlib 1.16"}}

haplotagging_s1{{"Description: haplotagging bams \n\n Main tools: WhatsHap \n\n Commands: whatshap haplotag \n\n Software versions: WhatsHap 2.3, Samtools 1.19, HTSlib 1.16"}}
haplotagging_s2{{"Description: haplotagging bams \n\n Main tools: WhatsHap \n\n Commands: whatshap haplotag \n\n Software versions: WhatsHap 2.3, Samtools 1.19, HTSlib 1.16"}}
haplotagging_s3{{"Description: haplotagging bams \n\n Main tools: WhatsHap \n\n Commands: whatshap haplotag \n\n Software versions: WhatsHap 2.3, Samtools 1.19, HTSlib 1.16"}}
haplotagging_s4{{"Description: haplotagging bams \n\n Main tools: WhatsHap \n\n Commands: whatshap haplotag \n\n Software versions: WhatsHap 2.3, Samtools 1.19, HTSlib 1.16"}}
haplotagging_s5{{"Description: haplotagging bams \n\n Main tools: WhatsHap \n\n Commands: whatshap haplotag \n\n Software versions: WhatsHap 2.3, Samtools 1.19, HTSlib 1.16"}}

sv_calling_s1{{"Description: structural variant calling \n\n Main tools: Sniffles2 and/or cuteSV \n\n Commands: sniffles and/or cuteSV \n\n Software versions: Sniffles2 2.3.3 or cuteSV 1.0.13, HTSlib 1.16"}}
sv_calling_s2{{"Description: structural variant calling \n\n Main tools: Sniffles2 and/or cuteSV \n\n Commands: sniffles and/or cuteSV \n\n Software versions: Sniffles2 2.3.3 or cuteSV 1.0.13, HTSlib 1.16"}}
sv_calling_s3{{"Description: structural variant calling \n\n Main tools: Sniffles2 and/or cuteSV \n\n Commands: sniffles and/or cuteSV \n\n Software versions: Sniffles2 2.3.3 or cuteSV 1.0.13, HTSlib 1.16"}}
sv_calling_s4{{"Description: structural variant calling \n\n Main tools: Sniffles2 and/or cuteSV \n\n Commands: sniffles and/or cuteSV \n\n Software versions: Sniffles2 2.3.3 or cuteSV 1.0.13, HTSlib 1.16"}}
sv_calling_s5{{"Description: structural variant calling \n\n Main tools: Sniffles2 and/or cuteSV \n\n Commands: sniffles and/or cuteSV \n\n Software versions: Sniffles2 2.3.3 or cuteSV 1.0.13, HTSlib 1.16"}}

ont_data_f1-.->merging_m1-.->alignment_s1-.->snp_indel_calling_s1-.->snp_indel_phasing_s1-.->haplotagging_s1-.->sv_calling_s1
ont_data_f2-.->merging_m1
ont_data_f5-.->alignment_s2-.->snp_indel_calling_s2-.->snp_indel_phasing_s2-.->haplotagging_s2-.->sv_calling_s2
ont_data_f6-.->alignment_s3-.->snp_indel_calling_s3-.->snp_indel_phasing_s3-.->haplotagging_s3-.->sv_calling_s3
ont_data_f7-.->merging_m3-.->alignment_s5-.->snp_indel_calling_s5-.->snp_indel_phasing_s5-.->haplotagging_s5-.->sv_calling_s5
ont_data_f8-.->merging_m3

pacbio_data_f3-.->merging_m2-.->alignment_s4-.->snp_indel_calling_s4-.->snp_indel_phasing_s4-.->haplotagging_s4-.->sv_calling_s4
pacbio_data_f4-.->merging_m2

alignment_s1-.->depth_s1
alignment_s2-.->depth_s2
alignment_s3-.->depth_s3
alignment_s4-.->depth_s4
alignment_s5-.->depth_s5

alignment_s1-.->haplotagging_s1
alignment_s2-.->haplotagging_s2
alignment_s3-.->haplotagging_s3
alignment_s4-.->haplotagging_s4
alignment_s5-.->haplotagging_s5

```

*Note. software versions used assumes the default [config/nextflow_pipeface.config](config/nextflow_pipeface.config) file is used*

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

*[See the list of software versions used by this version of pipeface](./docs/main_software_versions.txt)*

## Run it!

See a walkthrough for how to [run pipeface on NCI](./docs/run_on_nci.md).
