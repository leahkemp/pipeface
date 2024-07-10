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

subgraph "Pipeface - overview"

input_data("Input data: \n\n - ONT fastq.gz \n - ONT fastq \n - ONT uBAM \n - pacbio HiFi uBAM")
merging{{"Processes: merge_runs \n\n Description: Merge runs (if needed) \n\n Main tools: Samtools or GNU coreutils \n\n Commands: samtools merge or cat"}}
alignment{{"Processes: minimap2 \n\n Description: bam to fastq conversion (if needed), alignment, sorting \n\n Main tools: Minimap2 and Samtools \n\n Commands: samtools fastq (if needed) minimap2 and samtools sort"}}
snp_indel_calling{{"Processes: clair3 or deepvariant \n\n Description: SNP/indel variant calling \n\n Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) \n\n Commands: run_clair3.sh or pbrun deepvariant"}}
snp_indel_phasing{{"Processes: whatshap_phase_clair3 or whatshap_phase_dv \n\n Description: SNP/indel phasing \n\n Main tools: WhatsHap \n\n Commands: whatshap phase"}}
haplotagging{{"Processes: whatshap_haplotag \n\n Description: Haplotagging bams \n\n Main tools: WhatsHap \n\n Commands: whatshap haplotag"}}
sv_calling{{"Processes: sniffles or cutesv \n\n Description: Structural variant calling \n\n Main tools: Sniffles2 and/or cuteSV \n\n Commands: sniffles and/or cuteSV"}}

input_data-.->merging-.->alignment-.->snp_indel_calling-.->snp_indel_phasing-.->haplotagging-.->sv_calling
alignment-.->haplotagging

end

```

### Detailed

```mermaid
%%{init: {'theme':'dark'}}%%
flowchart LR

subgraph "Pipeface - detailed"

ont_data_f1("Input data: \n\n ONT fastq.gz \n\n (sample 1)")
ont_data_f2("Input data: \n\n ONT fastq.gz \n\n (sample 1)")
ont_data_f3("Input data: \n\n ONT fastq \n\n (sample 2)")
ont_data_f4("Input data: \n\n ONT uBAM \n\n (sample 3)")

pacbio_data_f5("Pacbio HiFi uBAM \n (sample 4)")
pacbio_data_f6("Pacbio HiFi uBAM \n (sample 4)")

merging_m1{{"Processes: merge_runs \n\n Description: Merge runs (if needed) \n\n Main tools: Samtools or GNU coreutils \n\n Commands: cat"}}
merging_m2{{"Processes: merge_runs \n\n Description: Merge runs (if needed) \n\n Main tools: Samtools or GNU coreutils \n\n Commands: samtools merge"}}

alignment_s1{{"Processes: minimap2 \n\n Description: bam to fastq conversion (if needed), alignment, sorting \n\n Main tools: Minimap2 and Samtools \n\n Commands: samtools fastq, minimap2 and samtools sort"}}
alignment_s2{{"Processes: minimap2 \n\n Description: bam to fastq conversion (if needed), alignment, sorting \n\n Main tools: Minimap2 and Samtools \n\n Commands: samtools fastq, minimap2 and samtools sort"}}
alignment_s3{{"Processes: minimap2 \n\n Description: bam to fastq conversion (if needed), alignment, sorting \n\n Main tools: Minimap2 and Samtools \n\n Commands: minimap2 and samtools sort"}}
alignment_s4{{"Processes: minimap2 \n\n Description: bam to fastq conversion (if needed), alignment, sorting \n\n Main tools: Minimap2 and Samtools \n\n Commands: minimap2 and samtools sort"}}

snp_indel_calling_s1{{"Processes: clair3 or deepvariant \n\n Description: SNP/indel variant calling \n\n Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) \n\n Commands: run_clair3.sh or pbrun deepvariant"}}
snp_indel_calling_s2{{"Processes: clair3 or deepvariant \n\n Description: SNP/indel variant calling \n\n Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) \n\n Commands: run_clair3.sh or pbrun deepvariant"}}
snp_indel_calling_s3{{"Processes: clair3 or deepvariant \n\n Description: SNP/indel variant calling \n\n Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) \n\n Commands: run_clair3.sh or pbrun deepvariant"}}
snp_indel_calling_s4{{"Processes: clair3 or deepvariant \n\n Description: SNP/indel variant calling \n\n Main tools: Clair3 or DeepVariant (NVIDIA Parabricks) \n\n Commands: run_clair3.sh or pbrun deepvariant"}}

snp_indel_phasing_s1{{"Processes: whatshap_phase_clair3 or whatshap_phase_dv \n\n Description: SNP/indel phasing \n\n Main tools: WhatsHap \n\n Commands: whatshap phase"}}
snp_indel_phasing_s2{{"Processes: whatshap_phase_clair3 or whatshap_phase_dv \n\n Description: SNP/indel phasing \n\n Main tools: WhatsHap \n\n Commands: whatshap phase"}}
snp_indel_phasing_s3{{"Processes: whatshap_phase_clair3 or whatshap_phase_dv \n\n Description: SNP/indel phasing \n\n Main tools: WhatsHap \n\n Commands: whatshap phase"}}
snp_indel_phasing_s4{{"Processes: whatshap_phase_clair3 or whatshap_phase_dv \n\n Description: SNP/indel phasing \n\n Main tools: WhatsHap \n\n Commands: whatshap phase"}}

haplotagging_s1{{"Processes: whatshap_haplotag \n\n Description: Haplotagging bams \n\n Main tools: WhatsHap \n\n Commands: whatshap haplotag"}}
haplotagging_s2{{"Processes: whatshap_haplotag \n\n Description: Haplotagging bams \n\n Main tools: WhatsHap \n\n Commands: whatshap haplotag"}}
haplotagging_s3{{"Processes: whatshap_haplotag \n\n Description: Haplotagging bams \n\n Main tools: WhatsHap \n\n Commands: whatshap haplotag"}}
haplotagging_s4{{"Processes: whatshap_haplotag \n\n Description: Haplotagging bams \n\n Main tools: WhatsHap \n\n Commands: whatshap haplotag"}}

sv_calling_s1{{"Processes: sniffles or cutesv \n\n Description: Structural variant calling \n\n Main tools: Sniffles2 and/or cuteSV \n\n Commands: sniffles and/or cuteSV"}}
sv_calling_s2{{"Processes: sniffles or cutesv \n\n Description: Structural variant calling \n\n Main tools: Sniffles2 and/or cuteSV \n\n Commands: sniffles and/or cuteSV"}}
sv_calling_s3{{"Processes: sniffles or cutesv \n\n Description: Structural variant calling \n\n Main tools: Sniffles2 and/or cuteSV \n\n Commands: sniffles and/or cuteSV"}}
sv_calling_s4{{"Processes: sniffles or cutesv \n\n Description: Structural variant calling \n\n Main tools: Sniffles2 and/or cuteSV \n\n Commands: sniffles and/or cuteSV"}}

ont_data_f1-.->merging_m1-.->alignment_s1-.->snp_indel_calling_s1-.->snp_indel_phasing_s1-.->haplotagging_s1-.->sv_calling_s1
ont_data_f2-.->merging_m1
ont_data_f3-.->alignment_s2-.->snp_indel_calling_s2-.->snp_indel_phasing_s2-.->haplotagging_s2-.->sv_calling_s2
ont_data_f4-.->alignment_s3-.->snp_indel_calling_s3-.->snp_indel_phasing_s3-.->haplotagging_s3-.->sv_calling_s3

pacbio_data_f5-.->merging_m2-.->alignment_s4-.->snp_indel_calling_s4-.->snp_indel_phasing_s4-.->haplotagging_s4-.->sv_calling_s4
pacbio_data_f6-.->merging_m2

alignment_s1-.->haplotagging_s1
alignment_s2-.->haplotagging_s2
alignment_s3-.->haplotagging_s3
alignment_s4-.->haplotagging_s4

end

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


