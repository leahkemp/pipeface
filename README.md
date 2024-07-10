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

subgraph "Pipeface"

input_data("ONT fastq.gz \n ONT fastq \n ONT uBAM \n pacbio HiFi uBAM")
minimap2{{"minimap2"}}
snp_indel_caller{{"clair3 OR deepvariant"}}
whatshap_phase{{"whatshap phase"}}
whatshap_haplotag{{"whatshap haplotag"}}
sv_caller{{"sniffles AND/OR cuteSV"}}

input_data-.->minimap2-.->snp_indel_caller-.->whatshap_phase-.->whatshap_haplotag-.->sv_caller
minimap2-.->whatshap_haplotag

end

```

### Detailed

```mermaid
%%{init: {'theme':'dark'}}%%
flowchart LR

subgraph "Pipeface"
ont_data_f1("ONT fastq.gz \n (sample 1)")
ont_data_f2("ONT fastq.gz \n (sample 1)")
ont_data_f3("ONT fastq \n (sample 2)")
ont_data_f4("ONT uBAM \n (sample 3)")

pacbio_data_f5("Pacbio HiFi uBAM \n (sample 4)")
pacbio_data_f6("Pacbio HiFi uBAM \n (sample 4)")

merge_m1{{"merge runs"}}
merge_m2{{"merge runs"}}

minimap2_s1{{"Minimap2"}}
minimap2_s2{{"minimap2"}}
minimap2_s3{{"minimap2 \n (bam to fastq on the fly)"}}
minimap2_s4{{"minimap2 \n (bam to fastq on the fly)"}}

snp_indel_caller_s1{{"clair3 OR deepvariant"}}
snp_indel_caller_s2{{"clair3 OR deepvariant"}}
snp_indel_caller_s3{{"clair3 OR deepvariant"}}
snp_indel_caller_s4{{"clair3 OR deepvariant"}}

whatshap_phase_s1{{"whatshap phase"}}
whatshap_phase_s2{{"whatshap phase"}}
whatshap_phase_s3{{"whatshap phase"}}
whatshap_phase_s4{{"whatshap phase"}}

whatshap_haplotag_s1{{"whatshap haplotag"}}
whatshap_haplotag_s2{{"whatshap haplotag"}}
whatshap_haplotag_s3{{"whatshap haplotag"}}
whatshap_haplotag_s4{{"whatshap haplotag"}}

sv_caller_s1{{"sniffles AND/OR cuteSV"}}
sv_caller_s2{{"sniffles AND/OR cuteSV"}}
sv_caller_s3{{"sniffles AND/OR cuteSV"}}
sv_caller_s4{{"sniffles AND/OR cuteSV"}}

ont_data_f1-.->merge_m1-.->minimap2_s1-.->snp_indel_caller_s1-.->whatshap_phase_s1-.->whatshap_haplotag_s1-.->sv_caller_s1
ont_data_f2-.->merge_m1
ont_data_f3-.->minimap2_s2-.->snp_indel_caller_s2-.->whatshap_phase_s2-.->whatshap_haplotag_s2-.->sv_caller_s2
ont_data_f4-.->minimap2_s3-.->snp_indel_caller_s3-.->whatshap_phase_s3-.->whatshap_haplotag_s3-.->sv_caller_s3

pacbio_data_f5-.->merge_m2-.->minimap2_s4-.->snp_indel_caller_s4-.->whatshap_phase_s4-.->whatshap_haplotag_s4-.->sv_caller_s4
pacbio_data_f6-.->merge_m2

minimap2_s1-.->whatshap_haplotag_s1
minimap2_s2-.->whatshap_haplotag_s2
minimap2_s3-.->whatshap_haplotag_s3
minimap2_s4-.->whatshap_haplotag_s4

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

- Aligned, sorted and haplotagged bam's
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


