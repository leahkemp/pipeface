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

input_data("ONT fastq.gz \n ONT fastq \n ONT uBAM \n pacbio HiFi fastq.gz \n pacbio HiFi fastq \n pacbio HiFi uBAM")
minimap2{{"minimap2"}}
snp_indel_caller{{"clair3 OR deepvariant"}}
whatshap{{"whatshap"}}
sv_caller{{"sniffles OR cuteSV"}}

input_data-.->minimap2-.->snp_indel_caller-.->whatshap-.->sv_caller
minimap2-.->whatshap

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

pacbio_data_f5("Pacbio HiFi fastq.gz \n (sample 4)")
pacbio_data_f6("Pacbio HiFi fastq \n (sample 5)")
pacbio_data_f7("Pacbio HiFi fastq \n (sample 5)")
pacbio_data_f8("Pacbio HiFi uBAM \n (sample 6)")
pacbio_data_f9("Pacbio HiFi uBAM \n (sample 6)")

merge_m1{{"merge runs"}}
merge_m2{{"merge runs"}}
merge_m3{{"merge runs"}}

bam_to_fastq_s3{{"bam to fastq"}}
bam_to_fastq_s6{{"bam to fastq"}}

minimap2_s1{{"minimap2"}}
minimap2_s2{{"minimap2"}}
minimap2_s3{{"minimap2"}}
minimap2_s4{{"minimap2"}}
minimap2_s5{{"minimap2"}}
minimap2_s6{{"minimap2"}}

snp_indel_caller_s1{{"clair3 OR deepvariant"}}
snp_indel_caller_s2{{"clair3 OR deepvariant"}}
snp_indel_caller_s3{{"clair3 OR deepvariant"}}
snp_indel_caller_s4{{"clair3 OR deepvariant"}}
snp_indel_caller_s5{{"clair3 OR deepvariant"}}
snp_indel_caller_s6{{"clair3 OR deepvariant"}}

whatshap_s1{{"whatshap"}}
whatshap_s2{{"whatshap"}}
whatshap_s3{{"whatshap"}}
whatshap_s4{{"whatshap"}}
whatshap_s5{{"whatshap"}}
whatshap_s6{{"whatshap"}}

sv_caller_s1{{""sniffles OR cuteSV""}}
sv_caller_s2{{""sniffles OR cuteSV""}}
sv_caller_s3{{""sniffles OR cuteSV""}}
sv_caller_s4{{""sniffles OR cuteSV""}}
sv_caller_s5{{""sniffles OR cuteSV""}}
sv_caller_s6{{""sniffles OR cuteSV""}}

ont_data_f1-.->merge_m1-.->minimap2_s1-.->snp_indel_caller_s1-.->whatshap_s1-.->sv_caller_s1
ont_data_f2-.->merge_m1
ont_data_f3-.->minimap2_s2-.->snp_indel_caller_s2-.->whatshap_s2-.->sv_caller_s2
ont_data_f4-.->bam_to_fastq_s3-.->minimap2_s3-.->snp_indel_caller_s3-.->whatshap_s3-.->sv_caller_s3

pacbio_data_f5-.->minimap2_s4-.->snp_indel_caller_s4-.->whatshap_s4-.->sv_caller_s4
pacbio_data_f6-.->merge_m2-.->minimap2_s5-.->snp_indel_caller_s5-.->whatshap_s5-.->sv_caller_s5
pacbio_data_f7-.->merge_m2
pacbio_data_f8-.->merge_m3-.->bam_to_fastq_s6-.->minimap2_s6-.->snp_indel_caller_s6-.->whatshap_s6-.->sv_caller_s6
pacbio_data_f9-.->merge_m3

minimap2_s1-.->whatshap_s1
minimap2_s2-.->whatshap_s2
minimap2_s3-.->whatshap_s3
minimap2_s4-.->whatshap_s4
minimap2_s5-.->whatshap_s5
minimap2_s6-.->whatshap_s6

end

```

## Main analyses

- ONT and pacbio HiFi data
- WGS and targeted
- hg38 or hs1 reference genome

## Main tools

- [minimap2](https://github.com/lh3/minimap2)
- [Clair3](https://github.com/HKU-BAL/Clair3) OR [deepvariant](https://github.com/google/deepvariant)
- [whatshap](https://github.com/whatshap/whatshap)
- [sniffles](https://github.com/fritzsedlazeck/Sniffles) OR [cuteSV](https://github.com/tjiangHIT/cuteSV)

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
- Clair3 *OR* deepvariant phased SNP/indel VCF file
- Clair3 *OR* deepvariant phase SNP/indel gVCF file
- Sniffles *OR* cuteSV phased SV VCF file

## Assumptions

- Running pipeline on Australia's [National Computational Infrastructure (NCI)](https://nci.org.au/)
- Access to if89 project on [National Computational Infrastructure (NCI)](https://nci.org.au/)
- Access to pipeline dependencies:
    - [Nextflow and it's java dependency](https://nf-co.re/docs/usage/installation). Validated to run on:
        - Nextflow 24.04.1
        - Java version 17.0.2

## Run it!

See a walkthrough for how to [run pipeface on NCI](./docs/run_on_nci.md).

