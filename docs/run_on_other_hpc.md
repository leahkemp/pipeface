# Run on other HPC

- [Run on other HPC](#run-on-other-hpc)
  - [Assumptions](#assumptions)
  - [1. Download variant databases (optional)](#1-download-variant-databases-optional)
    - [VEP cache](#vep-cache)
    - [REVEL](#revel)
    - [gnomAD](#gnomad)
    - [ClinVar](#clinvar)
    - [CADD](#cadd)
    - [spliceAI](#spliceai)
    - [AlphaMissense](#alphamissense)
  - [2. Modify nextflow\_pipeface\_container.config](#2-modify-nextflow_pipeface_containerconfig)
  - [3. Get pipeline dependencies](#3-get-pipeline-dependencies)
  - [4. Run pipeface](#4-run-pipeface)
  - [Information](#information)

## Assumptions

- Running on a HPC
- You have access to appropriate GPU's if running DeepVariant/DeepTrio

## 1. Download variant databases (optional)

Download the variant databases if you wish to run variant annotation.

> **_Note:_** variant annotation is only available for hg38

### VEP cache

Get a local copy of the VEP cache

```bash
curl -O https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_merged_vep_112_GRCh38.tar.gz
```

```bash
tar -xzf homo_sapiens_merged_vep_112_GRCh38.tar.gz
```

*Add expected md5sum**

### REVEL

TODO

### gnomAD

TODO

### ClinVar

Get a local copy of the ClinVar database

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20240825.vcf.gz
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20240825.vcf.gz.tbi
```

Check download was successful by checking md5sum

```bash
md5sum clinvar_20240825.vcf.gz
md5sum clinvar_20240825.vcf.gz.tbi
```

Expected md5sums

```bash
e05111f8e6418ce2898d78f68d39a019  clinvar_20240825.vcf.gz
74fd2cbee0c03af7809a4e9d2960c157  clinvar_20240825.vcf.gz.tbi
```

### CADD

Get a local copy of the CADD database

```bash
curl -O https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz
curl -O https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi
```

Check download was successful by checking md5sum

```bash
md5sum whole_genome_SNVs.tsv.gz
md5sum whole_genome_SNVs.tsv.gz.tbi
md5sum gnomad.genomes.r4.0.indel.tsv.gz
md5sum gnomad.genomes.r4.0.indel.tsv.gz.tbi
```

Expected md5sums

```bash
88577a55f1cd519d44e0f415ba248eb9  whole_genome_SNVs.tsv.gz
347df8fac17ea374c4598f4f44c7ce8b  whole_genome_SNVs.tsv.gz.tbi
4b9c685c96d396af4d001c2f7dd9d8f9  gnomad.genomes.r4.0.indel.tsv.gz
85f3d2daa9202c5915c0ce0f1c749a66  gnomad.genomes.r4.0.indel.tsv.gz.tbi
```

### spliceAI

TODO

### AlphaMissense

Get a local copy of the AlphaMissense database

```bash
curl -O https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg38.tsv.gz
```

*Grab pre-processed index*

## 2. Modify nextflow_pipeface_container.config

Specify the paths to your local copies of the variant databases. Eg:

```txt
params.vep_db = '/path/to/vep/grch38/'
params.revel_db = '/path/to/new_tabbed_revel_grch38.tsv.gz'
params.gnomad_db = '/path/to/gnomad.joint.v4.1.sites.chrall.vcf.gz'
params.clinvar_db = '/path/to/clinvar_20240825.vcf.gz'
params.cadd_snv_db = '/path/to/whole_genome_SNVs.tsv.gz'
params.cadd_indel_db = '/path/to/gnomad.genomes.r4.0.indel.tsv.gz'
params.cadd_sv_db = '/path/to/1000G_phase3_SVs.tsv.gz'
params.spliceai_snv_db = '/path/to/spliceai_scores.raw.snv.hg38.vcf.gz'
params.spliceai_indel_db = '/path/to/spliceai_scores.raw.indel.hg38.vcf.gz'
params.alphamissense_db = '/path/to/AlphaMissense_hg38.tsv.gz'
```

Modify the rest of the `nextflow_pipeface_container.config` for your specific HPC/job sheduler.

> **_Note:_** the 'deepvariant_call_variants' and 'deeptrio_call_variants' processes require a queue with appropriate GPU's available

## 3. Get pipeline dependencies

You'll need access to nextflow and singularity. Tested on:

- nextflow version 24.04.5
- singularity version 3.11.3

## 4. Run pipeface

For example:

```bash
nextflow run pipeface.nf -params-file ./config/parameters_pipeface.json -config ./config/nextflow_pipeface_container.config
```

## Information

Please keep in mind that some datasets will require modifications to the default resources (particularly memory, disk usage, walltime). For example WGS data with greater than typical (~30x) sequencing depth.

