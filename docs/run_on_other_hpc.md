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

Download the variant databases if you wish to run the variant annotation component of the pipeline.

> **_Note:_** variant annotation is only available for hg38

### VEP cache

Get a local copy of the VEP cache

```bash
curl -O https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_merged_vep_112_GRCh38.tar.gz
```

Check download was successful by checking md5sum

```bash
md5sum homo_sapiens_merged_vep_112_GRCh38.tar.gz
```

Expected md5sum

```bash
51f6fc181a41f7386c85b40219694427  homo_sapiens_merged_vep_112_GRCh38.tar.gz
```

Un-tar

```bash
tar -xzf homo_sapiens_merged_vep_112_GRCh38.tar.gz
```

### REVEL

Get a local copy of the REVEL database

```bash
curl -o revel-v1.3_all_chromosomes.zip https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip?download=1
```

Check download was successful by checking md5sum

```bash
md5sum revel-v1.3_all_chromosomes.zip
```

Expected md5sum

```bash
3ea2bc33e6b5455fc7e9899da863b5fe  revel-v1.3_all_chromosomes.zip
```

Unzip

```bash
unzip revel-v1.3_all_chromosomes.zip
```

Format

```bash
cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
bgzip new_tabbed_revel.tsv
zgrep -h -v ^#chr new_tabbed_revel.tsv.gz | awk '$3 != "." ' | sort -k1,1 -k3,3n - | cat - > new_tabbed_revel_grch38.tsv
bgzip new_tabbed_revel_grch38.tsv
tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_grch38.tsv.gz
```

### gnomAD

Get a local copy of the gnomAD database. Eg:

```bash
for i in {1..22} X Y; do curl -O https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr${i}.vcf.bgz; done
```

Check download was successful by checking md5sums

```bash
for i in {1..22} X Y; do md5sum gnomad.joint.v4.1.sites.chr${i}.vcf.bgz; done
```

Expected md5sums

```bash
11c62331b0a654fce6a9cd43838de648  gnomad.joint.v4.1.sites.chr1.vcf.bgz
563a8fe6f148621169b0215ac9f19602  gnomad.joint.v4.1.sites.chr2.vcf.bgz
c44f1661bafc15685f1eee593b4886ea  gnomad.joint.v4.1.sites.chr3.vcf.bgz
e6d438b84539c6adad5bc67d6febb33b  gnomad.joint.v4.1.sites.chr4.vcf.bgz
d226db0e0055b5e87e72f4b63158664a  gnomad.joint.v4.1.sites.chr5.vcf.bgz
2820f13a2439ebbdf55066a0c320cdb5  gnomad.joint.v4.1.sites.chr6.vcf.bgz
d1d50e4fa082246a5787eee76c036189  gnomad.joint.v4.1.sites.chr7.vcf.bgz
195f2825e94c9b5e43b34bb2b1ab5c7b  gnomad.joint.v4.1.sites.chr8.vcf.bgz
1c739fb01fd9de816e3fddc958668627  gnomad.joint.v4.1.sites.chr9.vcf.bgz
e2f174f150b5d709d5d7349ac241c438  gnomad.joint.v4.1.sites.chr10.vcf.bgz
b8651f2e5a0aafa23d7fc3406b35bc69  gnomad.joint.v4.1.sites.chr11.vcf.bgz
62795bafd326eae566ef49781a26bc91  gnomad.joint.v4.1.sites.chr12.vcf.bgz
92244327bee6d45973f6077a8134ccd9  gnomad.joint.v4.1.sites.chr13.vcf.bgz
f7a8344b03a4162cb71cdb628dd1e15b  gnomad.joint.v4.1.sites.chr14.vcf.bgz
40c8ab829f973688d2ef891ce5acabb7  gnomad.joint.v4.1.sites.chr15.vcf.bgz
58a5f920fc191b2069126c41278cc077  gnomad.joint.v4.1.sites.chr16.vcf.bgz
aa48657f45d8db7711c69fbf71a25cdc  gnomad.joint.v4.1.sites.chr17.vcf.bgz
80dd729bc61be464c964d3a3bfb0f41a  gnomad.joint.v4.1.sites.chr18.vcf.bgz
1853ca4993ceb25bd6f3a4554173f7cf  gnomad.joint.v4.1.sites.chr19.vcf.bgz
09263d3c29b822760c61607c6398f5c4  gnomad.joint.v4.1.sites.chr20.vcf.bgz
09263d3c29b822760c61607c6398f5c4  gnomad.joint.v4.1.sites.chr21.vcf.bgz
df15a5ea8ae2e3090eae112f548c74ef  gnomad.joint.v4.1.sites.chr22.vcf.bgz
a5288ced0c2fe893fcfae4d2022b9cd9  gnomad.joint.v4.1.sites.chrX.vcf.bgz
7b882f00919d582139acbc116a7a559f  gnomad.joint.v4.1.sites.chrY.vcf.bgz
```

Merge into a single file

```bash
zcat gnomad.joint.v4.1.sites.chr1.vcf.bgz | head -n1000 | grep '#' > gnomad.joint.v4.1.sites.chrall.vcf
for i in {1..22} X Y; do zgrep -v '#' gnomad.joint.v4.1.sites.chr${i}.vcf.gz >> gnomad.joint.v4.1.sites.chrall.vcf; done
bgzip gnomad.joint.v4.1.sites.chrall.vcf
tabix gnomad.joint.v4.1.sites.chrall.vcf.gz
```

### ClinVar

Get a local copy of the ClinVar database

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20240825.vcf.gz
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20240825.vcf.gz.tbi
```

Check download was successful by checking md5sums

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

Get a local copy of the CADD databases

```bash
curl -O https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz
curl -O https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi
curl -O https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz
curl -O https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz.tbi
```

Check download was successful by checking md5sums

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

Get a local copy of the spliceAI database

Manually download from Illumina basespace (https://basespace.illumina.com/s/otSPW8hnhaZR). See [the VEP spliceAI plugin documentation](https://asia.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#spliceai) for more detail).

### AlphaMissense

Get a local copy of the AlphaMissense database

```bash
curl -O https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
```

Check download was successful by checking md5sum

```bash
md5sum AlphaMissense_hg38.tsv.gz
```

Expected md5sums

```bash
9fd167735f16a1b87da6eb3e4c25fcb5  AlphaMissense_hg38.tsv.gz
```

Index

```bash
tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg38.tsv.gz
```

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

> **_Note:_** the 'deepvariant_call_variants' and 'deeptrio_call_variants' processes require access to appropriate GPU's

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

