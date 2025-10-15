# Setup

- [Setup](#setup)
  - [1. Get pipeline](#1-get-pipeline)
  - [2. Get pipeline inputs](#2-get-pipeline-inputs)
    - [Reference genome](#reference-genome)
      - [hg38](#hg38)
      - [hs1](#hs1)
    - [Somalier sites file (if running relatedness check)](#somalier-sites-file-if-running-relatedness-check)
      - [hg38](#hg38-1)
      - [hs1](#hs1-1)
    - [Clair3 models (if running clair3)](#clair3-models-if-running-clair3)
      - [ONT](#ont)
      - [Pacbio HiFi revio](#pacbio-hifi-revio)
  - [3. Modify in\_data.csv](#3-modify-in_datacsv)
    - [Singleton mode](#singleton-mode)
    - [Duo/Trio mode](#duotrio-mode)
  - [4. Modify parameters\_pipeface.json](#4-modify-parameters_pipefacejson)

## 1. Get pipeline

```bash
git clone https://github.com/leahkemp/pipeface.git
cd pipeface
```

## 2. Get pipeline inputs

### Reference genome

> **_Note:_** Variant annotation is only available for hg38

#### hg38

Get a copy of the hg38 reference genome

```bash
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz .
```

Check download was successful by checking md5sum

```bash
md5sum hg38.analysisSet.fa.gz
```

Expected md5sum

```txt
6d3c82e1e12b127d526395294526b9c8  hg38.analysisSet.fa.gz
```

gunzip and build index

```bash
gunzip hg38.analysisSet.fa.gz
samtools faidx hg38.analysisSet.fa
```

#### hs1

Get a copy of the hs1 reference genome

```bash
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz .
```

Check download was successful by checking md5sum

```bash
md5sum hs1.fa.gz
```

Expected md5sum

```txt
a493d5402cc86ecc3f54f6346d980036  hs1.fa.gz
```

gunzip and build index

```bash
gunzip hs1.fa.gz
samtools faidx hs1.fa
```

### Somalier sites file (if running relatedness check)

> **_Note:_** checking relatedness is only available for duo/trio mode

#### hg38

```bash
wget -O sites.hg38.v0.2.19.vcf.gz https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz
```

#### hs1

```bash
wget -O sites.chm13v2.T2T.v0.2.19.vcf.gz https://github.com/brentp/somalier/files/9954286/sites.chm13v2.T2T.vcf.gz
```

### Clair3 models (if running clair3)

#### ONT

Clone the Rerio github repository

```bash
git clone https://github.com/nanoporetech/rerio
```

Get a copy of the clair3 models

```bash
python3 rerio/download_model.py --clair3
```

#### Pacbio HiFi revio

Get a copy of the clair3 models

```bash
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/hifi_revio.tar.gz
```

Untar

```bash
tar -xvf hifi_revio.tar.gz
```

## 3. Modify in_data.csv

### Singleton mode

Specify the sample ID, family ID (optional), file path to the data, data type, file path to regions of interest bed file (optional) and file path to clair3 model (if running Clair3) for each data to be processed. Eg:

```csv
sample_id,family_id,family_position,file,data_type,regions_of_interest,clair3_model
sample_01,,,/path/to/PGXXXX240090.fastq.gz,ont,/path/to/regions.bed,/path/to/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_01,,,/path/to/PGXXXX240091.fastq.gz,ont,/path/to/regions.bed,/path/to/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_02,,,/path/to/PGXXXX240092.fastq,ont,/path/to/regions.bed,/path/to/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_03,,,/path/to/PGXXOX240065.bam,ont,NONE,/path/to/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_04,,,/path/to/m84088_240403_023825_s1.hifi_reads.bc2034.bam,pacbio,NONE,/path/to/clair3_models/hifi_revio/
sample_04,,,/path/to/m84088_240403_043745_s2.hifi_reads.bc2035.bam,pacbio,NONE,/path/to/clair3_models/hifi_revio/
```

> **_Note:_** In singleton mode, `family_id` will only used to organise the output files into subdirectories of `family_id` (if provided)

### Duo/Trio mode

Specify the sample ID, family ID, family position, file path to the data, data type, file path to regions of interest bed file (optional) and file path to clair3 model (if running Clair3) for each data to be processed. Eg:

```csv
sample_id,family_id,family_position,file,data_type,regions_of_interest,clair3_model
sample_01,family01,proband,/path/to/PGXXOX240065.bam,ont,NONE,NONE
sample_01,family01,proband,/path/to/PGXXOX240066.bam,ont,NONE,NONE
sample_02,family01,father,/path/to/PGXXOX240067.bam,ont,NONE,NONE
sample_03,family01,mother,/path/to/PGXXOX240068.bam,ont,NONE,NONE
sample_04,family02,proband,/path/to/PGXXOX240069.bam,ont,NONE,NONE
sample_05,family02,father,/path/to/PGXXOX240070.bam,ont,NONE,NONE
sample_04,family02,mother,/path/to/PGXXOX240071.bam,ont,NONE,NONE
```

> **_Note:_** In duo/trio mode, `family_id` and `family_position` are used to define the joint SNP/indel calling/merging

> **_Note:_** In duo mode, a `proband`, and either a `father` or `mother` must be defined in the `family_position` column for every `family_id`

> **_Note:_** In trio mode, a `proband`, `father` and `mother` must be defined in the `family_position` column for every `family_id`

> **_Note:_** Files with the same value in the `sample_id` column will be merged, this is used to handle multiple sequencing runs of the same sample

Requirements:

- leave `family_id` and `family_position` empty if not required
- please provide all entries for a given `sample_id` the same `family_id` (this is currently not error checked)
- set `regions_of_interest` to 'NONE' if not required
- similarly, set `clair3_model` to 'NONE' if not required (ie. if you have not selected clair3 as the SNP/indel caller)
- provide full file paths
- multiple entries for a given `sample_id` are required to have the same file extension in the `file` column (eg. '.bam', '.fastq.gz' or '.fastq')
- for entries in the `file` column, the file extension must be either '.bam', '.fastq.gz' or '.fastq' (as appropriate)
- for entries in the `file` column, files containing methylation data should be provided in uBAM format (and not FASTQ format)
- entries in the `data_type` column must be either 'ont' or 'pacbio' (as appropriate)

## 4. Modify parameters_pipeface.json

Specify the path to `in_data.csv`. Eg:

```json
    "in_data": "/path/to/in_data.csv",
```

Specify the input data format ('ubam_fastq'). Eg:

```json
    "in_data_format": "ubam_fastq",
```

Specify the path to the reference genome and it's index. Eg:

```json
    "ref": "/path/to/hg38.fa",
    "ref_index": "/path/to/hg38.fa.fai",
```

Optionally turn on haploid-aware mode (for XY samples only). Eg:

```json
    "haploidaware": "yes",
    "sex": "XY",
    "parbed": "/path/to/par.bed",
```

*OR*

```json
    "haploidaware": "no",
    "sex": "NONE",
    "parbed": "NONE"
```

Optionally specify the path to the tandem repeat bed file. Set to 'NONE' if not required. Eg:

```json
    "tandem_repeat": "/path/to/tandem_repeat.bed",
```

Specify the mode to run the pipeline in ('singleton', 'duo' or 'trio'). Eg:


```json
    "mode": "singleton",
```

Specify the SNP/indel caller to use ('clair3', 'deepvariant' or 'deeptrio'). Eg:

```json
    "snp_indel_caller": "deepvariant",
```

> **_Note:_** Running DeepVariant/DeepTrio on ONT data assumes r10 data

> **_Note:_** In singleton mode, Clair3 and DeepVariant is available

> **_Note:_** In duo mode, only DeepVariant is available

> **_Note:_** In trio mode, only DeepTrio is available

Specify the SV caller to use ('sniffles', 'cutesv' or 'both'). Eg:

```json
    "sv_caller": "sniffles",
```

Specify whether variant annotation should be carried out ('yes' or 'no'). Eg:

```json
    "annotate": "yes",
```

> **_Note:_** variant annotation is only available for hg38

Specify whether alignment depth should be calculated ('yes' or 'no'). Eg:

```json
    "calculate_depth": "yes",
```

Specify whether base modifications should be analysed ('yes' or 'no'). Eg:

```json
    "analyse_base_mods": "yes",
```

> **_Note:_** processing base modifications assume base modifications are present in the input data and the input data is in unaligned BAM (uBAM) format

Optionally run relatedness checks and specify the path to an appropriate somalier sites file. Set to 'NONE' if not required. Eg:

```json
    "check_relatedness": "yes",
    "sites": "/path/to/sites.hg38.vcf.gz",
```

*OR*

```json
    "check_relatedness": "no",
    "sites": "NONE"
```

> **_Note:_** checking relatedness is only available for duo/trio mode

Specify the directory in which to write the pipeline outputs (please provide a full path). Eg:

```json
    "outdir": "/path/to/results/"
```

