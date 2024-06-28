# Run pipeface on NCI

- [Run pipeface on NCI](#run-pipeface-on-nci)
  - [1. Get pipeline](#1-get-pipeline)
  - [2. Get pipeline inputs](#2-get-pipeline-inputs)
    - [Reference genome](#reference-genome)
      - [hg38](#hg38)
      - [hs1](#hs1)
    - [Clair3 models](#clair3-models)
      - [ONT](#ont)
      - [Pacbio HiFi revio](#pacbio-hifi-revio)
  - [3. Modify input.csv](#3-modify-inputcsv)
  - [4. Modify nextflow\_pipeface.config](#4-modify-nextflow_pipefaceconfig)
  - [5. Modify parameters\_pipeface.json](#5-modify-parameters_pipefacejson)
  - [6. Get pipeline dependencies](#6-get-pipeline-dependencies)
  - [7. Stub (dry) run](#7-stub-dry-run)
  - [8. Launch pipeline](#8-launch-pipeline)
  - [Advanced](#advanced)

## 1. Get pipeline

```bash
git clone https://github.com/leahkemp/pipeface.git
cd pipeface
```

## 2. Get pipeline inputs

*Please keep in mind that, while hs1 is new, smancy and exciting, hg38 is still the latest GRCh assembly and is better annotated by most projects*

### Reference genome

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
module load samtools/1.19
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
module load samtools/1.19
samtools faidx hs1.fa
```

### Clair3 models

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

## 3. Modify input.csv

Specify the sample ID, data file path, data type, file path to regions of interest bed file (optional) and file path to clair3 model (optional) for each data to be analysed. Eg:

```csv
sample_id,file,data_type,regions_of_interest,clair3_model
sample_01,/g/data/kr68/test_data/PGXXXX240090_minimal.fastq.gz,ont,/g/data/kr68/genome/ReadFish_v9_gene_targets.collapsed.hg38.bed,/g/data/kr68/clair3_models/ont/
sample_01,/g/data/kr68/test_data/PGXXXX240091_minimal.fastq.gz,ont,/g/data/kr68/genome/ReadFish_v9_gene_targets.collapsed.hg38.bed,/g/data/kr68/clair3_models/ont/
sample_02,/g/data/kr68/test_data/PGXXXX240092_minimal.fastq,ont,/g/data/kr68/genome/ReadFish_v9_gene_targets.collapsed.hg38.bed,/g/data/kr68/clair3_models/ont/
sample_03,/g/data/kr68/test_data/PGXXOX240065_minimal.bam,ont,NONE,/g/data/kr68/clair3_models/ont/
sample_04,/g/data/kr68/test_data/m84088_240403_023825_s1.hifi_reads.bc2034_minimal.fastq,pacbio,NONE,/g/data/kr68/clair3_models/hifi_revio/
sample_04,/g/data/kr68/test_data/m84088_240403_043745_s2.hifi_reads.bc2035_minimal.fastq,pacbio,NONE,/g/data/kr68/clair3_models/hifi_revio/
```

Requirements:

- set `regions_of_interest` to 'NONE' if not required
- similarly, set `clair3_model` to 'NONE' if not required (ie. if you have not selected clair3 as the SNP/indel caller)
- provide full file paths
- multiple entries for a given `sample_id` are required to have the same file extension in the `file` column (eg. '.bam', '.fastq.gz' or '.fastq')
- for entries in the `file` column, the file extension must be either '.bam', '.fastq.gz' or '.fastq' (as appropriate)
- entries in the `data_type` column must be either 'ont' or 'pacbio' (as appropriate)

## 4. Modify nextflow_pipeface.config

Modify the NCI project to which to charge the analysis. Eg:

```txt
    project = 'kr68'
```

Modify access to project specific directories. Eg:

```txt
    storage = 'gdata/if89+scratch/kr68+gdata/kr68'
```

> **_Note:_** Don't remove access to if89 gdata (`gdata/if89`). This is required to access environmental modules used in the pipeline

## 5. Modify parameters_pipeface.json

Specify the path to `in_data.csv`. Eg:

```json
    "input": "./config/in_data.csv",
```

Specify the path to the reference genome and it's index. Eg:

```json
    "ref": "./hg38.fa",
    "ref_index": "./hg38.fa.fai",
```

*OR*

```json
    "ref": "./hs1.fa",
    "ref_index": "./hs1.fa.fai",
```

Optionally specify the path to the tandem repeat bed file. Eg:

```json
    "tandem_repeat": "./*.trf.bed",
```

*OR*

```json
    "tandem_repeat": "NONE"
```

Specify the SNP/indel caller to use (either 'clair3' or 'deepvariant'). Eg:

```json
    "snp_indel_caller": "clair3",
```

*OR*

```json
    "snp_indel_caller": "deepvariant",
```

Specify the SV caller to use (either 'sniffles' or 'cutesv'). Eg:

```json
    "sv_caller": "sniffles",
```

*OR*

```json
    "sv_caller": "cutesv",
```

Specify the directory in which to write the pipeline outputs. Eg:

```json
    "outdir": "./results"
```

## 6. Get pipeline dependencies

You may use the centrally installed nextflow environmental module available on NCI to access the nextflow and java dependencies

```bash
module load nextflow/24.04.1
```

## 7. Stub (dry) run

```bash
nextflow run pipeface.nf -stub -params-file ./config/parameters_pipeface.json -config ./config/nextflow_pipeface.config
```

## 8. Launch pipeline

In the same directory, run the pipeline

```bash
nextflow run pipeface.nf -params-file ./config/parameters_pipeface.json -config ./config/nextflow_pipeface.config -with-timeline -with-dag -with-report
```

## Advanced

The resources requested and the queue job's are submitted to may be modified by modifying `./config/nextflow_pipeface.config`.

Similarly, with some coding skills, the environmental modules used by each job in the pipeline may be modified. This means you're able to sub in different versions of software used by the pipeline. However, keep in mind that the pipeline doesn't account for differences in parameterisation between software versions.

This also means this pipeline is adaptable to other HPC's if appropriate environmental modules are included in `./config/nextflow_pipeface.config` (or if you get around to creating a nextflow configuration file pointing to containers for appropriate software before I do) and modify the SGE job scheduler specific configuration if needed.

