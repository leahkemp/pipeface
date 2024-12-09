# Run pipeface on NCI

- [Run pipeface on NCI](#run-pipeface-on-nci)
  - [1. Get pipeline](#1-get-pipeline)
  - [2. Get pipeline inputs](#2-get-pipeline-inputs)
    - [Reference genome](#reference-genome)
      - [hg38](#hg38)
      - [hs1](#hs1)
    - [Clair3 models (if running clair3)](#clair3-models-if-running-clair3)
      - [ONT](#ont)
      - [Pacbio HiFi revio](#pacbio-hifi-revio)
    - [mosdepth binary (if running depth calculation)](#mosdepth-binary-if-running-depth-calculation)
    - [pb-CpG-tools binary (if processing pacbio data)](#pb-cpg-tools-binary-if-processing-pacbio-data)
  - [3. Modify in\_data.csv](#3-modify-in_datacsv)
  - [4. Modify nextflow\_pipeface.config](#4-modify-nextflow_pipefaceconfig)
  - [5. Modify parameters\_pipeface.json](#5-modify-parameters_pipefacejson)
  - [6. Get pipeline dependencies](#6-get-pipeline-dependencies)
  - [7. Stub (dry) run](#7-stub-dry-run)
  - [8. Full run](#8-full-run)
  - [Advanced](#advanced)

## 1. Get pipeline

```bash
git clone https://github.com/leahkemp/pipeface.git
cd pipeface
```

## 2. Get pipeline inputs

### Reference genome

> **_Note:_** SNP/indel variant annotation is only available for hg38

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

### mosdepth binary (if running depth calculation)

Get a local copy of the mosdepth v0.3.9 binary

```bash
wget https://github.com/brentp/mosdepth/releases/download/v0.3.9/mosdepth -O mosdepth_0.3.9
chmod +x mosdepth_0.3.9
```

### pb-CpG-tools binary (if processing pacbio data)

Get a local copy of the pb-CpG-tools v2.3.2 binary

```bash
wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v2.3.2/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz
tar -xzf pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz
```

## 3. Modify in_data.csv

Specify the sample ID, file path to the data, data type, file path to regions of interest bed file (optional) and file path to clair3 model (if running Clair3) for each data to be analysed. Eg:

```csv
sample_id,file,data_type,regions_of_interest,clair3_model
sample_01,/g/data/kr68/test_data/PGXXXX240090_minimal.fastq.gz,ont,/g/data/kr68/genome/ReadFish_v9_gene_targets.collapsed.hg38.bed,/g/data/kr68/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_01,/g/data/kr68/test_data/PGXXXX240091_minimal.fastq.gz,ont,/g/data/kr68/genome/ReadFish_v9_gene_targets.collapsed.hg38.bed,/g/data/kr68/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_02,/g/data/kr68/test_data/PGXXXX240092_minimal.fastq,ont,/g/data/kr68/genome/ReadFish_v9_gene_targets.collapsed.hg38.bed,/g/data/kr68/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_03,/g/data/kr68/test_data/PGXXOX240065_minimal.bam,ont,NONE,/g/data/kr68/clair3_models/ont/r1041_e82_400bps_sup_v420/
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
    storage = 'gdata/if89+gdata/xy86+scratch/kr68+gdata/kr68'
```

> **_Note:_** Don't remove access to if89 gdata (`gdata/if89`). This is required to access environmental modules used in the pipeline

> **_Note:_** Don't remove access to xy86 gdata (`gdata/xy86`) if running variant annotation. This is required to access variant annotation databases used in the pipeline

## 5. Modify parameters_pipeface.json

Specify the path to `in_data.csv`. Eg:

```json
    "in_data": "./config/in_data.csv",
```

Specify the input data format ('ubam_fastq'). Eg:

```json
    "in_data_format": "ubam_fastq",
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

Optionally specify the path to the tandem repeat bed file. Set to 'NONE' if not required. Eg:

```json
    "tandem_repeat": "./*.trf.bed",
```

*OR*

```json
    "tandem_repeat": "NONE"
```

Specify the SNP/indel caller to use ('clair3' or 'deepvariant'). Eg:


```json
    "snp_indel_caller": "clair3",
```

*OR*

```json
    "snp_indel_caller": "deepvariant",
```

> **_Note:_** Running DeepVariant on ONT data assumes r10 data

Specify the SV caller to use ('sniffles', 'cutesv' or 'both'). Eg:

```json
    "sv_caller": "sniffles",
```

*OR*

```json
    "sv_caller": "cutesv",
```

*OR*

```json
    "sv_caller": "both",
```

Specify whether variant annotation should be carried out ('yes' or 'no'). Eg:

```json
    "annotate": "no",
```

*OR*

```json
    "annotate": "yes",
```

> **_Note:_** SNP/indel variant annotation is only available for hg38

Specify whether alignment depth should be calculated ('yes' or 'no'). Eg:

```json
    "calculate_depth": "no",
```

*OR*

```json
    "calculate_depth": "yes",
```

Specify the directory in which to write the pipeline outputs (please provide a full path). Eg:

```json
    "outdir": "/g/data/ox63/results"
```

Specify the path to the mosdepth binary (if running depth calculation). Eg:

```json
    "mosdepth_binary": "./mosdepth_0.3.9"
```

*OR*

```json
    "mosdepth_binary": "NONE"
```

Specify the path to the pb-CpG-tools binary (if processing pacbio data). Eg:

```json
    "pbcpgtools_binary": "./pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/"
```

*OR*

```json
    "pbcpgtools_binary": "NONE"
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

## 8. Full run

```bash
nextflow run pipeface.nf -params-file ./config/parameters_pipeface.json -config ./config/nextflow_pipeface.config -with-timeline -with-dag -with-report
```

## Advanced

The resources requested and the queue each process is submitted to may be modified by modifying `./config/nextflow_pipeface.config`.

Similarly, with some coding skills, the environmental modules used by each process in the pipeline may be modified. This means you're able to substitute in different versions of software used by the pipeline. However, keep in mind that the pipeline doesn't account for differences in parameterisation between software versions.

This also means this pipeline is adaptable to other HPC's if appropriate environmental modules are included in `./config/nextflow_pipeface.config` (or if you get around to creating a nextflow configuration file pointing to appropriate containerised software before I do) and modify the job scheduler specific configuration if needed. If you wish to use the variant annotation component of the pipeline, you'll additionally need to create local copies of the variant annotation databases used by the pipeline.

