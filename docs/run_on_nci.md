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
    - [Singleton mode](#singleton-mode)
    - [Cohort mode](#cohort-mode)
  - [4. Modify nextflow\_pipeface.config](#4-modify-nextflow_pipefaceconfig)
  - [5. Modify parameters\_pipeface.json](#5-modify-parameters_pipefacejson)
  - [6. Start persistent session (optional)](#6-start-persistent-session-optional)
  - [7. Get pipeline dependencies](#7-get-pipeline-dependencies)
  - [8. Run pipeface](#8-run-pipeface)
  - [Advanced](#advanced)

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

### Singleton mode

Specify the sample ID, family ID (optional), file path to the data, data type, file path to regions of interest bed file (optional) and file path to clair3 model (if running Clair3) for each data to be analysed. Eg:

```csv
sample_id,family_id,family_position,file,data_type,regions_of_interest,clair3_model
sample_01,,,/g/data/kr68/test_data/PGXXXX240090_minimal.fastq.gz,ont,/g/data/kr68/genome/ReadFish_v9_gene_targets.collapsed.hg38.bed,/g/data/kr68/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_01,,,/g/data/kr68/test_data/PGXXXX240091_minimal.fastq.gz,ont,/g/data/kr68/genome/ReadFish_v9_gene_targets.collapsed.hg38.bed,/g/data/kr68/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_02,,,/g/data/kr68/test_data/PGXXXX240092_minimal.fastq,ont,/g/data/kr68/genome/ReadFish_v9_gene_targets.collapsed.hg38.bed,/g/data/kr68/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_03,,,/g/data/kr68/test_data/PGXXOX240065_minimal.bam,ont,NONE,/g/data/kr68/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_04,,,/g/data/kr68/test_data/m84088_240403_023825_s1.hifi_reads.bc2034_minimal.fastq,pacbio,NONE,/g/data/kr68/clair3_models/hifi_revio/
sample_04,,,/g/data/kr68/test_data/m84088_240403_043745_s2.hifi_reads.bc2035_minimal.fastq,pacbio,NONE,/g/data/kr68/clair3_models/hifi_revio/
```

> **_Note:_** In singleton mode, `family_id` will only used to organise the output files into subdirectories of `family_id` (if provided)

### Cohort mode

Specify the sample ID, family ID, family position, file path to the data, data type, file path to regions of interest bed file (optional) and file path to clair3 model (if running Clair3) for each data to be analysed. Eg:

```csv
sample_id,family_id,family_position,file,data_type,regions_of_interest,clair3_model
sample_01,family01,proband,/g/data/kr68/PGXXOX240065.bam,ont,NONE,NONE
sample_01,family01,proband,/g/data/kr68/PGXXOX240066.bam,ont,NONE,NONE
sample_02,family01,father,/g/data/kr68/PGXXOX240067.bam,ont,NONE,NONE
sample_03,family01,mother,/g/data/kr68/PGXXOX240068.bam,ont,NONE,NONE
sample_04,family02,proband,/g/data/kr68/PGXXOX240069.bam,ont,NONE,NONE
sample_05,family02,father,/g/data/kr68/PGXXOX240070.bam,ont,NONE,NONE
sample_04,family02,mother,/g/data/kr68/PGXXOX240071.bam,ont,NONE,NONE
```

> **_Note:_** In cohort mode, `family_id` and `family_position` are used to define the joint SNP/indel calling

> **_Note:_** In cohort mode, a `proband`, `father` and `mother` must be defined in the `family_position` column for every `family_id`

> **_Note:_** Files with the same value in the `sample_id` column will be merged before analysis, this is used to handle multiple sequencing runs of the same sample

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

## 4. Modify nextflow_pipeface.config

Modify the NCI project to which to charge the analysis. Eg:

```txt
    project = 'kr68'
```

Modify access to project specific directories. Eg:

```txt
    storage = 'gdata/if89+gdata/xy86+scratch/kr68+gdata/kr68+gdata/ox63'
```

> **_Note:_** Don't remove access to if89 gdata (`gdata/if89`) and xy86 gdata (`gdata/xy86`). These are required to access software installs and variant annotation databases used in the pipeline

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

Optionally turn on haploid-aware mode (for XY samples only). Eg:

```json
{
  "haploidaware": "yes",
  "sex": "XY",
  "parbed": "/path/to/par.bed",
}
```

*OR no haploid-aware mode:*

```json
{
  "haploidaware": "no",
  "sex": "NONE",
  "parbed": "NONE"
}
```

Optionally specify the path to the tandem repeat bed file. Set to 'NONE' if not required. Eg:

```json
    "tandem_repeat": "./hg38.analysisSet.trf.bed",
```

*OR*

```json
    "tandem_repeat": "NONE"
```

Specify the SNP/indel caller to use ('clair3', 'deepvariant' or 'deeptrio'). Eg:


```json
    "snp_indel_caller": "clair3",
```

*OR*

```json
    "snp_indel_caller": "deepvariant",
```

*OR*

```json
    "snp_indel_caller": "deeptrio",
```

> **_Note:_** Running DeepVariant/DeepTrio on ONT data assumes r10 data

> **_Note:_** Selecting DeepTrio as the SNP/indel caller initates cohort analysis

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
    "annotate": "yes",
```

*OR*

```json
    "annotate": "no",
```

> **_Note:_** variant annotation is only available for hg38

Specify whether alignment depth should be calculated ('yes' or 'no'). Eg:

```json
    "calculate_depth": "yes",
```

*OR*

```json
    "calculate_depth": "no",
```

Specify whether base modifications should be analysed ('yes' or 'no'). Eg:

```json
    "analyse_base_mods": "yes",
```

*OR*

```json
    "analyse_base_mods": "no",
```

> **_Note:_** these analyses assume base modifications are present in the input data and the input data is in unaligned BAM (uBAM) format

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

## 6. Start persistent session (optional)

Pipeface can be run within a [persistent session](https://opus.nci.org.au/spaces/Help/pages/241927941/Persistent+Sessions...)

## 7. Get pipeline dependencies

You may use the centrally installed nextflow environmental module available on NCI to access the nextflow and java dependencies

```bash
module load nextflow/24.04.4
```

## 8. Run pipeface

```bash
nextflow run pipeface.nf -params-file ./config/parameters_pipeface.json -config ./config/nextflow_pipeface.config -with-timeline -with-dag -with-report
```

## Advanced

The resources requested and the queue each process is submitted to may be modified by modifying [nextflow_pipeface.config](https://github.com/leahkemp/pipeface/blob/main/config/nextflow_pipeface.config). Please keep in mind that some datasets will require modifications to the default resources (particularly memory, disk usage, walltime). For example WGS data with greater than typical (~30x) sequencing depth.

Similarly, with some coding skills, the software installs used by each process in the pipeline may be modified. This means you're able to substitute in different software installs or different versions of software used by the pipeline. However, keep in mind that the pipeline doesn't account for differences in parameterisation between software versions.

This also means this pipeline is portable to other HPC's if appropriate environmental modules are included in [nextflow_pipeface.config](https://github.com/leahkemp/pipeface/blob/main/config/nextflow_pipeface.config) (or if you get around to creating a nextflow configuration file pointing to appropriate containerised software before I do) and modify the job scheduler specific configuration if needed. If you wish to use the variant annotation component of the pipeline, you'll additionally need to create local copies of the variant annotation databases used by the pipeline.

