# Pipeface setup

- [Pipeface setup](#pipeface-setup)
  - [1. Get pipeline](#1-get-pipeline)
  - [2. Get pipeline inputs](#2-get-pipeline-inputs)
    - [Reference genome](#reference-genome)
      - [hg38](#hg38)
      - [hs1](#hs1)
    - [Tandem repeat call regions file (if running tandem repeat calling)](#tandem-repeat-call-regions-file-if-running-tandem-repeat-calling)
      - [hg38](#hg38-1)
      - [hs1](#hs1-1)
    - [Somalier sites file (if running relatedness check)](#somalier-sites-file-if-running-relatedness-check)
      - [hg38](#hg38-2)
      - [hs1](#hs1-2)
    - [Clair3 models (if running clair3)](#clair3-models-if-running-clair3)
      - [ONT](#ont)
      - [Pacbio HiFi revio](#pacbio-hifi-revio)
  - [3. Modify in\_data\_pipeface.csv](#3-modify-in_data_pipefacecsv)
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

> [!NOTE]
> Variant annotation is only available for hg38.

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

### Tandem repeat call regions file (if running tandem repeat calling)

> [!NOTE]
> You can create a BED file defining the tandem repeats regions you wish to call, alternatively you can use the catalogs below.

#### hg38

Get a copy of the Broad Institute tandem repeat catalog

```bash
wget https://github.com/broadinstitute/tandem-repeat-catalog/releases/download/v1.0.2/variation_clusters_and_isolated_TRs_v1.0.2.hg38.TRGT.bed.gz
```

Check download was successful by checking md5sum

```bash
md5sum variation_clusters_and_isolated_TRs_v1.0.2.hg38.TRGT.bed.gz
```

Expected md5sum

```txt
d50345a1967c507bcdd3cf35c4db27d0  variation_clusters_and_isolated_TRs_v1.0.2.hg38.TRGT.bed.gz
```

gunzip

```bash
gunzip variation_clusters_and_isolated_TRs_v1.0.2.hg38.TRGT.bed.gz
```

Prepare file for LongTR

```bash
cat variation_clusters_and_isolated_TRs_v1.0.2.hg38.TRGT.bed | sed 's/ID.*MOTIFS=//' | sed 's/;.*//' | awk 'length($4) > 1' | awk '$3 - $2 <= 1000' > variation_clusters_and_isolated_TRs_v1.0.2.hg38.TRGT.longtr.bed
```

#### hs1

Get a copy of the Broad Institute tandem repeat catalog

```bash
curl https://zenodo.org/records/14597629/files/chm13.v2.bed.gz?download=1 -o chm13.v2.bed.gz
```

Check download was successful by checking md5sum

```bash
md5sum chm13.v2.bed.gz
```

Expected md5sum

```txt
51d118a5c70b63692f560077cbc0fa10  chm13.v2.bed.gz
```

gunzip

```bash
gunzip chm13.v2.bed.gz
```

Prepare file for LongTR

```bash
cut -f1-4 chm13.v2.bed | awk 'length($4) > 1' | awk '$3 - $2 <= 1000' > chm13.v2.longtr.bed
```

### Somalier sites file (if running relatedness check)

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

## 3. Modify in_data_pipeface.csv

### Singleton mode

Specify the sample ID, family ID, family position, file path to the data, data type, file path to regions of interest BED file and file path to clair3 model for each data to be processed. Eg:

```csv
sample_id,family_id,family_position,file,data_type,regions_of_interest,clair3_model
sample_01,NONE,NONE,/path/to/sample_01_1.fastq.gz,ont,/path/to/regions.bed,/path/to/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_01,NONE,NONE,/path/to/sample_01_2.fastq.gz,ont,/path/to/regions.bed,/path/to/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_02,NONE,NONE,/path/to/sample_02.fastq,ont,/path/to/regions.bed,/path/to/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_03,NONE,NONE,/path/to/sample_03.bam,ont,NONE,/path/to/clair3_models/ont/r1041_e82_400bps_sup_v420/
sample_04,NONE,NONE,/path/to/sample_04_1.bam,pacbio,NONE,/path/to/clair3_models/hifi_revio/
sample_04,NONE,NONE,/path/to/sample_04_2.bam,pacbio,NONE,/path/to/clair3_models/hifi_revio/
```

> [!NOTE]
> - In singleton mode, `family_id` is only used to define the output directory structure.
> - Files with the same value in the `sample_id` column will be merged, this is used to handle multiple sequencing runs of the same sample.

Requirements:

- entries in the `data_type` column must be either 'ont' or 'pacbio' (as appropriate)
- if `in_data_format` is `ubam_fastq`, entries in the `file` column must have a file extension of '.bam', '.fastq.gz' or '.fastq', and multiple entries for a given `sample_id` must share the same extension
- if `in_data_format` is `aligned_bam`, entries in the `file` column must be indexed BAM files (a `.bai` index must exist alongside each BAM)
- for entries in the `file` column, files containing methylation data should be provided in uBAM/aligned BAM format (and not FASTQ format)
- set `family_id` to 'NONE' if not required
- `family_position` can be any value (set to 'NONE' if not required)
- set `regions_of_interest` to 'NONE' if not required
- set `clair3_model` to the path of an appropriate Clair3 model when clair3 is selected as the SNP/indel caller, otherwise set to 'NONE'

### Duo/Trio mode

Specify the sample ID, family ID, family position, file path to the data, data type, file path to regions of interest BED file and file path to clair3 model for each data to be processed. Eg:

```csv
sample_id,family_id,family_position,file,data_type,regions_of_interest,clair3_model
sample_01,family01,proband,/path/to/sample_01_1.bam,ont,NONE,NONE
sample_01,family01,proband,/path/to/sample_01_2.bam,ont,NONE,NONE
sample_02,family01,father,/path/to/sample_02.bam,ont,NONE,NONE
sample_03,family01,mother,/path/to/sample_03.bam,ont,NONE,NONE
sample_04,family02,proband,/path/to/sample_04.bam,ont,NONE,NONE
sample_05,family02,father,/path/to/sample_05.bam,ont,NONE,NONE
sample_06,family02,mother,/path/to/sample_06.bam,ont,NONE,NONE
```

> [!NOTE]
> - In duo/trio mode, `family_id` and `family_position` are required for defining the joint SNP/indel gVCF merging, the SV VCF merging and the joint somalier relatedness/quality control checks.
> - Files with the same value in the `sample_id` column will be merged, this is used to handle multiple sequencing runs of the same sample.

Requirements:

- `family_id` must not be 'NONE'
- entries in the `data_type` column must be either 'ont' or 'pacbio' (as appropriate) and must be the same for a given `family_id`
- if `in_data_format` is `ubam_fastq`, entries in the `file` column must have a file extension of '.bam', '.fastq.gz' or '.fastq', and multiple entries for a given `sample_id` must share the same extension
- if `in_data_format` is `aligned_bam`, entries in the `file` column must be indexed BAM files (a `.bai` index must exist alongside each BAM)
- for entries in the `file` column, files containing methylation data should be provided in uBAM/aligned BAM format (and not FASTQ format)
- provide all entries for a given `sample_id` the same `family_id`
- in duo mode, exactly 2 unique `sample_id` values are required per `family_id`, with a `proband` and either a `father` or `mother` in the `family_position` column
- in trio mode, exactly 3 unique `sample_id` values are required per `family_id`, with a `proband`, `father` and `mother` in the `family_position` column
- set `regions_of_interest` to 'NONE' if not required
- set `clair3_model` to the path of an appropriate Clair3 model when clair3 is selected as the SNP/indel caller, otherwise set to 'NONE'

## 4. Modify parameters_pipeface.json

Specify the path to `in_data_pipeface.csv`. Eg:

```json
    "in_data": "/path/to/in_data_pipeface.csv",
```

Specify the input data format ('ubam_fastq' or 'aligned_bam'). Eg:

```json
    "in_data_format": "ubam_fastq",
```

> [!NOTE]
> - If you provide an aligned BAM and set `in_data_format` to `aligned_bam`, the pipeline will start from post-alignment processes.
> - If you provide an aligned BAM but set `in_data_format` to `ubam_fastq`, the data will start from the beginning and the aligned BAM will be re-aligned.
> - Providing an aligned BAM assumes that the file was generated with minimap2 and the minimap2 `-Y` flag was used (soft clipping for supplementary alignments).

Specify the path to the reference genome and its index. Eg:

```json
    "ref": "/path/to/hg38.fa",
    "ref_index": "/path/to/hg38.fa.fai",
```

Optionally turn on haploid-aware mode. Eg:

```json
    "haploidaware": "yes",
    "sex": "XY",
    "parbed": "/path/to/par.bed",
```

*OR*

```json
    "haploidaware": "no",
    "sex": "NONE",
    "parbed": "NONE",
```

> [!NOTE]
> - Haploid-aware mode is only available for singleton XY samples.
> - Haploid-aware mode requires both chrX and chrY to be present in the reference genome and, if provided, in the `regions_of_interest` file.

Optionally specify the path to the tandem repeat bed file (used by the SV caller to improve SV calling in tandem repeat regions). Set to 'NONE' if not required. Eg:

```json
    "tandem_repeat": "/path/to/tandem_repeat.bed",
```

*OR*

```json
    "tandem_repeat": "NONE",
```

Specify the mode to run the pipeline in ('singleton', 'duo' or 'trio') and the SNP/indel caller to use ('clair3', 'deepvariant' or 'deeptrio'). Eg:

```json
    "mode": "singleton",
    "snp_indel_caller": "deepvariant",
```

*OR*

```json
    "mode": "singleton",
    "snp_indel_caller": "clair3",
```

*OR*

```json
    "mode": "duo",
    "snp_indel_caller": "deepvariant",
```

*OR*

```json
    "mode": "duo",
    "snp_indel_caller": "clair3",
```

*OR*

```json
    "mode": "trio",
    "snp_indel_caller": "deeptrio",
```

*OR*

```json
    "mode": "trio",
    "snp_indel_caller": "clair3",
```

> [!NOTE]
> - Running DeepVariant/DeepTrio on ONT data assumes r10 data.
> - In singleton and duo mode, the SNP/indel caller must be 'clair3' or 'deepvariant'.
> - In trio mode, the SNP/indel caller must be 'clair3' or 'deeptrio'.

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

Optionally specify a threshold for the mapping quality (MAPQ) filter for structural variant calls. Set to 'NONE' to use default thresholds. Maximum value is 60. Eg:

```json
    "sv_mapq": "NONE",
```

*OR*

```json
    "sv_mapq": "60",
```

> [!NOTE]
> If you intend to merge the output SV VCF's with many samples in popface, it's recommended to use MAPQ 60 to allow the SV merging in popface to scale to a large number of samples (for example 500-1000 samples).

Specify whether variant annotation should be carried out ('yes' or 'no'). Eg:

```json
    "annotate": "yes",
```

*OR*

```json
    "annotate": "no",
```

> [!NOTE]
> Variant annotation is only available for hg38.

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

> [!NOTE]
> Processing base modifications assumes base modifications are present in the input data and the input data is in unaligned BAM (uBAM) format.

Optionally run tandem repeat calling and specify the path to an appropriate tandem repeat regions bed file (used by TRGT and LongTR to define the tandem repeat regions to genotype). Set to 'NONE' if not required. Eg:

```json
    "tr_calling": "yes",
    "tr_call_regions": "/path/to/variation_clusters_and_isolated_TRs_v1.0.2.hg38.TRGT.longtr.bed",
```

*OR*

```json
    "tr_calling": "no",
    "tr_call_regions": "NONE",
```

Optionally run relatedness checks and specify the path to an appropriate somalier sites file. Set to 'NONE' if not required. Eg:

```json
    "check_relatedness": "yes",
    "sites": "/path/to/sites.hg38.vcf.gz",
```

*OR*

```json
    "check_relatedness": "no",
    "sites": "NONE",
```

> [!NOTE]
> - In singleton mode, checking relatedness will produce a somalier extracted file.
> - In duo/trio mode, checking relatedness will additionally run joint relatedness and quality control checks.

Specify the directory in which to write the pipeline outputs. Eg:

```json
    "outdir": "/path/to/results/"
```
