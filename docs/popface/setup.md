# Popface setup

- [Popface setup](#popface-setup)
  - [1. Get pipeline](#1-get-pipeline)
  - [2. Modify in\_data\_popface.csv](#2-modify-in_data_popfacecsv)
  - [3. Modify parameters\_popface.json](#3-modify-parameters_popfacejson)

## 1. Get pipeline

```bash
git clone https://github.com/leahkemp/pipeface.git
cd pipeface
```

## 2. Modify in_data_popface.csv

Specify the population ID, sample ID, file path to the gVCF file, file path to the aligned BAM file, file path to the somalier extracted file and the data type for each data to be processed. Eg:

```csv
pop_id,sample_id,gvcf,bam,somalier_file,data_type
pop_01,sample_01,/path/to/sample_01.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_01.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_01.somalier,ont
pop_01,sample_02,/path/to/sample_02.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_02.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_02.somalier,ont
pop_01,sample_03,/path/to/sample_03.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_03.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_03.somalier,ont
pop_01,sample_04,/path/to/sample_04.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_04.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_04.somalier,ont
pop_02,sample_05,/path/to/sample_05.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_05.hg38.minimap2.whatshap.sorted.haplotagged.bam,NONE,pacbio
pop_02,sample_06,/path/to/sample_06.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_06.hg38.minimap2.whatshap.sorted.haplotagged.bam,NONE,pacbio
pop_02,sample_07,/path/to/sample_07.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_07.hg38.minimap2.whatshap.sorted.haplotagged.bam,NONE,pacbio
pop_02,sample_08,/path/to/sample_08.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_08.hg38.minimap2.whatshap.sorted.haplotagged.bam,NONE,pacbio
pop_02,sample_09,/path/to/sample_09.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_09.hg38.minimap2.whatshap.sorted.haplotagged.bam,NONE,pacbio
pop_02,sample_10,/path/to/sample_10.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_10.hg38.minimap2.whatshap.sorted.haplotagged.bam,NONE,pacbio
pop_02,sample_11,/path/to/sample_11.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_11.hg38.minimap2.whatshap.sorted.haplotagged.bam,NONE,pacbio
pop_02,sample_12,/path/to/sample_12.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_12.hg38.minimap2.whatshap.sorted.haplotagged.bam,NONE,pacbio
pop_02,sample_13,/path/to/sample_13.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_13.hg38.minimap2.whatshap.sorted.haplotagged.bam,NONE,pacbio
pop_02,sample_14,/path/to/sample_14.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_14.hg38.minimap2.whatshap.sorted.haplotagged.bam,NONE,pacbio
```

> **_Note:_** `pop_id` is used to define the SNP/indel gVCF merging and somalier relatedness/quality control checks

> **_Note:_** `gvcf`, `bam` and `somalier_file` are all optional (provide 'NONE' if not required)

Requirements:

- all entries in the `sample_id` column must be unique
- all entries in the `sample_id` column must match the sample ID's reocrded in the headers of the associated gVCF's and aligned BAM files
- all entries in the `gvcf`, `bam` and `somalier_file` columns for a given `pop_id` must be either all real files or all set to 'NONE'
- an aligned BAM file must be provided in the `bam` column if a gVCF file is provided in the `gvcf` column
- entries in the `data_type` column must be either 'ont' or 'pacbio' (as appropriate)

## 3. Modify parameters_popface.json

Specify the path to `in_data_popface.csv`. Eg:

```json
    "in_data": "/path/to/in_data_popface.csv",
```

Specify the path to the reference genome and it's index. Eg:

```json
    "ref": "/path/to/hg38.fa",
    "ref_index": "/path/to/hg38.fa.fai",
```

Specify whether variant annotation should be carried out ('yes' or 'no'). Eg:

```json
    "annotate": "yes",
```

> **_Note:_** variant annotation is only available for hg38

Optionally run tandem repeat calling and specify the path to an appropriate tandem repeat regions bed file. Set to 'NONE' if not required. Eg:

```json
    "tr_calling": "yes",
    "tr_call_regions": "/path/to/variation_clusters_and_isolated_TRs_v1.0.2.hg38.TRGT.longtr.bed",
```

*OR*

```json
    "tr_calling": "no",
    "tr_call_regions": "NONE"
```

Specify the directory in which to write the pipeline outputs (please provide a full path). Eg:

```json
    "outdir": "/path/to/results/"
```

