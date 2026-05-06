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

Specify the population ID, sample ID, relatedness, family position, file path to the gVCF file, file path to the aligned BAM file, file path to the sniffles SV VCF file, file path to the cuteSV SV VCF file, file path to the somalier extracted file and the data type for each data to be processed. Eg:

```csv
pop_id,sample_id,related,family_position,gvcf,bam,sniffles,cutesv,somalier,data_type
pop_01,sample_01,yes,proband,/path/to/sample_01.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_01.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_01.hg38.sniffles.sv.phased.vcf.gz,NONE,/path/to/sample_01.somalier,ont
pop_01,sample_02,yes,father,/path/to/sample_02.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_02.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_02.hg38.sniffles.sv.phased.vcf.gz,NONE,/path/to/sample_02.somalier,ont
pop_01,sample_03,yes,mother,/path/to/sample_03.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_03.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_03.hg38.sniffles.sv.phased.vcf.gz,NONE,/path/to/sample_03.somalier,ont
pop_01,sample_04,yes,grandfather,/path/to/sample_04.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_04.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_04.hg38.sniffles.sv.phased.vcf.gz,NONE,/path/to/sample_04.somalier,ont
pop_02,sample_05,no,NONE,/path/to/sample_05.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_05.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_05.hg38.sniffles.sv.phased.vcf.gz,/path/to/sample_05.hg38.cutesv.sv.vcf.gz,NONE,pacbio
pop_02,sample_06,no,NONE,/path/to/sample_06.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_06.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_06.hg38.sniffles.sv.phased.vcf.gz,/path/to/sample_06.hg38.cutesv.sv.vcf.gz,NONE,pacbio
pop_02,sample_07,no,NONE,/path/to/sample_07.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_07.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_07.hg38.sniffles.sv.phased.vcf.gz,/path/to/sample_07.hg38.cutesv.sv.vcf.gz,NONE,pacbio
pop_02,sample_08,no,NONE,/path/to/sample_08.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_08.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_08.hg38.sniffles.sv.phased.vcf.gz,/path/to/sample_08.hg38.cutesv.sv.vcf.gz,NONE,pacbio
pop_02,sample_09,no,NONE,/path/to/sample_09.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_09.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_09.hg38.sniffles.sv.phased.vcf.gz,/path/to/sample_09.hg38.cutesv.sv.vcf.gz,NONE,pacbio
pop_02,sample_10,no,NONE,/path/to/sample_10.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_10.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_10.hg38.sniffles.sv.phased.vcf.gz,/path/to/sample_10.hg38.cutesv.sv.vcf.gz,NONE,pacbio
pop_02,sample_11,no,NONE,/path/to/sample_11.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_11.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_11.hg38.sniffles.sv.phased.vcf.gz,/path/to/sample_11.hg38.cutesv.sv.vcf.gz,NONE,pacbio
pop_02,sample_12,no,NONE,/path/to/sample_12.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_12.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_12.hg38.sniffles.sv.phased.vcf.gz,/path/to/sample_12.hg38.cutesv.sv.vcf.gz,NONE,pacbio
pop_02,sample_13,no,NONE,/path/to/sample_13.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_13.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_13.hg38.sniffles.sv.phased.vcf.gz,/path/to/sample_13.hg38.cutesv.sv.vcf.gz,NONE,pacbio
pop_02,sample_14,no,NONE,/path/to/sample_14.hg38.deepvariant.snp_indel.g.vcf.gz,/path/to/sample_14.hg38.minimap2.whatshap.sorted.haplotagged.bam,/path/to/sample_14.hg38.sniffles.sv.phased.vcf.gz,/path/to/sample_14.hg38.cutesv.sv.vcf.gz,NONE,pacbio
```

> **_Note:_** `pop_id` is used to define the SNP/indel gVCF merging, the SV VCF merging and the somalier relatedness/quality control checks.

> **_Note:_** `gvcf`, `bam`, `sniffles`, `cutesv`, and `somalier` are all optional (provide 'NONE' if not required).

Requirements:

- `pop_id` must not be 'NONE'
- `sample_id` must not be 'NONE'
- all `sample_id` values must be unique across the entire CSV
- `related` must be 'yes' or 'no'
- `data_type` must be 'ont' or 'pacbio'
- when `related` is 'no', `family_position` must be 'NONE'
- when `related` is 'yes', `family_position` must not be 'NONE'
- when `gvcf` is provided, `bam` must also be provided
- when `sniffles` or `cutesv` is provided, `bam` must also be provided
- `bam` files must have an associated `.bai` index in the same directory
- `sniffles` and `cutesv` VCF files must have an associated `.tbi` index in the same directory
- each `pop_id` must have 2 or more samples
- `data_type` must be identical for all samples within a `pop_id`
- `related` must be identical for all samples within a `pop_id`
- `gvcf`, `bam`, `sniffles`, `cutesv`, and `somalier` entries for a given `pop_id` must be either all real files or all 'NONE' (not mixed)
- when `related` is 'yes', exactly one sample per `pop_id` must have `family_position` = 'proband'
- when `related` is 'yes', the proband must be listed first in the CSV for that `pop_id`
- when `gvcf` files are provided, `snp_indel_caller` must not be 'NONE'
- when all `gvcf` entries are 'NONE', `snp_indel_caller` must be 'NONE'
- when `tr_calling` is 'yes', at least one `bam` must not be 'NONE'

## 3. Modify parameters_popface.json

Specify the path to `in_data_popface.csv`. Eg:

```json
    "in_data": "/path/to/in_data_popface.csv",
```

Specify the path to the reference genome and its index. Eg:

```json
    "ref": "/path/to/hg38.fa",
    "ref_index": "/path/to/hg38.fa.fai",
```

Specify the SNP/indel caller used to generate the gVCF files ('clair3', 'deepvariant', or 'NONE' if no gVCF files are provided). Eg:

```json
    "snp_indel_caller": "deepvariant",
```

*OR*

```json
    "snp_indel_caller": "clair3",
```

*OR*

```json
    "snp_indel_caller": "NONE",
```

Specify whether variant annotation should be carried out ('yes' or 'no'). Eg:

```json
    "annotate": "yes",
```

*OR*

```json
    "annotate": "no",
```

> **_Note:_** Variant annotation is only available for hg38.

Optionally run tandem repeat calling and specify the path to an appropriate tandem repeat regions bed file. Set to 'NONE' if not required. Eg:

```json
    "tr_calling": "yes",
    "tr_call_regions": "/path/to/variation_clusters_and_isolated_TRs_v1.0.2.hg38.TRGT.longtr.bed",
```

*OR*

```json
    "tr_calling": "no",
    "tr_call_regions": "NONE",
```

Specify the directory in which to write the pipeline outputs. Eg:

```json
    "outdir": "/path/to/results/"
```
