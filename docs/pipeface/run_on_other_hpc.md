# Run pipeface on other HPC

- [Run pipeface on other HPC](#run-pipeface-on-other-hpc)
  - [Assumptions](#assumptions)
  - [1. Download variant databases (optional)](#1-download-variant-databases-optional)
    - [hg38](#hg38)
      - [VEP cache](#vep-cache)
      - [REVEL](#revel)
      - [gnomAD](#gnomad)
      - [ClinVar](#clinvar)
      - [CADD](#cadd)
      - [SpliceAI](#spliceai)
      - [AlphaMissense](#alphamissense)
    - [chm13](#chm13)
      - [VEP GFF](#vep-gff)
      - [gnomAD](#gnomad-1)
      - [ClinVar](#clinvar-1)
      - [SpliceAI](#spliceai-1)
      - [AlphaMissense](#alphamissense-1)
    - [Dfam](#dfam)
  - [2. Modify nextflow\_pipeface.config](#2-modify-nextflow_pipefaceconfig)
  - [3. Get pipeline dependencies](#3-get-pipeline-dependencies)
  - [4. Run pipeface](#4-run-pipeface)
  - [Information](#information)

## Assumptions

- Running on a HPC
- You have access to appropriate GPUs if running DeepVariant/DeepTrio

## 1. Download variant databases (optional)

Download the variant databases if you wish to run the variant annotation component of the pipeline.

> [!NOTE]
> Variant annotation is only available for hg38 and chm13

> [!IMPORTANT]
> These variant annotation databases are third-party resources, each distributed under its own license. Review and comply with each database's license before downloading, using or redistributing them.

### hg38

#### VEP cache

Get a local copy of the VEP cache

```bash
curl -O https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_merged_vep_112_GRCh38.tar.gz
```

Expected md5sum

```bash
51f6fc181a41f7386c85b40219694427  homo_sapiens_merged_vep_112_GRCh38.tar.gz
```

Un-tar

```bash
tar -xzf homo_sapiens_merged_vep_112_GRCh38.tar.gz
```

#### REVEL

Get a local copy of the REVEL database

```bash
curl -o revel-v1.3_all_chromosomes.zip https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip?download=1
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

#### gnomAD

Get a local copy of the gnomAD database

```bash
for i in {1..22} X Y; do curl -O https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr${i}.vcf.bgz; done
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
for chr in {1..22} X Y; do
    echo gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz
done > vcf_list.txt
bcftools concat --naive --file-list vcf_list.txt --output-type z --threads 24 --output gnomad.joint.v4.1.sites.chrall.vcf.gz
```

#### ClinVar

Get a local copy of the ClinVar database

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20240825.vcf.gz
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20240825.vcf.gz.tbi
```

Expected md5sums

```bash
e05111f8e6418ce2898d78f68d39a019  clinvar_20240825.vcf.gz
74fd2cbee0c03af7809a4e9d2960c157  clinvar_20240825.vcf.gz.tbi
```

#### CADD

Get a local copy of the CADD databases

```bash
curl -O https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz
curl -O https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi
curl -O https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz
curl -O https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz.tbi
curl -O https://kircherlab.bihealth.org/download/CADD-SV/v1.1/1000G_phase3_SVs.tsv.gz
curl -O https://kircherlab.bihealth.org/download/CADD-SV/v1.1/1000G_phase3_SVs.tsv.gz.tbi
```

Expected md5sums

```bash
88577a55f1cd519d44e0f415ba248eb9  whole_genome_SNVs.tsv.gz
347df8fac17ea374c4598f4f44c7ce8b  whole_genome_SNVs.tsv.gz.tbi
4b9c685c96d396af4d001c2f7dd9d8f9  gnomad.genomes.r4.0.indel.tsv.gz
85f3d2daa9202c5915c0ce0f1c749a66  gnomad.genomes.r4.0.indel.tsv.gz.tbi
1149f1df5490db9daaf520f957f59f23  1000G_phase3_SVs.tsv.gz
6ab7ffaead8e9f4246953bbd3b35dbd6  1000G_phase3_SVs.tsv.gz.tbi
```

#### SpliceAI

Get a local copy of the SpliceAI database

Manually download from Illumina basespace (https://basespace.illumina.com/s/otSPW8hnhaZR). See [the VEP SpliceAI plugin documentation](https://asia.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#spliceai) for more detail.

#### AlphaMissense

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

### chm13

#### VEP GFF

Get a local copy of the VEP GFF

```bash
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz.tbi
```

Expected md5sums

```bash
00fb576ec456d2d97fbb75fc197a4d88  chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz
5422632fe30edb57c64988b3ed60d406  chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz.tbi
```

#### gnomAD

> [!IMPORTANT]
> The gnomAD data that was lifted over to chm13 is made available by the Genome Aggregation Database consortium under the Open Data Commons Open Database License (ODbL) v1.0. See the bucket [NOTICE.txt](https://s3.ap-southeast-2.wasabisys.com/pipeface-anno/NOTICE.txt) for the full terms, attribution and the modifications made.

Get a local copy of the gnomAD joint (genomes and exomes) database (gnomAD v4.1 joint sites lifted over from hg38 to chm13)

```bash
wget https://s3.ap-southeast-2.wasabisys.com/pipeface-anno/gnomad.joint.v4.1.sites.chm13t2t.vcf.gz
wget https://s3.ap-southeast-2.wasabisys.com/pipeface-anno/gnomad.joint.v4.1.sites.chm13t2t.vcf.gz.tbi
```

Expected md5sums

```bash
5fa2ee21536b85ff400539b275471855  gnomad.joint.v4.1.sites.chm13t2t.vcf.gz
3c796ba0a5e17bff57cd89b94d589500  gnomad.joint.v4.1.sites.chm13t2t.vcf.gz.tbi
```

#### ClinVar

Get a local copy of the ClinVar database

```bash
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/liftover/chm13v2.0_ClinVar20220313.vcf.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/liftover/chm13v2.0_ClinVar20220313.vcf.gz.tbi
```

Expected md5sums

```bash
4cea2750c5990f0e8e0fa11eab995fcd  chm13v2.0_ClinVar20220313.vcf.gz
83333b91de82407e503f7cb1c3dd01a4  chm13v2.0_ClinVar20220313.vcf.gz.tbi
```

#### SpliceAI

> [!IMPORTANT]
> The SpliceAI scores that were lifted over to chm13 are made available by Illumina for academic and not-for-profit research use only. See the bucket [NOTICE.txt](https://s3.ap-southeast-2.wasabisys.com/pipeface-anno/NOTICE.txt) for the full terms, attribution and the modifications made.

Get local copies of the SpliceAI SNV and indel databases (Illumina SpliceAI v1.3 scores lifted over from hg38 to chm13)

```bash
wget https://s3.ap-southeast-2.wasabisys.com/pipeface-anno/spliceai_scores.raw.snv.chm13.vcf.gz
wget https://s3.ap-southeast-2.wasabisys.com/pipeface-anno/spliceai_scores.raw.snv.chm13.vcf.gz.tbi
wget https://s3.ap-southeast-2.wasabisys.com/pipeface-anno/spliceai_scores.raw.indel.chm13.vcf.gz
wget https://s3.ap-southeast-2.wasabisys.com/pipeface-anno/spliceai_scores.raw.indel.chm13.vcf.gz.tbi
```

Expected md5sums

```bash
8e33f4286e6335a96bfd6d027bda7057  spliceai_scores.raw.snv.chm13.vcf.gz
50ba4be4aa4f073e5d1f80cf75e97422  spliceai_scores.raw.snv.chm13.vcf.gz.tbi
b9105deba6662ae980676b612f119207  spliceai_scores.raw.indel.chm13.vcf.gz
8c285451029811cfde55108778d8d500  spliceai_scores.raw.indel.chm13.vcf.gz.tbi
```

#### AlphaMissense

> [!IMPORTANT]
> The AlphaMissense predictions that were lifted over to chm13 are made available by Google DeepMind under the Creative Commons Attribution 4.0 International (CC BY 4.0) license. See the bucket [NOTICE.txt](https://s3.ap-southeast-2.wasabisys.com/pipeface-anno/NOTICE.txt) for the full terms, attribution and the modifications made.

Get a local copy of the AlphaMissense database (lifted over from hg38 to chm13)

```bash
wget https://s3.ap-southeast-2.wasabisys.com/pipeface-anno/AlphaMissense_chm13.tsv.gz
wget https://s3.ap-southeast-2.wasabisys.com/pipeface-anno/AlphaMissense_chm13.tsv.gz.tbi
```

Expected md5sums

```bash
669474b95b93f247ebf35ba4373c7a21  AlphaMissense_chm13.tsv.gz
3c57d7ce3f7a699cb9471d8c9ce5601c  AlphaMissense_chm13.tsv.gz.tbi
```

### Dfam

Get a local copy of the Dfam 3.9 Mammalia partition (used for both hg38 and chm13) and put it in a directory of its own. Eg.

```bash
mkdir -p /path/to/dfam
cd /path/to/dfam
curl -O https://www.dfam.org/releases/Dfam_3.9/families/FamDB/dfam39_full.7.h5.gz
gunzip dfam39_full.7.h5.gz
```

## 2. Modify nextflow_pipeface.config

Specify the paths to your local copies of the variant databases. Eg:

```txt
// hg38 annotation databases
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

// chm13 annotation databases
params.vep_gff = '/path/to/chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz'
params.spliceai_snv_chm13_db = '/path/to/spliceai_scores.raw.snv.chm13.vcf.gz'
params.spliceai_indel_chm13_db = '/path/to/spliceai_scores.raw.indel.chm13.vcf.gz'
params.alphamissense_chm13_db = '/path/to/AlphaMissense_chm13.tsv.gz'
params.gnomad_chm13_db = '/path/to/gnomad.joint.v4.1.sites.chm13t2t.vcf.gz'
params.clinvar_chm13_db = '/path/to/chm13v2.0_ClinVar20220313.vcf.gz'

// sv repeat annotation database
params.dfam_db = '/path/to/dfam/'
```

Modify the rest of the `nextflow_pipeface.config` for your specific HPC/job scheduler: the `executor`, `queue`, `project` and `storage` settings in the `process` block, the `cacheDir` the software containers are pulled to, and per-process resources where needed. Alternatively keep `nextflow_pipeface.config` untouched and put your settings in a small config that starts with `includeConfig 'nextflow_pipeface.config'`, as `nextflow_pipeface_nci.config` does for NCI.

> [!NOTE]
> The 'deepvariant_call_variants' and 'deeptrio_call_variants' processes require access to appropriate GPUs

## 3. Get pipeline dependencies

You'll need access to nextflow and singularity. Tested on:

- nextflow version 25.10.3
- singularity version 3.11.3

## 4. Run pipeface

Run the pipeline. Eg:

```bash
nextflow run pipeface.nf -params-file ./config/parameters_pipeface.json -config ./config/nextflow_pipeface.config
```

Or run a dry run to validate parameters without executing processes. Eg:

```bash
nextflow run pipeface.nf -stub -params-file ./config/parameters_pipeface.json -config ./config/nextflow_pipeface.config
```

If you need to resume a pipeline run, use the `-resume` flag. Eg:

```bash
nextflow run pipeface.nf -resume -params-file ./config/parameters_pipeface.json -config ./config/nextflow_pipeface.config
```

## Information

Please keep in mind that some datasets will require modifications to the default resources (particularly memory, disk usage, walltime). For example WGS data with greater than typical (~30x) sequencing depth.

