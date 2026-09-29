# How the chm13 annotation databases were made

- [How the chm13 annotation databases were made](#how-the-chm13-annotation-databases-were-made)
  - [Common setup](#common-setup)
  - [gnomAD](#gnomad)
  - [SpliceAI](#spliceai)
  - [AlphaMissense](#alphamissense)
  - [Retention](#retention)

The chm13 (T2T-CHM13v2.0) annotation databases served from the `pipeface-anno` bucket were lifted over from their hg38 releases. This page records how, so the files can be checked or regenerated. Licences and attribution are in the bucket [NOTICE.txt](https://s3.ap-southeast-2.wasabisys.com/pipeface-anno/NOTICE.txt).

## Common setup

All liftovers use [bcftools](https://github.com/samtools/bcftools) 1.23 built with plugins, the [bcftools +liftover plugin](https://github.com/freeseek/score), the hg38 no-alt reference `hg38.analysisSet.fa`, the chm13 reference `chm13v2.0.fa` (both indexed with `samtools faidx`) and the UCSC hg38 to hs1 (chm13) chain file

```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz
```

## gnomAD

Source: gnomAD v4.1 joint (genomes and exomes) sites release, per chromosome

```bash
for c in chr{1..22} chrX chrY; do
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.${c}.vcf.bgz
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.${c}.vcf.bgz.tbi
done
```

Lift over each chromosome, keeping the rejected records so retention can be checked

```bash
for c in chr{1..22} chrX chrY; do
    bcftools view -Ou gnomad.joint.v4.1.sites.${c}.vcf.bgz \
    | bcftools +liftover -Ou -- -s hg38.analysisSet.fa -f chm13v2.0.fa -c hg38ToHs1.over.chain.gz --reject gnomad.joint.v4.1.sites.${c}.rejected.vcf.gz --reject-type z --write-reject \
    | bcftools sort -Oz -o gnomad.joint.v4.1.sites.${c}.chm13.vcf.gz
    bcftools index -t gnomad.joint.v4.1.sites.${c}.chm13.vcf.gz
done
```

Concatenate in chm13 reference contig order, sort and index

```bash
bcftools concat --naive-force -Oz -o gnomad.joint.v4.1.sites.chm13t2t.unsorted.vcf.gz $(for c in $(cut -f1 chm13v2.0.fa.fai); do ls gnomad.joint.v4.1.sites.${c}.chm13.vcf.gz 2>/dev/null; done)
bcftools sort -m 480G -T ./sort_tmp -Oz -o gnomad.joint.v4.1.sites.chm13t2t.vcf.gz gnomad.joint.v4.1.sites.chm13t2t.unsorted.vcf.gz
tabix -p vcf gnomad.joint.v4.1.sites.chm13t2t.vcf.gz
```

The whole-genome sort needs a machine with several hundred GB of memory, or a smaller `-m` with more temporary files. Only genomic coordinates change. The INFO fields pipeface annotates with (`AF_joint`, `AF_exomes`, `AF_genomes`, `nhomalt_joint`, `nhomalt_exomes`, `nhomalt_genomes`) are carried through unchanged, and records that could not be placed on chm13 are dropped.

## SpliceAI

Source: Illumina SpliceAI v1.3 raw scores, `spliceai_scores.raw.snv.hg38.vcf.gz` and `spliceai_scores.raw.indel.hg38.vcf.gz`, downloaded manually from [Illumina basespace](https://basespace.illumina.com/s/otSPW8hnhaZR). These files need fixing before they can be lifted over: the header does not declare the alt and random contigs that appear in records, and contig names lack the `chr` prefix that the chain file and chm13 reference use.

Rebuild the header from the hg38 reference index (contig names without `chr`, to match the records), reheader, then add the `chr` prefix to every record

```bash
awk '{sub(/^chr/, "", $1); print "##contig=<ID="$1",length="$2">"}' hg38.analysisSet.fa.fai > contigs.txt
bcftools view -h spliceai_scores.raw.snv.hg38.vcf.gz | grep -v '^##contig' | awk -v f=contigs.txt '/^#CHROM/{while((getline l < f) > 0) print l} {print}' > header.txt
bcftools reheader -h header.txt spliceai_scores.raw.snv.hg38.vcf.gz | bcftools view -H | awk 'BEGIN{OFS="\t"} {$1="chr"$1; print}' > records.txt
{ sed 's/ID=\([^,]*\)/ID=chr\1/' header.txt; cat records.txt; } | bgzip > spliceai_scores.raw.snv.hg38.fixed.vcf.gz
tabix -p vcf spliceai_scores.raw.snv.hg38.fixed.vcf.gz
```

Find positions whose REF allele does not match hg38 (per chromosome, in parallel) and merge them into one exclusion list

```bash
for c in chr{1..22} chrX chrY; do
    bcftools norm -r $c -f hg38.analysisSet.fa --check-ref w spliceai_scores.raw.snv.hg38.fixed.vcf.gz 2>&1 >/dev/null | awk '/REF_MISMATCH/{print $NF}' | tr ':' '\t' > mismatches_$c.txt &
done; wait
cat mismatches_chr*.txt | sort -u > mismatch_positions.txt
```

Lift over each chromosome excluding those positions, keeping the rejected records, then concatenate and sort

```bash
for c in chr{1..22} chrX chrY; do
    bcftools view -r $c -T ^mismatch_positions.txt spliceai_scores.raw.snv.hg38.fixed.vcf.gz \
    | bcftools +liftover -Ou -- -s hg38.analysisSet.fa -f chm13v2.0.fa -c hg38ToHs1.over.chain.gz --reject $c.rejected.vcf.gz --reject-type z --write-reject \
    | bcftools sort -Oz -o $c.lifted.vcf.gz -W=tbi
done
bcftools concat --naive-force -Oz -o spliceai_scores.raw.snv.chm13.unsorted.vcf.gz $(for c in $(cut -f1 chm13v2.0.fa.fai); do ls $c.lifted.vcf.gz 2>/dev/null; done)
bcftools sort -Oz -o spliceai_scores.raw.snv.chm13.vcf.gz spliceai_scores.raw.snv.chm13.unsorted.vcf.gz
tabix -p vcf spliceai_scores.raw.snv.chm13.vcf.gz
```

The indel file goes through the same steps to produce `spliceai_scores.raw.indel.chm13.vcf.gz`.

## AlphaMissense

Source: `AlphaMissense_hg38.tsv.gz` from [Google DeepMind](https://github.com/google-deepmind/alphamissense)

```bash
wget https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
```

AlphaMissense is a TSV with one row per transcript, so the same CHROM, POS, REF and ALT can appear on several rows. Convert it to a VCF with one record per variant, collapsing the per-transcript columns (`uniprot_id`, `transcript_id`, `protein_variant`, `am_pathogenicity`, `am_class`) into comma-delimited INFO fields, then bgzip and tabix it as `AlphaMissense_hg38.vcf.gz`.

Drop REF-mismatching records, lift over, sort

```bash
bcftools norm -f hg38.analysisSet.fa --check-ref x AlphaMissense_hg38.vcf.gz \
| bcftools +liftover -Ou -- -s hg38.analysisSet.fa -f chm13v2.0.fa -c hg38ToHs1.over.chain.gz --reject AlphaMissense_chm13.rejected.vcf.gz --reject-type z --write-reject \
| bcftools sort -Oz -o AlphaMissense_chm13.vcf.gz -W=tbi
```

Convert back to the original TSV layout, expanding the comma-delimited INFO fields to one row per transcript again, then compress and index for the VEP plugin

```bash
bgzip AlphaMissense_chm13.tsv
tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_chm13.tsv.gz
```

Only genomic coordinates change. Prediction values are untouched.

## Retention

| Database | Input records | Lifted | Retention |
|---|---|---|---|
| SpliceAI SNV | 3,433,386,137 | 3,424,209,286 | 99.73% |
| SpliceAI indel | 9,155,875,164 | 9,128,689,946 | 99.68% |
| AlphaMissense | 71,034,269 | 70,998,360 | 99.95% |

Losses are REF-mismatching records and records the chain file could not place on chm13. gnomAD retention has not been tallied yet.
