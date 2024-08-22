# Software versions

*Assuming the default [nextflow_pipeface.config](../config/nextflow_pipeface.config) file is used*

## For entire pipeline

### Main tools:

- Samtools: 1.19
- Minimap2: 2.28-r1209
- Clair3: 1.0.9
- Clara Parabricks: 4.2.1-1.beta4 (equivalent to - DeepVariant 1.5.0)
- WhatsHap: 2.3
- Sniffles2: 2.3.3
- cuteSV: 1.0.13

### Other tools

- GNU coreutils: 8.30
- HTSlib: 1.16
- GNU Awk: 4.2.1

## Per process

- scrape_settings: GNU coreutils 8.30
- scrape_bam_header: Samtools 1.19
- merge_runs: GNU coreutils 8.30, HTSlib 1.16, Samtools 1.19
- minimap2: Samtools 1.19, Minimap2 2.28-r1209
- depth: Samtools 1.19, GNU coreutils 8.30, GNU Awk 4.2.1
- clair3: Clair3 1.0.9, GNU coreutils 8.30
- deepvariant: Clara Parabricks 4.2.1-1.beta4, HTSlib 1.16
- whatshap_phase_clair3, whatshap_phase_dv, whatshap_haplotag: WhatsHap 2.3, Samtools 1.19, HTSlib 1.16
- sniffles: Sniffles2 2.3.3
- cutesv: cuteSV 1.0.13, HTSlib 1.16
- publish_settings, publish_bam_header, publish_minimap2, publish_clair3, publish_deepvariant, publish_whatshap_phase_clair3, publish_whatshap_phase_dv, publish_whatshap_haplotag, publish_sniffles, publish_cutesv: none
