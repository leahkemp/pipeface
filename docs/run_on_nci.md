# Run on NCI

- [Run on NCI](#run-on-nci)
  - [Assumptions](#assumptions)
  - [1. Modify nextflow\_pipeface.config](#1-modify-nextflow_pipefaceconfig)
  - [2. Start persistent session (optional)](#2-start-persistent-session-optional)
  - [3. Get pipeline dependencies](#3-get-pipeline-dependencies)
  - [4. Run pipeface](#4-run-pipeface)
  - [Information](#information)

## Assumptions

- Running pipeline on Australia's [National Computational Infrastructure (NCI)](https://nci.org.au/)
- Access to if89 project (to access software installs used by pipeface)
- Access to xy86 project (to access variant databases used by pipeface, only required if running variant annotation)

## 1. Modify nextflow_pipeface.config

Modify the NCI project to which to charge the analysis. Eg:

```txt
    project = 'kr68'
```

Modify access to project specific directories. Eg:

```txt
    storage = 'gdata/if89+gdata/xy86+scratch/kr68+gdata/kr68+gdata/ox63'
```

> **_Note:_** Don't remove access to if89 gdata (`gdata/if89`) and xy86 gdata (`gdata/xy86`). These are required to access software installs and variant annotation databases used in the pipeline

## 2. Start persistent session (optional)

Pipeface can be run within a [persistent session](https://opus.nci.org.au/spaces/Help/pages/241927941/Persistent+Sessions...)

## 3. Get pipeline dependencies

You may use the centrally installed nextflow environmental module available on NCI to access the nextflow dependency. Eg:

```bash
module load nextflow/24.04.5
```

## 4. Run pipeface

```bash
nextflow run pipeface.nf -params-file ./config/parameters_pipeface.json -config ./config/nextflow_pipeface.config
```

## Information

Please keep in mind that some datasets will require modifications to the default resources (particularly memory, disk usage, walltime). For example WGS data with greater than typical (~30x) sequencing depth.

