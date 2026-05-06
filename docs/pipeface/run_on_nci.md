# Run pipeface on NCI

- [Run pipeface on NCI](#run-pipeface-on-nci)
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
    storage = 'gdata/if89+gdata/xy86+scratch/kr68+gdata/kr68'
```

> **_Note:_** Don't remove access to if89 gdata (`gdata/if89`). This is required to access software installs used in the pipeline. `gdata/xy86` is only required if running variant annotation and can be omitted from `storage` if not.

## 2. Start persistent session (optional)

Pipeface can be run in a screen within a [persistent session](https://opus.nci.org.au/spaces/Help/pages/241927941/Persistent+Sessions...).

## 3. Get pipeline dependencies

You may use the centrally installed nextflow environmental module available on NCI to access the nextflow dependency. Eg:

```bash
module load nextflow/25.10.3
```

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
