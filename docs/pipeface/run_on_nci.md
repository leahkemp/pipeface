# Run pipeface on NCI

- [Run pipeface on NCI](#run-pipeface-on-nci)
  - [Assumptions](#assumptions)
  - [1. Modify nextflow\_pipeface\_nci.config](#1-modify-nextflow_pipeface_nciconfig)
  - [2. Start persistent session (optional)](#2-start-persistent-session-optional)
  - [3. Get pipeline dependencies](#3-get-pipeline-dependencies)
  - [4. Run pipeface](#4-run-pipeface)
  - [Information](#information)

## Assumptions

- Running pipeline on Australia's [National Computational Infrastructure (NCI)](https://nci.org.au/)
- Access to if89 project (to access the VEP cache used by pipeface, only required if running variant annotation)
- Access to xy86 project (to access variant databases used by pipeface, only required if running variant annotation)

## 1. Modify nextflow_pipeface_nci.config

`nextflow_pipeface_nci.config` layers the NCI settings on top of the generic `nextflow_pipeface.config` (software containers, process resources).

Modify the NCI project to which to charge the analysis. Eg:

```txt
    project = 'kr68'
```

Modify access to project specific directories. Eg:

```txt
    storage = 'gdata/if89+gdata/xy86+scratch/kr68+gdata/kr68'
```

Modify the directory the software containers are pulled to (the first run pulls them, later runs re-use them). Eg:

```txt
    cacheDir = '/g/data/kr68/install/singularity_cache'
```

> [!NOTE]
> Add the project holding the container cache directory to `storage` (eg. `gdata/kr68` above). `gdata/if89` and `gdata/xy86` are only required if running variant annotation and can be omitted from `storage` if not.

## 2. Start persistent session (optional)

Pipeface can be run in a screen within a [persistent session](https://opus.nci.org.au/spaces/Help/pages/241927941/Persistent+Sessions...).

## 3. Get pipeline dependencies

You can use the nextflow and singularity environmental modules available on NCI. Eg:

```bash
module load nextflow/25.10.3 singularity
```

## 4. Run pipeface

Run the pipeline. Eg:

```bash
nextflow run pipeface.nf -params-file ./config/parameters_pipeface.json -config ./config/nextflow_pipeface_nci.config
```

If you need to resume a pipeline run, use the `-resume` flag. Eg:

```bash
nextflow run pipeface.nf -resume -params-file ./config/parameters_pipeface.json -config ./config/nextflow_pipeface_nci.config
```

## Information

Please keep in mind that some datasets will require modifications to the default resources (particularly memory, disk usage, walltime). For example WGS data with greater than typical (~30x) sequencing depth.
