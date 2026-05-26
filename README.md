# Pipeface/popface

## Overview

Nextflow pipelines to process long read [ONT](https://nanoporetech.com/) and/or [pacbio](https://www.pacb.com/) HiFi data.

Pipeface takes unaligned BAMs, FASTQs or aligned BAMs and generates variant calls for singletons, duos or trios. Popface takes the outputs of pipeface to generate quad to population variant calls for up to 1000 individuals.

See more information about [pipeface](./docs/pipeface/info.md) and [popface](./docs/popface/info.md).

<p align="center">
    <img src="./images/pipeface.png">

## Run it!

See a walkthrough for [setting up pipeface](./docs/pipeface/setup.md) and [running it on NCI](./docs/pipeface/run_on_nci.md) or [another HPC](./docs/pipeface/run_on_other_hpc.md). Similarly here is a walkthrough for [setting up popface](./docs/popface/setup.md) and [running it on NCI](./docs/popface/run_on_nci.md) or [another HPC](./docs/popface/run_on_other_hpc.md).

## Credit

This is a highly collaborative project, with many contributions from the [Deveson Lab](https://www.garvan.org.au/research/labs-groups/deveson-lab). Notably, Dr Andre Reis and Dr Ira Deveson are closely involved in the development of this pipeline. Optimisations involving DeepVariant and DeepTrio have been contributed by Dr Kisaru Liyanage and Dr Matthew Downton from the [National Computational Infrastructure](https://nci.org.au), with support from Australian BioCommons as part of the Workflow Commons project. Haploid-aware mode has been contributed by Dr Hardip Patel & Kirat Alreja from the [National Centre for Indigenous Genomics](https://ncig.anu.edu.au). The installation and hosting of software used in this pipeline has and continues to be supported by the [Australian BioCommons Tools and Workflows project (if89)](https://australianbiocommons.github.io/ables/if89/).
