# :snake: hydra-genetics/sentieon

#### sentieon tools

![Lint](https://github.com/hydra-genetics/sentieon/actions/workflows/lint.yaml/badge.svg?branch=develop)
![Snakefmt](https://github.com/hydra-genetics/sentieon/actions/workflows/snakefmt.yaml/badge.svg?branch=develop)
![snakemake dry run](https://github.com/hydra-genetics/sentieon/actions/workflows/snakemake-dry-run.yaml/badge.svg?branch=develop)
![integration test](https://github.com/hydra-genetics/sentieon/actions/workflows/integration1.yaml/badge.svg?branch=develop)

![pycodestyle](https://github.com/hydra-genetics/sentieon/actions/workflows/pycodestyl.yaml/badge.svg?branch=develop)
![pytest](https://github.com/hydra-genetics/sentieon/actions/workflows/pytest.yaml/badge.svg?branch=develop)

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/gpl-3.0.html)

## :speech_balloon: Introduction

The module consists of tools from Sentieon, including alignment, duplicate reads removal, indel realignment, recalibration table calculation and DNAScope and TNScope SNV/structural variants callers.

## :heavy_exclamation_mark: Dependencies

In order to use this module, the following dependencies are required:

[![hydra-genetics](https://img.shields.io/badge/hydragenetics-v0.9.1-blue)](https://github.com/hydra-genetics/)
[![pandas](https://img.shields.io/badge/pandas-1.3.1-blue)](https://pandas.pydata.org/)
[![python](https://img.shields.io/badge/python-3.8-blue)
[![snakemake](https://img.shields.io/badge/snakemake-6.8.0-blue)](https://snakemake.readthedocs.io/en/stable/)
[![singularity](https://img.shields.io/badge/singularity-3.0.0-blue)](https://sylabs.io/docs/)

## :school_satchel: Preparations

### Sample data

Input data should be added to [`samples.tsv`](https://github.com/hydra-genetics/sentieon/blob/develop/config/samples.tsv)
and [`units.tsv`](https://github.com/hydra-genetics/sentieon/blob/develop/config/units.tsv).
The following information need to be added to these files:

| Column Id | Description |
| --- | --- |
| **`samples.tsv`** |
| sample | unique sample/patient id, one per row |
| **`units.tsv`** |
| sample | same sample/patient id as in `samples.tsv` |
| type | data type identifier (one letter), can be one of **T**umor, **N**ormal, **R**NA |
| platform | type of sequencing platform, e.g. `NovaSeq` |
| machine | specific machine id, e.g. NovaSeq instruments have `@Axxxxx` |
| flowcell | identifer of flowcell used |
| lane | flowcell lane number |
| barcode | sequence library barcode/index, connect forward and reverse indices by `+`, e.g. `ATGC+ATGC` |
| fastq1/2 | absolute path to forward and reverse reads |
| adapter | adapter sequences to be trimmed, separated by comma |

## :white_check_mark: Testing

The workflow repository contains a small test dataset `.tests/integration` which can be run like so:

```bash
$ cd .tests/integration
$ snakemake -s ../../Snakefile -j1 --use-singularity
```

## :rocket: Usage

To use this module in your workflow, follow the description in the
[snakemake docs](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).
Add the module to your `Snakefile` like so:

```bash
module sentieon:
    snakefile:
        github(
            "hydra-genetics/sentieon",
            path="workflow/Snakefile",
            tag="add_tools",
        )
    config:
        config


use rule * from sentieon as sentieon_*
```

### Output files

The following output files should be targeted via another rule:

| File | Description |
|---|---|
| `sentieon/realign/{sample}_{type}_REALIGNED.bam` | Aligned data, duplicates removed and indels realigned |
| `sentieon/qualcal/{sample}_{type}_RECAL_DATA.TABLE` | Recalibration table |
| `sentieon/tnscope/{sample}_TNscope_tn_ML.vcf` | Output SNVs/indels vcf for tumor with matched normal |
| `sentieon/dnascope/{sample}_{type}_DNAscope_modelfiltered.vcf` | Output SNVs/indels vcf for germline |



## :judge: Rule Graph
