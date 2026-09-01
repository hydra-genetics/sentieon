# :snake: hydra-genetics/sentieon

#### sentieon tools

![CI](https://github.com/hydra-genetics/sentieon/actions/workflows/ci.yaml/badge.svg?branch=develop)

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/gpl-3.0.html)

## :speech_balloon: Introduction

The module consists of tools from Sentieon, including alignment, duplicate reads removal, indel realignment, recalibration table calculation and DNAScope and TNScope SNV/structural variants callers.

## :heavy_exclamation_mark: Dependencies

In order to use this module, the following dependencies are required:

[![hydra-genetics](https://img.shields.io/badge/hydragenetics-snakemake9%20branch-blue)](https://github.com/hydra-genetics/hydra-genetics/tree/migrate-to-snakemake9-python3.12)
[![pandas](https://img.shields.io/badge/pandas-1.3.1-blue)](https://pandas.pydata.org/)
[![python](https://img.shields.io/badge/python-3.12-blue)](https://www.python.org/)
[![snakemake](https://img.shields.io/badge/snakemake-9.0.0-blue)](https://snakemake.readthedocs.io/en/stable/)
[![apptainer](https://img.shields.io/badge/apptainer-1.4.5-blue)](https://apptainer.org/)

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
$ snakemake -s ../../workflow/Snakefile -j1 --software-deployment-method apptainer
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
