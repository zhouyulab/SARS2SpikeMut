# SARS2SpikeMut

## This repository contains the following files:

- Snakefile: main workflow to pre-process SARS-CoV-2 GISAID spike protein sequences.
- scripts: analysis scripts for this project.
- data: GISAID sequencing data, reference SARS-CoV-2 Spike protein sequence and experimental data for this project.

## How to use:

1. Download SARS-CoV-2 Spike protein sequence from GISAID database (data/prot/spikeprot0705/README.txt).

2. Install the environment:
  * python packages:
    - pysam (0.15.4)
    - snakemake (5.4.4)
  * R packages:
    - readxl (1.3.1)
    - readr (1.3.1)
    - ggplot2 (3.3.0)
    - Biostrings (2.54.0)
    - BiocParallel (1.20.1)
    - ComplexHeatmap (2.2.0)
    - plotrix (3.8-2)

3. Use `snakemake` program to run the `Snakefile` file to perform the sequence preprocessing steps.

4. Use `R` program to run the `scripts/plot_variants.R` to plot the results.

## Developer

- Dehe Wang

