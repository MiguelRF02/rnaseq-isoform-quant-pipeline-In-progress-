# RNA-Seq Isoform Quantification Pipeline

**Status: In progress**

## Purpose

A reproducible RNA-seq pipeline, orchestrated with **Snakemake**, that goes from raw reads to transcript/isoform-level quantification using **Salmon**. The goal is to build a foundation for isoform-level analysis (as opposed to gene-level only) — connecting quality control, alignment-free quantification, and (eventually) differential transcript usage in a single reproducible workflow.

This project is being developed to build practical experience with workflow management systems and isoform-level RNA-seq analysis, complementing prior work on gene-expression-based classification (see other repos in this profile).

## Planned data source

- **[rnaseqDTU](https://github.com/mikelove/rnaseqDTU)** (Love, Soneson & Patro) — a Bioconductor workflow for differential transcript usage (DTU) following Salmon quantification. Uses small simulated RNA-seq reads.
