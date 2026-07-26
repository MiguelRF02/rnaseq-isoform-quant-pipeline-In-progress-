# RNA-Seq Isoform Quantification Pipeline

**Status: In progress**

## Purpose

A reproducible RNA-seq pipeline, orchestrated with **Snakemake**, that goes
from raw reads to transcript/isoform-level quantification using **Salmon**.
The goal is to build a foundation for isoform-level analysis (as opposed to
gene-level only) — connecting quality control, alignment-free quantification,
and (eventually) differential transcript usage in a single reproducible
workflow.

This project is being developed to build practical experience with workflow
management systems and isoform-level RNA-seq analysis, complementing prior
work on gene-expression-based classification (see other repos in this
profile).

## Pipeline

```
                 ┌───────────────────┐
   data_source → │ simulate (toy)    │
      "toy"      │   or               │
      "real"     │ download (real)   │
                 └─────────┬─────────┘
                           │
                  ┌────────┴────────┐
                  ▼                 ▼
              FastQC          Salmon index
                  │                 │
                  │                 ▼
                  │           Salmon quant (per sample)
                  │                 │
                  └────────┬────────┘
                           ▼
                       MultiQC report

  (por separado, usando los quant.sf ya generados)
  tximport → DRIMSeq → stageR → dtu_results.tsv + isoform_proportions.pdf
```

## Project structure

```
├── Snakefile              # pipeline rules
├── config.yaml            # toy/real switch + all parameters
├── environment.yml        # base env: just Snakemake
├── envs/                  # per-rule conda envs
│   ├── simulate.yaml
│   ├── fastqc.yaml
│   └── salmon.yaml
├── workflow/
│   └── scripts/
│       └── simulate_toy_data.py
├── data/                  # generated/downloaded inputs (gitignored)
└── results/               # QC reports + Salmon quantifications (gitignored)
```

## Quickstart (toy dataset)

The default `config.yaml` has `data_source: "toy"`, which generates a small
synthetic dataset (a handful of genes, each with two isoforms related by
exon skipping, with isoform usage proportions that differ between two
conditions) instead of downloading anything. This is meant for testing the
pipeline logic quickly, end to end, before scaling up.

```bash
# 1. create the base environment (just needs Snakemake)
conda env create -f environment.yml
conda activate rnaseq-isoform-pipeline

# 2. run the whole pipeline; --use-conda lets each rule create/use its
#    own tool-specific environment (fastqc, salmon, etc.) automatically
snakemake --cores 4 --use-conda

# 3. see what would run without actually running it
snakemake -n --use-conda

# 4. visualize the DAG
snakemake --dag | dot -Tpng > dag.png
```

Outputs land in `results/qc/` (FastQC reports) and `results/quant/<sample>/quant.sf`
(one Salmon quantification file per sample, ready for `tximport`).

## Planned data source (for scaling beyond the toy dataset)

- **[rnaseqDTU](https://github.com/mikelove/rnaseqDTU)** (Love, Soneson &
  Patro) — a Bioconductor workflow for differential transcript usage (DTU)
  following Salmon quantification. The simulated FASTQ reads for this
  workflow are hosted on Zenodo (not in the GitHub repo itself), in three
  batches of eight samples each, simulated against the GENCODE v28 human
  transcriptome.

To switch from the toy dataset to this real dataset:

1. Resolve the exact Zenodo record(s) for the FASTQ batches and for the
   GENCODE v28 transcriptome FASTA/GTF (see the *Data availability* section
   of the [F1000Research paper](https://f1000research.com/articles/7-952)).
2. Fill in the `real:` section of `config.yaml` (`transcriptome_fasta_url`,
   `samples`, `reads_urls`).
3. Set `data_source: "real"` in `config.yaml`.
4. Re-run `snakemake --cores <n> --use-conda` — no changes to the
   `Snakefile` logic are needed, since `fastqc`, `salmon_index`, and
   `salmon_quant` are already data-source-agnostic.

Note: the real dataset is at realistic scale (24 samples, ~31–38 million
paired-end reads each), so budget disk space, bandwidth, and compute time
accordingly — this is not something to run casually on a laptop without
planning for it.

## Roadmap

- [x] Project scaffold + Snakemake skeleton
- [x] Toy dataset generation for fast local iteration
- [x] FastQC + Salmon index + Salmon quant rules
- [x] MultiQC aggregation step
- [x] Differential transcript usage on the toy dataset (tximport / DRIMSeq / stageR)
- [ ] Switch to real rnaseqDTU data (Zenodo + GENCODE v28)
- [ ] Downstream differential transcript usage (DRIMSeq / DEXSeq / stageR),
      following the rnaseqDTU workflow
