"""
Generate a small synthetic isoform-quantification toy dataset.

Invoked by Snakemake's `script:` directive (see rule `simulate_toy_data`
in the Snakefile), so a `snakemake` object with .config/.output/.log is
injected automatically -- this is not meant to be run standalone.

Produces:
  - a tiny transcriptome FASTA: N genes, each with 2 isoforms related
    by a simple exon-skipping event
        isoform 1 = exonA + exonB + exonC   (full length)
        isoform 2 = exonA + exonC           (exon B skipped)
  - a transcript-to-gene mapping table (tx2gene.tsv)
  - a samples.csv (sample_id, condition)
  - gzipped paired-end FASTQ reads per sample, where isoform usage
    proportions differ by condition -> a toy differential transcript
    usage (DTU) signal for later analysis

Deliberately dependency-light (Python standard library only) so it
runs identically on any machine with no extra installs.
"""

import gzip
import os
import random
import sys

config = snakemake.config          # noqa: F821 (injected by Snakemake)
out = snakemake.output             # noqa: F821

toy_cfg = config["toy"]
seed = config.get("seed", 42)
random.seed(seed)

N_GENES = toy_cfg["n_genes"]
N_DTU_GENES = toy_cfg.get("n_dtu_genes", N_GENES)  # rest are "null" genes
EXON_LEN = toy_cfg["exon_length"]
READ_LEN = toy_cfg["read_length"]
READS_PER_SAMPLE = toy_cfg["reads_per_sample"]
SAMPLES = toy_cfg["samples"]
USAGE = toy_cfg["condition_isoform_usage"]

BASES = "ACGT"
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def random_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def mutate(seq, error_rate=0.01):
    """Inject a small amount of per-base sequencing error."""
    seq = list(seq)
    for i, base in enumerate(seq):
        if random.random() < error_rate:
            seq[i] = random.choice([b for b in BASES if b != base])
    return "".join(seq)


def simulate_fragment(seq, read_len):
    """Return (read1, read2) sequences for one simulated paired-end fragment."""
    frag_len = min(len(seq), read_len * 2 + random.randint(0, 50))
    start = 0 if len(seq) <= frag_len else random.randint(0, len(seq) - frag_len)
    frag = seq[start:start + frag_len]
    r1 = frag[:read_len]
    r2 = revcomp(frag[-read_len:])
    return mutate(r1), mutate(r2)


# ------------------------------------------------------------------
# 1. Build the tiny transcriptome
# ------------------------------------------------------------------
transcripts = {}   # transcript_id -> sequence
tx2gene = {}       # transcript_id -> gene_id

for g in range(1, N_GENES + 1):
    gene_id = f"gene{g}"
    exon_a, exon_b, exon_c = (random_seq(EXON_LEN) for _ in range(3))

    iso1, iso2 = f"{gene_id}_iso1", f"{gene_id}_iso2"
    transcripts[iso1] = exon_a + exon_b + exon_c
    transcripts[iso2] = exon_a + exon_c
    tx2gene[iso1] = gene_id
    tx2gene[iso2] = gene_id

os.makedirs(os.path.dirname(out.fasta), exist_ok=True)
with open(out.fasta, "w") as fh:
    for tx_id, seq in transcripts.items():
        fh.write(f">{tx_id}\n{seq}\n")

with open(out.tx2gene, "w") as fh:
    fh.write("transcript_id\tgene_id\n")
    for tx_id, gene_id in tx2gene.items():
        fh.write(f"{tx_id}\t{gene_id}\n")

# Ground truth: which genes were simulated with a real DTU signal vs
# which are "null" (constant 50/50 usage in both conditions). Useful
# to sanity-check the DRIMSeq/stageR output against.
with open(out.ground_truth, "w") as fh:
    fh.write("gene_id\ttrue_dtu\n")
    for g in range(1, N_GENES + 1):
        fh.write(f"gene{g}\t{'TRUE' if g <= N_DTU_GENES else 'FALSE'}\n")

# ------------------------------------------------------------------
# 2. samples.csv
# ------------------------------------------------------------------
os.makedirs(os.path.dirname(out.samples_csv), exist_ok=True)
with open(out.samples_csv, "w") as fh:
    fh.write("sample_id,condition\n")
    for s in SAMPLES:
        fh.write(f"{s['id']},{s['condition']}\n")

# ------------------------------------------------------------------
# 3. Simulate paired-end reads per sample
# ------------------------------------------------------------------
reads_by_sample = {s["id"]: {1: [], 2: []} for s in SAMPLES}

for s in SAMPLES:
    sample_id, condition = s["id"], s["condition"]
    dtu_props = USAGE[condition]  # e.g. [0.8, 0.2] for A, [0.2, 0.8] for B
    null_props = [0.5, 0.5]       # constant across conditions -> no DTU signal

    for g in range(1, N_GENES + 1):
        gene_id = f"gene{g}"
        iso1, iso2 = f"{gene_id}_iso1", f"{gene_id}_iso2"
        props = dtu_props if g <= N_DTU_GENES else null_props

        gene_reads = READS_PER_SAMPLE // N_GENES
        n_iso1 = int(gene_reads * props[0])
        n_iso2 = gene_reads - n_iso1

        for tx_id, n_reads in [(iso1, n_iso1), (iso2, n_iso2)]:
            seq = transcripts[tx_id]
            for i in range(n_reads):
                r1, r2 = simulate_fragment(seq, READ_LEN)
                read_name = f"{sample_id}:{tx_id}:{i}"
                reads_by_sample[sample_id][1].append((read_name, r1))
                reads_by_sample[sample_id][2].append((read_name, r2))

# Map each declared output path back to (sample_id, read_mate)
reads_out = {}
for path in out.reads:
    fname = os.path.basename(path)
    sample_id, mate = fname.replace(".fastq.gz", "").rsplit("_R", 1)
    reads_out[(sample_id, int(mate))] = path

os.makedirs(os.path.dirname(out.reads[0]), exist_ok=True)
for (sample_id, mate), path in reads_out.items():
    with gzip.open(path, "wt") as fh:
        for name, seq in reads_by_sample[sample_id][mate]:
            qual = "I" * len(seq)
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")

print(
    f"[simulate_toy_data] {len(SAMPLES)} samples, {N_GENES} genes x 2 isoforms "
    f"({N_DTU_GENES} with true DTU signal, {N_GENES - N_DTU_GENES} null), "
    f"seed={seed}.",
    file=sys.stderr,
)
