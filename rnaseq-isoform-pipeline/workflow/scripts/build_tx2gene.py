"""
Build tx2gene.tsv (transcript_id -> gene_id) by parsing a GTF file.
Invoked by Snakemake's `script:` directive (rule `build_tx2gene`).
"""
import re

gtf_path = snakemake.input["gtf"]      # noqa: F821 (injected by Snakemake)
out_path = snakemake.output["tx2gene"]  # noqa: F821

gene_re = re.compile(r'gene_id "([^"]+)"')
tx_re = re.compile(r'transcript_id "([^"]+)"')

seen = set()
n_written = 0
with open(gtf_path) as fin, open(out_path, "w") as fout:
    fout.write("transcript_id\tgene_id\n")
    for line in fin:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "transcript":
            continue
        attrs = fields[8]
        gm, tm = gene_re.search(attrs), tx_re.search(attrs)
        if gm and tm and tm.group(1) not in seen:
            seen.add(tm.group(1))
            fout.write(f"{tm.group(1)}\t{gm.group(1)}\n")
            n_written += 1

print(f"[build_tx2gene] wrote {n_written} transcript->gene mappings")
