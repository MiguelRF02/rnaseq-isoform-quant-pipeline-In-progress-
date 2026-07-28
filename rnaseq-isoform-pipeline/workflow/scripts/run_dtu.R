# ============================================================
# Differential Transcript Usage (DTU) analysis
# ============================================================
# Invoked by Snakemake's `script:` directive (rule `dtu_analysis`),
# so a `snakemake` S4 object is injected automatically with
# snakemake@input, @output, @log, @config -- not meant to be run
# standalone.
#
# Pipeline: tximport (read Salmon quant.sf) -> DRIMSeq (test isoform
# proportion differences per gene) -> stageR (two-stage testing to
# control the false discovery rate across gene- and transcript-level
# tests). Same overall approach as the rnaseqDTU workflow, adapted to
# our small toy dataset.

# Redirect all console output to the Snakemake log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")

suppressMessages({
  library(tximport)
  library(DRIMSeq)
  library(stageR)
  library(dplyr)
  library(readr)
})

quant_files <- unlist(snakemake@input[["quant"]])
tx2gene_path <- snakemake@input[["tx2gene"]]
samples_path <- snakemake@input[["samples_csv"]]
out_table <- snakemake@output[["table"]]
out_plot <- snakemake@output[["plot"]]
seed <- snakemake@config[["seed"]]

set.seed(seed)


# 1. Load sample metadata and transcript-to-gene mapping

# quant.sf paths look like .../quant/<sample_id>/quant.sf -- recover
# the sample_id from the parent directory name
sample_ids <- basename(dirname(quant_files))
names(quant_files) <- sample_ids

samples <- read_csv(samples_path, show_col_types = FALSE)
samples <- samples[match(sample_ids, samples$sample_id), ]
stopifnot(all(samples$sample_id == sample_ids))

tx2gene <- read_tsv(tx2gene_path, show_col_types = FALSE)

cat(sprintf(
  "Loaded %d samples, %d transcripts (tx2gene), conditions: %s\n",
  length(sample_ids), nrow(tx2gene),
  paste(unique(samples$condition), collapse = ", ")
))


# 2. tximport: transcript-level counts from Salmon's quant.sf

txi <- tximport(
  quant_files,
  type = "salmon",
  txOut = TRUE,                      # keep transcript-level (not gene-level)
  countsFromAbundance = "scaledTPM", # recommended by DRIMSeq for count-based testing
  dropInfReps = TRUE                 # we didn't run salmon with --numBootstraps,
                                      # so there are no inferential replicates to import
)

counts_df <- data.frame(
  feature_id = rownames(txi$counts),
  gene_id = tx2gene$gene_id[match(rownames(txi$counts), tx2gene$transcript_id)],
  txi$counts,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# 3. DRIMSeq: model isoform usage proportions per gene

pd <- data.frame(sample_id = samples$sample_id, condition = samples$condition)

d <- dmDSdata(counts = counts_df, samples = pd)


d <- dmFilter(
  d,
  min_samps_gene_expr = nrow(pd),   # gene expressed in every sample
  min_samps_feature_expr = 3,       # isoform expressed in >=3 samples
  min_gene_expr = 10,
  min_feature_expr = 5
)

cat(sprintf("%d genes / %d transcripts left after filtering.\n",
            length(unique(counts(d)$gene_id)), nrow(counts(d))))

design_full <- model.matrix(~condition, data = DRIMSeq::samples(d))

d <- dmPrecision(d, design = design_full)
d <- dmFit(d, design = design_full)
d <- dmTest(d, coef = colnames(design_full)[ncol(design_full)])

res_gene <- DRIMSeq::results(d)
res_tx <- DRIMSeq::results(d, level = "feature")


# 4. stageR: two-stage testing (screen genes, then confirm transcripts)

res_tx_clean <- res_tx[!is.na(res_tx$pvalue), ]

pScreen <- res_gene$pvalue
names(pScreen) <- res_gene$gene_id

pConfirmation <- matrix(res_tx_clean$pvalue, ncol = 1)
rownames(pConfirmation) <- res_tx_clean$feature_id

tx2gene_stager <- data.frame(
  txID = res_tx_clean$feature_id,
  geneID = res_tx_clean$gene_id
)

stageRObj <- stageRTx(
  pScreen = pScreen,
  pConfirmation = pConfirmation,
  pScreenAdjusted = FALSE,
  tx2gene = tx2gene_stager
)
stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu", alpha = 0.05)
stager_res <- getAdjustedPValues(stageRObj, order = TRUE, onlySignificantGenes = FALSE)


# 5. Write results table

merged <- merge(res_tx_clean, stager_res, by.x = "feature_id", by.y = "txID")
merged <- merged[order(merged$transcript), ]
write_tsv(merged, out_table)

cat(sprintf("Wrote DTU results table: %s (%d rows)\n", out_table, nrow(merged)))


# 6. Proportions plot per gene (visualize the isoform-usage swap)

gene_ids <- unique(counts(d)$gene_id)
pdf(out_plot, width = 7, height = 4)
for (g in gene_ids) {
  print(plotProportions(d, gene_id = g, group_variable = "condition"))
}
dev.off()

cat(sprintf("Wrote proportions plot: %s (%d genes)\n", out_plot, length(gene_ids)))

sink(type = "message")
sink()
