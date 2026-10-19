#!/usr/bin/env Rscript

# Load required packages
library(data.table)
library(matrixStats)

# ------------------------------------------------------------------------------
# 1. Read the computeMatrix output matrix (decompressed file)
#    The file has a '#' comment line at top, then a header line, then data rows.
#    Using comment.char="#" skips the first line automatically.
# ------------------------------------------------------------------------------
mat_file <- "Output.mat"   # change to your actual file path if needed
mat <- fread(mat_file, comment.char = "#", header = TRUE, sep = "\t")
rownames(mat) <- mat[[1]]          # first column: transcript name
mat <- as.matrix(mat[, -1])         # keep only numeric columns

# ------------------------------------------------------------------------------
# 2. Define sample groups and number of bins
#    Expected order: noTreat1, noTreat2, RNase1, RNase2, each with 20 bins
# ------------------------------------------------------------------------------
nBins <- 20
sampleNames <- c("noTreat1", "noTreat2", "RNase1", "RNase2")
nSamples <- length(sampleNames)

# Check that the matrix dimensions match expectations
stopifnot(ncol(mat) == nSamples * nBins)

# Split matrix into per-sample bin matrices (rows = transcripts, cols = bins)
splitCols <- lapply(seq_len(nSamples), function(i) {
  cols <- ((i - 1) * nBins + 1):(i * nBins)
  mat[, cols]
})
names(splitCols) <- sampleNames

# ------------------------------------------------------------------------------
# 3. Compute intra-group Pearson correlation coefficients
#    (noTreat1 vs noTreat2; RNase1 vs RNase2)
# ------------------------------------------------------------------------------
intra_pcc_noTreat <- sapply(seq_len(nrow(mat)), function(r) {
  cor(splitCols[["noTreat1"]][r, ], splitCols[["noTreat2"]][r, ], method = "pearson")
})
intra_pcc_RNase <- sapply(seq_len(nrow(mat)), function(r) {
  cor(splitCols[["RNase1"]][r, ], splitCols[["RNase2"]][r, ], method = "pearson")
})

# ------------------------------------------------------------------------------
# 4. Compute inter-group PCC (between mean vectors of the two groups)
# ------------------------------------------------------------------------------
mean_noTreat <- (splitCols[["noTreat1"]] + splitCols[["noTreat2"]]) / 2
mean_RNase   <- (splitCols[["RNase1"]] + splitCols[["RNase2"]]) / 2

inter_pcc <- sapply(seq_len(nrow(mat)), function(r) {
  cor(mean_noTreat[r, ], mean_RNase[r, ], method = "pearson")
})

# ------------------------------------------------------------------------------
# 5. Apply filtering criteria
#    a) At least 3 non???NA bins (already guaranteed by skipZeros, but kept for safety)
#    b) Both intra???group PCC > 0.8
#    c) Inter???group PCC < minimum of the two intra???group PCC values
# ------------------------------------------------------------------------------
valid_bins <- rowSums(!is.na(mat)) >= 3
cond_intra <- (intra_pcc_noTreat > 0.8) & (intra_pcc_RNase > 0.8)
min_intra  <- pmin(intra_pcc_noTreat, intra_pcc_RNase)
cond_inter <- inter_pcc < min_intra

candidates <- valid_bins & cond_intra & cond_inter
cat("Number of transcripts passing basic filters:", sum(candidates), "\n")

# ------------------------------------------------------------------------------
# 6. Compute Digested Proportion (RNase???sensitive bin fraction)
#    For each candidate transcript:
#      - Rank bins by mean read count in untreated group
#      - If rank drops >= 30% after RNase treatment, mark bin as sensitive
#      - Digested proportion = sensitive bins / total bins
# ------------------------------------------------------------------------------
digest_prop <- rep(NA, nrow(mat))
for (r in which(candidates)) {
  before_rank <- rank(mean_noTreat[r, ], ties.method = "average")
  after_rank  <- rank(mean_RNase[r, ], ties.method = "average")
  rank_drop   <- (before_rank - after_rank) / before_rank
  sensitive   <- sum(rank_drop >= 0.3, na.rm = TRUE)
  digest_prop[r] <- sensitive / nBins
}

# ------------------------------------------------------------------------------
# 7. Build final result table
# ------------------------------------------------------------------------------
result <- data.frame(
  transcript           = rownames(mat)[candidates],
  intra_PCC_noTreat    = intra_pcc_noTreat[candidates],
  intra_PCC_RNase      = intra_pcc_RNase[candidates],
  inter_PCC            = inter_pcc[candidates],
  digested_proportion  = digest_prop[candidates],
  stringsAsFactors     = FALSE
)

# Sort by digested proportion (descending)
result <- result[order(result$digested_proportion, decreasing = TRUE), ]

# ------------------------------------------------------------------------------
# 8. Save results to file
# ------------------------------------------------------------------------------
out_file <- "tneRNA_candidates.txt"
write.table(result, file = out_file,
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Candidate list saved to", out_file, "\n")