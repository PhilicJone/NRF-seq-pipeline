# NRF-seq-pipeline
NRF-seq (Nuclear RNase A Footprinting Sequencing) analysis pipeline for identifying trans-nuclear-envelope RNAs (tneRNAs). Includes bin-wise Pearson correlation, rank-difference scoring, candidate selection, and FDR estimation. Reproduces the discovery of KLNC as a prototypical tneRNA.
1. Quality Control and Preprocessing
Raw paired‑end FASTQ files are inspected with FastQC.
Adapter trimming and removal of low‑quality bases (Phred < 30) are performed with Trim Galore.
A second FastQC run confirms clean reads.
2. Alignment
Clean reads are aligned to the GENCODE GRCh38.v49​ reference genome using STAR​ (two‑pass mode for variant calling).
Only uniquely mapped reads are retained for downstream analysis.
4. Quantification
Gene‑level raw counts: featureCounts (exon‑based, strand‑specific).
Transcript‑level expression: StringTie (TPM).
Coverage tracks: bamCoverage (deepTools, bin size = 50, CPM normalization) produces BigWig files for visualization in IGV.

5. Candidate tneRNA Screening – Bin‑wise Pearson Correlation Coefficient (PCC) Screening
a. Transcripts (rather than genes) are split into 20 bins. Transcripts with fewer than three bins are excluded.
b. Bin‑wise read counts are extracted using computeMatrix (deepTools).
c. Intra‑group PCC​ is calculated between replicates (untreated nuclear vs. untreated nuclear; RNase A‑treated vs. RNase A‑treated). Only transcripts with intra‑group PCC > 0.8 are retained.
d. Inter‑group PCC​ is calculated between RNase A‑treated and untreated nuclear samples. Candidates must have inter‑group PCC lower than the corresponding intra‑group PCC.
e. Digested proportion: Within each transcript, bins are ranked by read count. A bin is called “RNase A‑sensitive” if its rank drops ≥30% after RNase A treatment. The digested proportion = sensitive bins / total bins.
f. Final ranking integrates a high rank in the background‑corrected bin‑wise PCC analysis.
