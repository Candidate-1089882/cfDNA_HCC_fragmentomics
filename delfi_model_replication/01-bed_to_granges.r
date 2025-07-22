# ---- Description ----
# This script processes a BED file, filters fragments based on size, 
# computes GC content using a reference cytosine file, and saves the 
# results as a GRanges object in an RDS file.

# ---- Install and load packages ----
required_packages <- c("getopt", "data.table", "ggplot2")
bioconductor_packages <- c("GenomicRanges", "Biostrings")

# Install CRAN packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  suppressMessages(library(pkg, character.only = TRUE))
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in bioconductor_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  suppressMessages(library(pkg, character.only = TRUE))
}

# ---- Define command-line options ----
spec <- matrix(c(
  "id",     "i", 1, "character",  # --id or -i for BED file
  "outdir", "o", 1, "character",  # --outdir or -o for output dir
  "refdir", "r", 1, "character",  # --refdir or -r for reference directory
  "figuredir", "f", 1, "character"  # --figuredir or -f for figure directory
), byrow = TRUE, ncol = 4)

# Parse the command line arguments
opt <- getopt(spec)


# ---- Validate input ----
if (is.null(opt$id) || is.null(opt$outdir) || is.null(opt$refdir)) {
  cat(getopt(spec, usage = TRUE))
  quit(status = 1)
}

bedfile <- opt$id
outdir <- opt$outdir
refdir <- opt$refdir
figuredir <- opt$figuredir 
sample_name <- gsub(".bed$", "", basename(bedfile))

cat("Sample:", sample_name, "\n")
cat("BED file:", bedfile, "\n")
cat("Reference directory:", refdir, "\n")
cat("Figure directory:", figuredir, "\n")

# ---- Load and process BED file ----
bed <- fread(bedfile)
setnames(bed, c("chr", "start", "end", "mapq"))

# bed files are 0 indexed at start position and 1 indexed at end.
# For consistency in R, add 1 to start position but leave end alone.
# Filter on mapq > 30 to retain only high quality alignments.
# Filter outliers. Very large fragments are likely alignment artifacts.
bed <- bed[, start := start + 1][mapq >= 30]
bed <- bed[(end - start) < 1000] 
frags <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)

# ---- Load cytosine reference and compute GC content ----
cytosines <- readRDS(file.path(refdir, "cytosine_ref_hg19.rds"))

frags$gc_count <- countOverlaps(frags, cytosines, ignore.strand = TRUE)
frags$w <- width(frags)
frags$gc <- frags$gc_count / frags$w

# ---- Save results ----
saveRDS(frags, file.path(outdir, paste0(sample_name, ".rds")))
cat("Fragment file saved to", outdir, "\n")

# ---- Sanity checks -----
# Get a summary of the GRanges object
cat("Summary of GRanges object:\n")
summary(frags) #general summary
cat("\nSummary of GC content distribution:\n")
summary(frags$gc)  # GC content distribution
cat("\nSummary of Fragment Width distribution:\n")
summary(frags$w)    # Fragment width distribution

# ---- Convert GRanges metadata to data.frame ----
frag_df <- as.data.frame(mcols(frags))

# include seqnames for chromosome plots
frag_df$seqnames <- as.character(seqnames(frags))

# ---- Downsample fragment data ----
set.seed(42)  # For reproducibility
n_sample <- 100000  # Sample 100,000 fragments

if (nrow(frag_df) > n_sample) {
  frag_df <- frag_df[sample(seq_len(nrow(frag_df)), n_sample), ]
  cat("Subsetted fragment dataframe to", n_sample, "rows for plotting.\n")
}

# ---- Create plots ----

# 1. Histogram of Fragment Widths
p1 <- ggplot(frag_df, aes(x = w)) +
  geom_histogram(binwidth = 10, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Fragment Width Distribution", x = "Fragment Width", y = "Count")
ggsave(filename = file.path(figuredir, paste0(sample_name, "_frag_width_histogram.png")), plot = p1, width = 8, height = 6)

# 2. Histogram of GC Content
p2 <- ggplot(frag_df, aes(x = gc)) +
  geom_histogram(binwidth = 0.02, fill = "lightgreen", color = "black") +
  theme_minimal() +
  labs(title = "GC Content Distribution", x = "GC Content", y = "Count")
ggsave(filename = file.path(figuredir, paste0(sample_name, "_gc_content_histogram.png")), plot = p2, width = 8, height = 6)

# 3. Boxplot of GC Content vs Fragment Width (after binning width)
frag_df$w_bin <- cut(frag_df$w, breaks = seq(0, 1000, by = 50)) # 50 bp bins

p3 <- ggplot(frag_df, aes(x = w_bin, y = gc)) +
  geom_boxplot(fill = "lightcoral", outlier.size = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6)) +
  labs(title = "GC Content vs Fragment Width", x = "Fragment Width Bin (bp)", y = "GC Content")
ggsave(filename = file.path(figuredir, paste0(sample_name, "_gc_vs_width_boxplot.png")), plot = p3, width = 10, height = 6)

# 4. Density plot of Fragment Widths
p4 <- ggplot(frag_df, aes(x = w)) +
  geom_density(fill = "blue", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Density Plot of Fragment Widths", x = "Fragment Width", y = "Density")
ggsave(filename = file.path(figuredir, paste0(sample_name, "_frag_width_density_plot.png")), plot = p4, width = 8, height = 6)

# 5. Density plot of GC Content
p5 <- ggplot(frag_df, aes(x = gc)) +
  geom_density(fill = "red", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Density Plot of GC Content", x = "GC Content", y = "Density")
ggsave(filename = file.path(figuredir, paste0(sample_name, "_gc_content_density_plot.png")), plot = p5, width = 8, height = 6)

# 6. Fragment Counts per Chromosome (sanity check)
p6 <- ggplot(frag_df, aes(x = seqnames)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Fragment Counts per Chromosome", x = "Chromosome", y = "Fragment Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = file.path(figuredir, paste0(sample_name, "_frag_locations_barplot.png")), plot = p6, width = 10, height = 6)

cat("Sanity check plots saved to", figuredir, "\n")