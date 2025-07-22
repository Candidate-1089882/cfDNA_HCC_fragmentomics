# ---- Description ----
# This script processes fragmentomics data from multiple CSV files in a specified directory.
# It computes standardized fragment size ratios and coverage metrics per bin and per sample,
# reshapes them into a wide format, and writes the result to an output file.

# ---- Install and load packages ----
required_packages <- c("getopt", "data.table", "dplyr", "tidyverse", "readr")

# Install CRAN packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  suppressMessages(library(pkg, character.only = TRUE))
}

# ---- Define needed functions -----

# ---- Define command-line options ----
spec <- matrix(c(
  "datadir", "d", 1, "character",  # --datadir or -d for directory with the data
  "outfile", "o", 1, "character"  # --outfile or -o for output file
), byrow = TRUE, ncol = 4)

# Parse the command line arguments
opt <- getopt(spec)

# ---- Validate input ----
if ( is.null(opt$datadir) || is.null(opt$outfile)) {
  cat(getopt(spec, usage = TRUE))
  quit(status = 1)
}

datadir <- opt$datadir
outfile <- opt$outfile

cat("Datadir", datadir, "\n")
cat("Outfile:", outfile, "\n")

# ---- Load sample files --------
filenames<-list.files(datadir,pattern=".csv",full.names=TRUE)
x <- rbindlist(lapply(filenames, fread))
x_f <- setDT(x %>% filter(chr != "chrX"))

# ---- Compute base features ----

# Set key by sample ID for faster operations
setkey(x_f, id)

# Total raw coverage
x_f[, cov := short + long]

# ---- Fragment size ratio (corrected) ----
x_f[, ratio.cor := short.cor / long.cor]

# Standardize ratio per sample
x_f[, ratio.scaled := scale(ratio.cor), by = id]
x_f[, ratiovar := factor(paste0("ratio_", bin), paste0("ratio_", 1:.N)), by = id]

# Reshape to wide format: one row per sample, one column per bin
features.ratios <- dcast(x_f, id ~ ratiovar, value.var = "ratio.scaled")

# ---- Coverage features ----

# Total corrected coverage
x_f[, cov.cor := short.cor + long.cor]

# Standardize coverage per sample
x_f[, cov.scaled := scale(cov.cor), by = id]
x_f[, covvar := factor(paste0("cov_", bin), paste0("cov_", 1:.N)), by = id]

# Reshape to wide format: one row per sample, one column per bin
features.covs <- dcast(x_f, id ~ covvar, value.var = "cov.scaled")

# ---- Merge coverage and ratio features ----

setkey(features.ratios, id)
setkey(features.covs, id)
features.full <- features.covs[features.ratios]

# ---- Write output ----
write_csv(features.full,outfile)