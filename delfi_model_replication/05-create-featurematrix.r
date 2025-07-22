# ---- Description ----

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
  "zfile", "z", 1, "character",  # --zfile or -z for file with z scores
  "covfile", "c", 1, "character",  # --covfile or -c for file with coverages
  "outfile", "o", 1, "character"  # --outfile or -o for output file
), byrow = TRUE, ncol = 4)

# Parse the command line arguments
opt <- getopt(spec)

# ---- Validate input ----
if (is.null(opt$zfile) || is.null(opt$covfile) || is.null(opt$outfile)) {
  cat(getopt(spec, usage = TRUE))
  quit(status = 1)
}

zfile <- opt$zfile
covfile <- opt$covfile
outfile <- opt$outfile

cat("Z file:", zfile, "\n")
cat("Cov file:", covfile, "\n")
cat("Outfile:", outfile, "\n")

# ---- Load files --------
zscores<- read_csv(zfile)
cov <- read_csv(covfile)

# --- Merge files on ID ---
setDT(zscores)
setDT(cov)

setkey(cov, id)
setkey(zscores, id)

full_features <- zscores[cov] # right join: all rows in 'cov' will be kept

# Define expected types for clinical columns
clinical_types <- list(
  clinical_nlratio = NA_real_,
  clinical_CRP = NA_real_,
  clinical_cfdna_conc = NA_real_,
  clinical_age = NA_real_,
  clinical_IL6 = NA_real_,
  clinical_YKL40 = NA_real_,
  clinical_CEA = NA_real_,
  clinical_bmi = NA_real_,
  clinical_packyears = NA_real_,
  clinical_smokingstatus = factor(NA, levels = c("never", "former", "current")),
  clinical_COPD = as.integer(NA),
  multinucratio = NA_real_
)

# Add columns with correct types
for (col in names(clinical_types)) {
  full_features[[col]] <- clinical_types[[col]]
}

# ---- Write output ----
write_csv(full_features,outfile)

