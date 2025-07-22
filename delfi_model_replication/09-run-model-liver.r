# ---- Description ----

# ---- Install and load packages ----
required_packages <- c("getopt", "recipes", "caret", "tidyverse", "here", "gbm")
# Install CRAN packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  suppressMessages(library(pkg, character.only = TRUE))
}

spec <- matrix(c(
  "featurefile", "f", 1, "character",  # --featurefile or -f for file with z scores
  "model", "m", 1, "character",  # --model or -m for file with coverages
  "outfile", "o", 1, "character"  # --outfile or -o for output file
), byrow = TRUE, ncol = 4)

# Parse the command line arguments
opt <- getopt(spec)

# ---- Validate input ----
if (is.null(opt$featurefile) || is.null(opt$model) || is.null(opt$outfile)) {
  cat(getopt(spec, usage = TRUE))
  quit(status = 1)
}

featurefile <- opt$featurefile
model_path <- opt$model
outfile <- opt$outfile

cat("Feature file:", featurefile, "\n")
cat("Model:", model_path, "\n")
cat("Outfile:", outfile, "\n")

# ---- Load features ----
features <- read_csv(featurefile)

# ---- Load model ----
model <- readRDS(model_path)

newfeatures <- features %>% select(-any_of(c(names(features)[startsWith(names(features), "cov_")], names(features)[startsWith(names(features), "clinical_")],"multinucratio")))
modelpreds <- predict(model, newdata=newfeatures, type="prob")
modelpreds <- modelpreds %>% mutate(id = newfeatures$id)

# ---- Save the results ----
write_csv(modelpreds, outfile)