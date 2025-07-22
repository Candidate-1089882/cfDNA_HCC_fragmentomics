# ---- Description ----
# This script processes genomic bin-level coverage data to compute normalized coverage and chromosome arm-level z-scores
# relative to a reference dataset. It performs GC correction, normalization for filtered bases, and outputs a wide-format
# z-score matrix for each sample that can be used in downstream ML analysis.

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
countAndNormalize <- function(bins, measure) {
    ### normalize for filtered bases
    ### weight is to get count per expected width of basepairs
    ## Maybe 
    bins[,weight:=(end - start + 1)/(end - start + 1 - filtered.bases)]
    bins[,counts.mult:=log2((get(measure) * weight + 1))]
    bins[,loess.pred:=gcCorrectLoess(counts.mult, gc),by=id]
    bins[,adjusted:=counts.mult - loess.pred]
    bins[]
}

getArmMeans <- function(bins) {
        bins[,adjusted.cent:=adjusted-median(adjusted, na.rm=TRUE), by=id]
        arms <- bins[,.(armmean=mean(adjusted.cent, na.rm=TRUE)),
                     by=.(id, arm)] 
        arms <- arms[,armmean:=armmean - mean(armmean, na.rm = TRUE)]
        arms[]
}

gcCorrectLoess <- function(counts.mult, gc) {
    lower <- 0
    upper <- quantile(counts.mult, .99, na.rm=TRUE)
    trend.points <- counts.mult > lower & counts.mult < upper
    trend.counts <- counts.mult[trend.points]
    trend.gc <- gc[trend.points]
    num.sample.points <- min(length(trend.counts), 10E3L)
    samp <- sample(1:length(trend.counts), num.sample.points)
    #    pad counts and GC
    include <- c(which(gc==min(gc)), which(gc==max(gc)))
    trend.counts <- c(counts.mult[include], trend.counts[samp])
    trend.gc <- c(gc[include], trend.gc[samp])
    initial.trend <- loess(trend.counts ~ trend.gc)
    i <- seq(min(gc, na.rm=TRUE), max(gc, na.rm=TRUE), by = 0.001)
    final.trend <- loess(predict(initial.trend, i) ~ i)
    counts.pred <- predict(final.trend, gc)
    return(counts.pred)
}

# ---- Define command-line options ----
spec <- matrix(c(
  "datadir", "d", 1, "character",  # --datadir or -d for directory with the data
  "outfile", "o", 1, "character",  # --outfile or -o for output file
  "refdir", "r", 1, "character"  # --refdir or -r for reference directory
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
refdir <- opt$refdir

cat("Datadir", datadir, "\n")
cat("Outfile:", outfile, "\n")
cat("Refdir:", refdir, "\n")

# ---- Load sample files --------
filenames<-list.files(datadir,pattern=".csv",full.names=TRUE)
x <- rbindlist(lapply(filenames, fread))
x <- setDT(x %>% filter(chr != "chrX"))
x[, cov := short + long]
setkey(x, chr, start)

# ----- Load Reference ----
# Construct full path
refbins <- file.path(refdir, "ref20_bins.rda")

# Load ref20_bins
if (file.exists(refbins)) {
  load(refbins)  # loads 'ref20_bins' into environment
  if (!exists("ref20_bins")) {
    stop("ref20_bins object not found after loading: ", refbins)
  }
  refbins <- as.data.table(ref20_bins)  # Convert to data.table
  setkey(refbins, chr, start)
  refbins[,cov:=short+long]
} else {
  stop("ref20_bins.rda not found at: ", refbins)
}

# ---- Compute normalized coverage and arm-level means ----
binsforzscores <- rbind(refbins, x, fill = TRUE)
binsforzscores <- countAndNormalize(binsforzscores, measure = "cov")
armmeans_all <- getArmMeans(binsforzscores)

# ---- Separate reference and sample data ----
ref_ids <- unique(refbins$id)
armmeans_ref <- armmeans_all %>% filter(id %in% ref_ids)
armmeans_samples <- armmeans_all %>% filter(!id %in% ref_ids)

# ---- Compute z-scores ----
ref_stats <- armmeans_ref %>%
  group_by(arm) %>%
  summarize(Mean = mean(armmean), STD = sd(armmean), .groups = "drop")

armmeans_samples <- armmeans_samples %>%
  left_join(ref_stats, by = "arm") %>%
  mutate(zscore = (armmean - Mean) / STD)

# ---- Format output ----
armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p","19q","20p","20q","21q","22q")

armmeans_samples$arm <- factor(armmeans_samples$arm, levels = armlevels)
armmeans_samples[, armvar := paste0("zscore_", arm)]

features_zscores <- dcast(armmeans_samples, id ~ armvar, value.var = "zscore")

# ---- Write output ----
write_csv(features_zscores, outfile)