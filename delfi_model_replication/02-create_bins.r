# ---- Description ----
# This script processes cell-free DNA fragment data for downstream analysis.
# It performs the following steps:
# 1. Loads required libraries and reference files
# 2. Loads and cleans the fragment file (removes blacklisted regions)
# 3. Applies GC-content correction based on a target reference
# 4. Bins the corrected fragments into genomic bins
# 5. Outputs the final binned summary to a CSV file
# Usage (from command line):
# Rscript script.R --id <fragment_file.rds> --outdir <output_directory> --refdir <reference_directory>

# ---- Install and load packages ----
required_packages <- c("getopt", "data.table")
bioconductor_packages <- c("GenomicRanges")

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

#Set data.table to use only 1 thread
setDTthreads(1)

# ---- Define needed functions -----
gcCorrectTarget <- function(fragments, ref, bychr=TRUE){
    fragments[, gc := round(gc, 2)]
    if(bychr) {
        DT.gc <- fragments[,.(n=.N), by=.(gc, chr)]
        DT.gc <- DT.gc[gc >= .20 & gc <= .80]
        DT.gc <- DT.gc[order(gc, chr)]
    } else {
        DT.gc <- fragments[,.(n=.N), by=gc]
        DT.gc <- DT.gc[gc >= .20 & gc <= .80]
        DT.gc <- DT.gc[order(gc)]
    }
# setkey(mediandt, gc, seqnames)

    if(bychr) {
        setkey(DT.gc, gc, chr)
        setkey(ref, gc, chr)
    } else {
        setkey(DT.gc, gc)
        setkey(ref, gc)
    }
#     DT.gc <- DT.gc[ref][order(chr, gc)]
    DT.gc <- DT.gc[ref]
    DT.gc[,w:=target/n]
    if(bychr) { 
        fragments[DT.gc, on= .(chr, gc), weight := i.w]
    }
    else fragments[DT.gc, on= .(gc), weight := i.w]
    fragments <- fragments[!is.na(weight)]
    fragments[]
}

## Be sure to setname(fragments, c("seqnames", "gc"), c("chr", "fraggc"))
binFrags <- function(fragments, bins, cutoff=150,
                     chromosomes=paste0("chr",c(1:22, "X"))) {
    setkey(bins, chr, start, end)
    fragbins <- foverlaps(fragments[chr %in% chromosomes],
                          bins, type="within", nomatch=NULL)
    bins2 <- fragbins[,.(arm=unique(arm), gc=gc[1], map=map[1],
                         short = sum(w >= 100 & w <= cutoff ),
                         long = sum(w > cutoff & w <= 250),
                         short.cor = sum(weight[w >= 100 & w <= cutoff]),
                         long.cor = sum(weight[w > cutoff & w <= 250]),
                         ultrashort = sum(w < 100),
                         ultrashort.cor = sum(weight[w < 100]),
                         multinucs = sum(w > 250),
                         multinucs.cor = sum(weight[w > 250]),
                         mediansize = as.integer(median(w)),
                         frag.gc = mean(fraggc)),
            by=.(chr, start, end)]

    setkey(bins2, chr, start, end)
    bins2 <- bins2[bins]
    bins2 <- bins2[is.na(i.gc), which(grepl("short|long|multi", colnames(bins2))):=0]
    bins2[,`:=`(gc=i.gc, map=i.map, arm=i.arm)]
    bins2[,which(grepl("^i.", colnames(bins2))):=NULL]
    bins2[, bin:=1:.N]
    setcolorder(bins2, c("chr", "start", "end", "bin"))
    bins2[]
}

# ---- Define command-line options ----
spec <- matrix(c(
  "id",     "i", 1, "character",  # --id or -i for frag file
  "outdir", "o", 1, "character",  # --outdir or -o for output dir
  "refdir", "r", 1, "character"  # --refdir or -r for reference directory
), byrow = TRUE, ncol = 4)

# Parse the command line arguments
opt <- getopt(spec)

# ---- Validate input ----
if (is.null(opt$id) || is.null(opt$outdir) || is.null(opt$refdir)) {
  cat(getopt(spec, usage = TRUE))
  quit(status = 1)
}

fragfile <- opt$id
outdir <- opt$outdir
refdir <- opt$refdir
sample_name <- gsub(".rds", "", basename(fragfile))

cat("Sample:", sample_name, "\n")
cat("FRAG file:", fragfile, "\n")
cat("Reference directory:", refdir, "\n")

# ---- Load reference files from refdir ----

# Construct full paths
filters_file <- file.path(refdir, "filters.hg19.rda")
bins_file <- file.path(refdir, "bins_5mb.csv")
target_file <- file.path(refdir, "target20.rda")

# Load filters
if (file.exists(filters_file)) {
  load(filters_file)  # loads 'filters.hg19' into environment
  if (!exists("filters.hg19")) {
    stop("filters.hg19 object not found after loading: ", filters_file)
  }
  filters <- as.data.table(filters.hg19)
  setnames(filters, "seqnames", "chr")
} else {
  stop("filters.hg19.rda not found at: ", filters_file)
}

# Load bins
if (file.exists(bins_file)) {
  bins <- fread(bins_file)
  bins <- bins[map >= 0.90 & gc >= 0.30]
} else {
  stop("bins_5mb.csv not found at: ", bins_file)
}

# Load target20
if (file.exists(target_file)) {
  load(target_file)  # loads 'target20' into environment
  if (!exists("target20")) {
    stop("target20 object not found after loading: ", target_file)
  }
  target <- as.data.table(target20)  # Convert to data.table
  setnames(target, c("seqnames", "gcmed"), c("chr", "target"))
} else {
  stop("target20.rda not found at: ", target_file)
}

# ---- Load fragment file and filter by blacklist ----

# Load fragments
if (!file.exists(fragfile)) {
  stop("Fragment file not found: ", fragfile)
}
fragments <- readRDS(fragfile)

# Keep only chr1–22 and chrX
fragments <- keepSeqlevels(fragments, paste0("chr", c(1:22, "X")), pruning.mode = "coarse")

# Convert to data.table and standardize naming
fragments <- as.data.table(fragments)
setnames(fragments, "seqnames", "chr")

# Filter blacklisted regions using filters.hg19
#If no overlap is found for a fragment, the blacklist columns (start, end, name, etc.) are NA in that row.
setkey(filters, chr, start, end)
fragdt <- foverlaps(fragments, filters, type = "any")

# Keep only fragments not overlapping blacklist (i.e. NA in 'start')
fragdt <- fragdt[is.na(start)]

# Clean up unused columns from blacklist join
fragdt[, `:=`(start = NULL, end = NULL, width = NULL,
              strand = NULL, name = NULL, score = NULL)]

# Restore fragment columns
setnames(fragdt, c("i.start", "i.end", "i.width", "i.strand"),
         c("start", "end", "width", "strand"))

# ---- GC Correction----

# Run GC correction
fragdt <- gcCorrectTarget(fragdt, target)

# ---- Binning ----
# Ensure chromosome levels match expected order
bins[, chr := factor(chr, paste0("chr", c(1:22, "X")))]
fragdt[, chr := factor(chr, paste0("chr", c(1:22, "X")))]

# Rename GC column in fragments
setnames(fragdt, "gc", "fraggc")

# Set keys for overlap
setkey(fragdt, chr, start, end)
setkey(bins, chr, start, end)

# Assign a unique bin number to each row
bins[, bin := 1:.N]

# Bin the corrected fragments using the defined genomic bins
bins2 <- binFrags(fragdt, bins)

#cleanup
bins2[, id := gsub("_downsamp", "", sample_name)]
setcolorder(bins2, c("id", "chr", "start", "end", "bin"))

# ---- Write output ----
filename <- file.path(outdir, paste0(sample_name, "_5mb.csv"))
fwrite(bins2, filename)

# ---- Exit R ----
q("no")
