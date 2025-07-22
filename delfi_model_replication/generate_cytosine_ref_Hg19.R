###------package installing and loading--------

bioc_packages <- c("GenomicRanges", "Biostrings", "BSgenome", "GenomeInfoDb", "BSgenome.Hsapiens.UCSC.hg19")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  suppressMessages(library(pkg, character.only = TRUE))
}

library(BSgenome.Hsapiens.UCSC.hg19)  
library(GenomicRanges)
library(Biostrings)

###--------main------

# Load the genome
genome <- BSgenome.Hsapiens.UCSC.hg19

# Find cytosine and guanine positions
find_cg <- function(seqnames, chunk_size = 1000000, all_seqlevels = NULL) {
  seq <- genome[[seqnames]]  # Get the actual sequence
  cat("Processing chromosome:", seqnames, "\n")  # Debug print

  seq_len <- length(seq)  # Length of the chromosome
  positions <- integer(0)  # Initialize an empty vector for positions
  
  # Loop through in chunks of size chunk_size
  for (start_pos in seq(1, seq_len, by = chunk_size)) {
    end_pos <- min(start_pos + chunk_size - 1, seq_len)
    seq_chunk <- as.character(seq[start_pos:end_pos])  # Get the chunk of the sequence
    
    # Split the chunk into individual bases
    seq_str <- strsplit(seq_chunk, "")[[1]]
    
    # Find positions of C and G
    c_pos <- which(seq_str == "C") + (start_pos - 1)  # Adjust positions based on chunk start
    g_pos <- which(seq_str == "G") + (start_pos - 1)  # Adjust positions based on chunk start
    
    positions <- c(positions, c_pos, g_pos)
  }
  
  cat("Found", length(positions), "C's and G's on chromosome", seqnames, "\n")  # Debug

  # Return GRanges object
  gr <- GRanges(seqnames = rep(seqnames, length(positions)),
          ranges = IRanges(start = positions, width = 1))
  
  if (!is.null(all_seqlevels)) {
    seqlevels(gr) <- all_seqlevels
  }
  return(gr)
}

# Apply to all chromosomes (e.g., chr1 to chr22, X, Y)
chromosomes <- names(genome)[grepl("^chr[0-9XY]+$", names(genome))]
print(chromosomes)

# Set the all_seqlevels (important for merging later without warnings)
all_seqlevels <- seqnames(genome)
cg_list <- lapply(chromosomes, function(chr) find_cg(chr, chunk_size = 1000000, all_seqlevels = all_seqlevels))
cytosine_ref <- do.call(c, cg_list)

warnings()
# Save it
saveRDS(cytosine_ref, file = "cytosine_ref_hg19.rds", compress = TRUE)