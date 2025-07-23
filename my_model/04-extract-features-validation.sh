#!/bin/bash
#SBATCH --job-name=extract_features_val
#SBATCH --output=/well/ludwig/users/cnr137/methylation_model/logs/binning_val/%A.out
#SBATCH --error=/well/ludwig/users/cnr137/methylation_model/logs/binning_val/%A.err
#SBATCH --time=96:00:00
#SBATCH --mem=600G


###------------------------------------------------- module loading 
# controls: /well/ludwig/users/ikb109/OAC_Immuno_Trial/TAPS_cfDNA/PerReadCalls
# cfDNA calls (validation): /well/ludwig/users/ikb109/Deliver_calls/1.2/PerReadCalls
# solid tissue samples: /well/ludwig/users/ikb109/Tissue_map/Results/1.2/PerReadCalls

#modules
module load BEDTools/2.31.0-GCC-12.3.0 
module load BEDOPS/2.4.41-foss-2023a

#activate conda env with fetchChromSize, python, statsmodels, matplotlib, seaborn and numpy installed
source /apps/eb/el8/2023a/skylake/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh
conda activate /users/ludwig/cnr137/.conda/envs/epigenetics_env

###------------------------------------------------- flag definition and default definition

while getopts "w:s:t:r:f:g:" flag; do
    case "${flag}" in
        w) WORKDIR="${OPTARG}" ;;       
        s) SAMPLEDIR="${OPTARG}" ;;
        t) TEMPDIR="${OPTARG}" ;;
        r) REFDIR="${OPTARG}" ;; 
        f) FEATUREDIR="${OPTARG}" ;;
        g) GC_CORRECTION="${OPTARG}" ;;      
    esac
done

# Set Default variables:

# Make sure a feature dir (output dir) exists
DEFAULT_FEATUREDIR="./featuredir"
if [[ -z "$FEATUREDIR" ]]; then
    FEATUREDIR="$DEFAULT_FEATUREDIR"
fi
mkdir -p "$FEATUREDIR"

# Make sure tempdir exists
# Use tempdir because can't write in sampledir
DEFAULT_TEMPDIR="./tempdir"
if [[ -z "$TEMPDIR" ]]; then
    TEMPDIR="$DEFAULT_TEMPDIR"
fi
mkdir -p "$TEMPDIR"

# Set default for GC correction (disabled by default)
DEFAULT_GC_CORRECTION="false"
if [[ -z "$GC_CORRECTION" ]]; then
    GC_CORRECTION="$DEFAULT_GC_CORRECTION"
fi

#If no workdir provided, make a workdir in the current directory
DEFAULT_WORKDIR="./workdir"
if [[ -z "$WORKDIR" ]]; then
    WORKDIR="$DEFAULT_WORKDIR"
fi
mkdir -p "$WORKDIR"
cd "$WORKDIR"


###------------------------------------------------------ processing pipeline

#----- Step 1: Create references (if needed) -------
# Step 1a: Create genome windows
WINDOWS_BED="$REFDIR"/"windows.bed"

if [[ ! -f $WINDOWS_BED ]]; then
    echo "Creating genome windows..."
    fetchChromSizes hg38 | grep -v '_' | \
      awk -v FS="\t" -v OFS="\t" '{ print $1, "0", $2 }' | \
      sort-bed - | bedops --chop 5000000 - > "$WINDOWS_BED"
fi

# Step 2b: Make sure GC content per window is available if GC correction is needed
GC_BED="$REFDIR/gc_content_windows.bed"

if [[ ! -f $GC_BED ]]; then
    echo "Calculating GC content..."
    bedtools nuc -fi "$REFDIR/hg38.fa" -bed "$WINDOWS_BED" | \
      awk 'NR > 1 {print $1, $2, $3, $5}' OFS="\t" > "$GC_BED"
fi

#------- Step 2: Convert BED files into desired format ----------
# org: chr	start	end	read_id	mapq	orientation	insert_size	flag	num_cpg	num_mod	mod_cps	unmod_cpgs	snp_cpgs
# goal: chr    start    end    num_CpGs    methylation_level
mkdir -p "$TEMPDIR/tmp/"
for file in "$SAMPLEDIR"/*.bed.gz; do
  base=$(basename "$file" .bed.gz)
  echo "Converting $base..."
  mkdir -p "$TEMPDIR/tmp/${base}/"
  CONVERTED="$TEMPDIR/tmp/${base}/${base}_converted.bed"

  zcat "$file" | awk -v FS="\t" -v OFS="\t" '
  BEGIN { OFS="\t" }
  !/^#/ {
    mapq = $5;
    num_CpGs = $9;
    if (mapq > 10 && num_CpGs ~ /^[0-9]+$/ && num_CpGs > 0) {
      methyl_level = $10/num_CpGs;
      print $1, $2, $3, num_CpGs, methyl_level
    }
  }' > "$CONVERTED"

done

#------- Step 3: Process BED files and calculate counts of reads and hypomethylated reads per window -------
for file in "$SAMPLEDIR"/*.bed.gz; do
  base=$(basename "$file" .bed.gz)
  mkdir -p "$TEMPDIR/tmp/${base}/"
  CONVERTED="${base}_converted.bed"
  cd "$TEMPDIR/tmp/${base}/" 

  # Step 3a: Filter reads with >=2 CpGs (col 4 = num_CpGs)
  FILTERED="${base}_filtered.bed"
  awk -v FS="\t" -v OFS="\t" '$4 >= 2 { print $1, $2, $3, ".", $4, $5 }' "$CONVERTED" | sort-bed - > "$FILTERED"

  # Step 3b: Select hypomethylated reads from BED (hypomethylated = methylation level =< 0.35)
  HYPO="${base}_hypo.bed"
  awk -v FS="\t" -v OFS="\t" '$4 >= 3 && $5 <= 0.35 { print $1, $2, $3, ".", $4, $5 }' "$CONVERTED" | sort-bed - > "$HYPO"

  # Step 3c: Count total reads in each window
  TOTAL_COUNT="${base}_total.count"
  bedmap --count "$WINDOWS_BED" "$FILTERED" > "$TOTAL_COUNT"

  # Step 3d: Count hypomethylated reads in each window
  HYPO_COUNT="${base}_hypo.count"
  bedmap --count "$WINDOWS_BED" "$HYPO" > "$HYPO_COUNT"

  # Step 3e: Calculate hypomethylation fraction per window 
  FRACTION="${base}_hypo_fraction.bed"
  paste "$WINDOWS_BED" "$HYPO_COUNT" "$TOTAL_COUNT" | \
    awk -v OFS="\t" '{ if ($5 == 0) frac = "NA"; else frac = $4 / $5; print $1, $2, $3, frac }' > "$FRACTION"

  mv "$FRACTION" "$FEATUREDIR"
  echo "Done: $FRACTION"
done
 
#--------- Step 4: peform GC correction if needed -----------
if [[ "$GC_CORRECTION" == "true" ]]; then
  echo "Performing GC correction..."

  for file in "$SAMPLEDIR"/*.bed.gz; do
    base=$(basename "$file" .bed.gz)
    cd "$TEMPDIR/tmp/${base}/" 
  
    #Step 4a: Join GC content with counts
    TOTAL_COUNT="${base}_total.count"
    HYPO_COUNT="${base}_hypo.count"
    paste "$GC_BED" "$HYPO_COUNT" "$TOTAL_COUNT" | \
      awk -v OFS="\t" '{chrom = $1; start = $2; end = $3; gc = $4; hypo = $5; total = $6; print chrom, start, end, gc, hypo, total}' > "${base}_gc_counts.tsv"
  
    #Step 4b: run Python script with GC correction based on linear reg
    python "$WORKDIR/GC_correction.py" "${base}_gc_counts.tsv" "${base}_corrected_hypo_fraction.bed"
    mv "${base}_corrected_hypo_fraction.bed" "$FEATUREDIR/corr/"

    # Check
    echo "GC-corrected output saved to $FEATUREDIR/corr/${base}_corrected_hypo_fraction.bed"
  
  done
fi

#cleanup
# rm -r $TEMPDIR

# Example usage:
# sbatch -p long 04-extract-features-validation.sh -w /users/ludwig/cnr137 -s /well/ludwig/users/ikb109/Deliver_calls/1.2/PerReadCalls -t /well/ludwig/users/cnr137/methylation_model/validation_samples/not_corr -r /well/ludwig/users/cnr137/references -f /well/ludwig/users/cnr137/methylation_model/validation_samples/not_corr/features -g false
# sbatch -p long 04-extract-features-validation.sh -w /users/ludwig/cnr137 -s /well/ludwig/users/ikb109/Deliver_calls/1.2/PerReadCalls -t /well/ludwig/users/cnr137/methylation_model/validation_samples/corr -r /well/ludwig/users/cnr137/references -f /well/ludwig/users/cnr137/methylation_model/validation_samples/corr/features_corr -g true