#!/bin/bash
#SBATCH --time=96:00:00            # Time limit
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --cpus-per-task=8          # Number of CPU cores per task
#SBATCH --mem=60G                   # Total memory per node
#SBATCH --ntasks=1                # Number of tasks: 59 BAM files
#SBATCH --array=1-59             # Job array with task indices (59 or 119 or 2 or 148)
#SBATCH --mail-type=ALL            # Notifications for all job events
#SBATCH --job-name=Post_alignment  # Job name
#SBATCH -o /well/ludwig/users/cnr137/CD_cohort/logs/post_alignment/post_alignment_%A_%a.out  # Standard output file with job ID and array task ID
#SBATCH -e /well/ludwig/users/cnr137/CD_cohort/logs/post_alignement/post_alignment_%A_%a.err  # Standard error file with job ID and array task ID

###------------------------------------------------- module loading
module purge
module load Miniforge3/24.1.2-0 
module load SAMtools/1.18-GCC-12.3.0

#activate conda env with sambamba and ucsc-liftover installed
source /apps/eb/el8/2023a/skylake/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh
conda activate /users/ludwig/cnr137/.conda/envs/fragmentomics

###------------------------------------------------- flag definition

#check for provided options and assign them to variables
while getopts "b:w:n:" flag; do
    case "${flag}" in
        b) BAMDIR="${OPTARG}" ;;
        w) WORKDIR="${OPTARG}" ;;
        n) num_samples="${OPTARG}" ;;
    esac
done

# Set Default variables:

#If no number of samples given: count number of BAM files in bamdir
if [[ -z "$num_samples" ]]; then
    if [[ -d "$BAMDIR" ]]; then
        num_samples=$(find "$BAMDIR" -maxdepth 1 -name "*.bam" | wc -l)
    else
        echo "Error: BAM directory '$BAMDIR' does not exist."
        exit 1
    fi
fi

#If no workdir provided, make a workdir in the current directory
DEFAULT_WORKDIR="./workdir"
if [[ -z "$WORKDIR" ]]; then
    WORKDIR="$DEFAULT_WORKDIR"
fi
mkdir -p "$WORKDIR"

#directory of BAM files:
#/well/ludwig/processed/Song_lab/HCC_cfDNA/TAPS_cfDNA/Results/Rastair/1.0/Alignments


###----------------------------------------------------- set other variables

BEDDIR="${WORKDIR}/BEDfiles" 
SCRATCHDIR="${WORKDIR}/tmp" 
mkdir -p $BEDDIR $SCRATCHDIR

###----------------------------------------------------- get sample names and check if already processed

samplepath=$(find $BAMDIR -maxdepth 1 -name "*.bam" | \
    sort | \
    head -n $SLURM_ARRAY_TASK_ID | \
    tail -n 1)

sample=$(basename $samplepath .bam)

if [ -f $BEDDIR/$sample.bed ]; then
   echo Sample $sample has been processed. Exiting.
   exit 0
fi

###------------------------------------------------------ processing pipeline

# Flagging and removing duplicates
# ulimit: set higher file limit for sambamba
echo "Flagging duplicates in $sample with sambamba..."
CORES=8
sambamba markdup -t $CORES -l 1 --tmpdir="$SCRATCHDIR" \
    --hash-table-size=10000000 \
    --overflow-list-size=10000000 \
    --sort-buffer-size=2048 \
    -r "$samplepath" "$SCRATCHDIR/${sample}_rmdup.bam"

# Step 1: Sort by read name
samtools sort -n \
  -@ $CORES \
  -m 2G \
  -T "$SCRATCHDIR/${sample}" \
  -o "$SCRATCHDIR/${sample}_sorted.bam" \
  "$SCRATCHDIR/${sample}_rmdup.bam" || { echo "Error: samtools sort failed."; exit 1; }

# Step 2: Fix mates and create a fixed BAM file
samtools fixmate -r "$SCRATCHDIR/${sample}_sorted.bam" "$SCRATCHDIR/${sample}_fixmate.bam" || { echo "Error: samtools fixmate failed."; exit 1; }

# Step 3: Convert to BEDPE format
bamToBed -bedpe -i "$SCRATCHDIR/${sample}_fixmate.bam" > "$SCRATCHDIR/${sample}.raw.bedpe" || { echo "Error: bamToBed conversion failed."; exit 1; }

# Step 4: liftOver from hg38 to hg19 (for downstream analysis purposes)
# make a BEDfile first for liftOver
# Get certain columns: Chromosome, Start of first read, End of second read, and Mapping quality
awk '($2 != -1) && ($5 != -1) && ($1 == $4)' "$SCRATCHDIR/${sample}.raw.bedpe" | \
cut -f 1,2,6,8 > "$SCRATCHDIR/${sample}.prelift.bed"

CHAIN_FILE="/well/ludwig/users/cnr137/tools/hg38ToHg19.over.chain.gz"
LIFTOVER_OUT="$SCRATCHDIR/${sample}.lifted.bed"
LIFTOVER_UNMAPPED="$SCRATCHDIR/${sample}.unmapped.bed"

# Pass the full BEDPE file to liftOver
echo "Lifting over $sample from hg38 to hg19..."
liftOver "$SCRATCHDIR/${sample}.prelift.bed" \
         "$CHAIN_FILE" \
         "$LIFTOVER_OUT" \
         "$LIFTOVER_UNMAPPED" || { echo "Error: liftOver failed."; exit 1; }


# Step 5: Sort the lifted BED file
sort -k1,1 -k2,2n "$LIFTOVER_OUT" > "$BEDDIR/$sample.bed" || { echo "Error: Sorting lifted BED failed."; exit 1; }

# Cleanup
#rm -r "$SCRATCHDIR/tmp"

echo "Successfully processed $sample, output saved to $BEDDIR/$sample.bed"

#usage: sbatch -p long post_alignment.sh -b /well/ludwig/processed/Song_lab/HCC_cfDNA/TAPS_cfDNA/Results/Rastair/1.0/Alignments -w /well/ludwig/users/cnr137/post_alignment_lifted -n 59

# sbatch -p long post_alignment.sh -b /well/ludwig/projects/processed/Lu_lab/OAC_immuno_Trial/TAPS_cfDNA/Results/1.6.1/Alignments -w /well/ludwig/users/cnr137/post_alignment_hightumour_lifted -n 119
# sbatch -p long post_alignment.sh -b /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references -w /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references -n 2
# sbatch -p long post_alignment.sh -b /well/ludwig/projects/processed/Lu_lab/OAC_immuno_Trial/TAPS_cfDNA/CD/Results/1.1/Alignments -w /well/ludwig/users/cnr137/CD_cohort/post_alginment -n 148
