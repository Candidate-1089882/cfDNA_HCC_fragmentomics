#!/bin/bash
#SBATCH --job-name=bed_to_granges
#SBATCH --output=%x_logs/%A_%a.out
#SBATCH --error=%x_logs/%A_%a.err
#SBATCH --time=96:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-59
#SBATCH --mail-type=ALL            # Notifications for all job events


###------------------------------------------------- module loading and define storage log files
mkdir -p "${SLURM_JOB_NAME}_logs"

module load R/4.3.2-gfbf-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2

###------------------------------------------------- flag definition

while getopts "b:w:n:f:r:" flag; do
    case "${flag}" in
        b) BEDDIR="${OPTARG}" ;;
        w) WORKDIR="${OPTARG}" ;;
        n) num_samples="${OPTARG}" ;;
        f) FRAGDIR="${OPTARG}" ;;
        r) REFDIR="${OPTARG}" ;;
    esac
done

# Set Default variables:

#If no workdir provided, make a workdir in the current directory
DEFAULT_WORKDIR="./workdir"
if [[ -z "$WORKDIR" ]]; then
    WORKDIR="$DEFAULT_WORKDIR"
fi
mkdir -p "$WORKDIR"

#If no number of samples given: count number of BED files in beddir
if [[ -z "$num_samples" ]]; then
    if [[ -d "$BEDDIR" ]]; then
        num_samples=$(find "$BEDDIR" -maxdepth 1 -name "*.bed" | wc -l)
    else
        echo "Error: BED directory '$BEDDIR' does not exist."
        exit 1
    fi
fi

#If no fragdir provided, make a fragdir in the workdir
DEFAULT_FRAGDIR="${WORKDIR}/fragdir" 
if [[ -z "$FRAGDIR" ]]; then
    FRAGDIR="$DEFAULT_FRAGDIR"
fi
mkdir -p "$FRAGDIR"


###----------------------------------------------------- get sample names and check if already processed

samplepath=$(find $BEDDIR -maxdepth 1 -name "*.bed" | \
    sort | \
    head -n $SLURM_ARRAY_TASK_ID | \
    tail -n 1)

sample=$(basename $samplepath .bed)

if [ -f $FRAGDIR/$sample.rds ]; then
   echo Sample $sample has been processed. Exiting.
   exit 0
fi

###----------------------------------------------------- set other variables
FIGUREDIR="${WORKDIR}/figures/${sample}"  # Define the figures directory within the output directory
mkdir -p "$FIGUREDIR"

bed_file="$samplepath"

###------------------------------------------------------ processing pipeline

Rscript 01-bed_to_granges.r --id "$bed_file" --outdir "$FRAGDIR" --refdir "$REFDIR" --figuredir "$FIGUREDIR"
echo "Finished creating fragments!"

#example usage: sbatch -p long bed_to_granges.sh -b /well/ludwig/users/cnr137/post_alignment_lifted/BEDfiles -w /well/ludwig/users/cnr137/bed_to_granges_lifted -n 59 -f /well/ludwig/users/cnr137/bed_to_granges/fragdir_lifted -r /well/ludwig/users/cnr137/references