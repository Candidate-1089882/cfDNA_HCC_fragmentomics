#!/bin/bash
#SBATCH --job-name=make_bins
#SBATCH --output=/well/ludwig/users/cnr137/CD_cohort/logs/bins/%A_%a.out
#SBATCH --error=/well/ludwig/users/cnr137/CD_cohort/logs/bins/%A_%a.err
#SBATCH --time=96:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-59


###------------------------------------------------- module loading and define storage log files

module load R/4.3.2-gfbf-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2

###------------------------------------------------- flag definition

while getopts "w:n:f:b:r:" flag; do
    case "${flag}" in
        w) WORKDIR="${OPTARG}" ;;
        n) num_samples="${OPTARG}" ;;
        f) FRAGDIR="${OPTARG}" ;;
        b) BINDIR="${OPTARG}" ;;
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

#If no number of samples given: count number of rds files in fragdir
if [[ -z "$num_samples" ]]; then
    if [[ -d "$FRAGDIR" ]]; then
        num_samples=$(find "$FRAGDIR" -maxdepth 1 -name "*.rds" | wc -l)
    else
        echo "Error: BED directory '$FRAGDIR' does not exist."
        exit 1
    fi
fi

#If no bindir provided, make a bindir in the workdir
DEFAULT_BINDIR="${WORKDIR}/bindir" 
if [[ -z "$BINDIR" ]]; then
    BINDIR="$DEFAULT_BINDIR"
fi
mkdir -p "$BINDIR"


###----------------------------------------------------- get sample names and check if already processed

samplepath=$(find $FRAGDIR -maxdepth 1 -name "*.rds" | \
    sort | \
    head -n $SLURM_ARRAY_TASK_ID | \
    tail -n 1)

sample=$(basename $samplepath | awk '{ gsub(".rds", "") ; print $0}')

if [ -f $BINDIR/${sample}_5mb.csv ]; then
   echo Sample $sample has been processed. Exiting.
   exit 0
fi

###----------------------------------------------------- set other variables
frag_file="$samplepath"

###------------------------------------------------------ processing pipeline

Rscript 02-create_bins.r --id "$frag_file" --outdir "$BINDIR" --refdir "$REFDIR" 
echo "Finished creating fragments!"

#example usage: sbatch -p long bins.sh -w /well/ludwig/users/cnr137/lifted  -n 59 -f /well/ludwig/users/cnr137/bed_to_granges/fragdir_lifted -b /well/ludwig/users/cnr137/lifted/bindir -r /well/ludwig/users/cnr137/references

# sbatch -p long bins.sh -w /well/ludwig/users/cnr137/CD_cohort  -n 148 -f /well/ludwig/users/cnr137/CD_cohort/fragdir -b /well/ludwig/users/cnr137/CD_cohort/bindir -r /well/ludwig/users/cnr137/references
