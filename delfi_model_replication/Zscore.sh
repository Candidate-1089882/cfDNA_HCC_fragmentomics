#!/bin/bash
#SBATCH --job-name=calculate_zscores
#SBATCH --output=/well/ludwig/users/cnr137/CD_cohort/logs/zscore/%A_%a.out
#SBATCH --error=/well/ludwig/users/cnr137/CD_cohort/logs/zscore/%A_%a.err
#SBATCH --time=96:00:00
#SBATCH --mem=200G


###------------------------------------------------- module loading and define storage log files

module load R/4.3.2-gfbf-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2

###------------------------------------------------- flag definition

while getopts "w:b:o:r:" flag; do
    case "${flag}" in
        w) WORKDIR="${OPTARG}" ;;
        b) BINDIR="${OPTARG}" ;;
        o) OUTFILE="${OPTARG}" ;;
        r) REFDIR="${OPTARG}" ;;
    esac
done

###------------------------------------------------------ processing pipeline

Rscript 03-get_Zscores.r --datadir "$BINDIR" --outfile "$OUTFILE" --refdir "$REFDIR"
echo "Finished calculating Zscores!"

#example usage: sbatch -p long Zscore.sh -w /well/ludwig/users/cnr137/lifted -b /well/ludwig/users/cnr137/lifted/bindir -o /well/ludwig/users/cnr137/lifted/zscores.csv -r  /well/ludwig/users/cnr137/references
#example usage: sbatch -p long Zscore.sh -w /well/ludwig/users/cnr137/post_alignment_hightumour_lifted -b /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/bindir -o /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/zscores_hightumour.csv -r /well/ludwig/users/cnr137/references
#sbatch -p long Zscore.sh -w /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references -b /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references/bindir -o /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references/zscores_references.csv -r /well/ludwig/users/cnr137/references
#sbatch -p long Zscore.sh -w /well/ludwig/users/cnr137/well/ludwig/users/cnr137/CD_cohort -b /well/ludwig/users/cnr137/CD_cohort/bindir -o /well/ludwig/users/cnr137/CD_cohort/zscores_CD.csv -r /well/ludwig/users/cnr137/references
