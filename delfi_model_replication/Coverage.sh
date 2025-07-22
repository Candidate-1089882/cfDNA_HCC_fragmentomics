#!/bin/bash
#SBATCH --job-name=calculate_coverage
#SBATCH --output=/well/ludwig/users/cnr137/CD_cohort/logs/coverage/%A.out
#SBATCH --error=/well/ludwig/users/cnr137/CD_cohort/logs/coverage/%A.err
#SBATCH --time=96:00:00
#SBATCH --mem=200G

###------------------------------------------------- module loading and define storage log files

module load R/4.3.2-gfbf-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2

###------------------------------------------------- flag definition

while getopts "w:b:o:" flag; do
    case "${flag}" in
        w) WORKDIR="${OPTARG}" ;;
        b) BINDIR="${OPTARG}" ;;
        o) OUTFILE="${OPTARG}" ;;
    esac
done

###------------------------------------------------------ processing pipeline

Rscript 04-get-coverage.r --datadir "$BINDIR" --outfile "$OUTFILE"
echo "Finished calculating Coverage!"

#example usage: sbatch -p long Coverage.sh -w /well/ludwig/users/cnr137/lifted -b /well/ludwig/users/cnr137/lifted/bindir -o /well/ludwig/users/cnr137/lifted/coverage.csv
#example usage: sbatch -p long Coverage.sh -w /well/ludwig/users/cnr137/post_alignment_hightumour_lifted -b /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/bindir -o /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/coverage_hightumour.csv
#sbatch -p long Coverage.sh -w /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references -b /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references/bindir -o /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references/coverage_references.csv
#sbatch -p long Coverage.sh -w /well/ludwig/users/cnr137/well/ludwig/users/cnr137/CD_cohort -b /well/ludwig/users/cnr137/CD_cohort/bindir -o /well/ludwig/users/cnr137/CD_cohort/coverage_CD.csv -r /well/ludwig/users/cnr137/references
