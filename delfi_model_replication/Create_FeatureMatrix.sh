#!/bin/bash
#SBATCH --job-name=create_featurematrix
#SBATCH --output=/well/ludwig/users/cnr137/CD_cohort/logs/features/%A.out
#SBATCH --error=/well/ludwig/users/cnr137/CD_cohort/logs/features/%A.err
#SBATCH --time=96:00:00
#SBATCH --mem=200G

###------------------------------------------------- module loading and define storage log files

module load R/4.3.2-gfbf-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2

###------------------------------------------------- flag definition

while getopts "w:z:c:o:" flag; do
    case "${flag}" in
        w) WORKDIR="${OPTARG}" ;;
        z) ZFILE="${OPTARG}" ;;
        c) COVFILE="${OPTARG}" ;;   
        o) OUTFILE="${OPTARG}" ;;
    esac
done

###------------------------------------------------------ processing pipeline

Rscript 05-create-featurematrix.r --zfile "$ZFILE" --covfile "$COVFILE" --outfile "$OUTFILE"
echo "Finished making the Feature Matrix!"

#example usage: sbatch -p long Create_FeatureMatrix.sh -w /well/ludwig/users/cnr137/lifted -z /well/ludwig/users/cnr137/lifted/zscores.csv -c /well/ludwig/users/cnr137/lifted/coverage.csv -o /well/ludwig/users/cnr137/lifted/allfeatures.csv
#example usage: sbatch -p long Create_FeatureMatrix.sh -w /well/ludwig/users/cnr137/post_alignment_hightumour_lifted -z /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/zscores_hightumour.csv -c /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/coverage_hightumour.csv -o /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/allfeatures_hightumour.csv
#sbatch -p long Create_FeatureMatrix.sh -w /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references -z /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references/zscores_references.csv -c /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references/coverage_references.csv -o /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references/allfeatures_references.csv
#sbatch -p long Create_FeatureMatrix.sh -w /well/ludwig/users/cnr137/well/ludwig/users/cnr137/CD_cohort -z /well/ludwig/users/cnr137/CD_cohort/zscores_CD.csv -c /well/ludwig/users/cnr137/CD_cohort/coverage_CD.csv -o /well/ludwig/users/cnr137/CD_cohort/allfeatures_CD.csv
