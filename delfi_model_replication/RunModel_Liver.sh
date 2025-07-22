#!/bin/bash
#SBATCH --job-name=run_model
#SBATCH --output=/well/ludwig/users/cnr137/CD_cohort/logs/model/%A.out
#SBATCH --error=/well/ludwig/users/cnr137/CD_cohort/logs/model/%A.err
#SBATCH --time=96:00:00
#SBATCH --mem=200G


###------------------------------------------------- module loading and define storage log files

module load R/4.3.2-gfbf-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2

###------------------------------------------------- flag definition

while getopts "w:f:m:o:" flag; do
    case "${flag}" in
        w) WORKDIR="${OPTARG}" ;;
        f) FEATUREFILE="${OPTARG}" ;;
        m) MODEL="${OPTARG}" ;;   
        o) OUTFILE="${OPTARG}" ;;
    esac
done

###------------------------------------------------------ processing pipeline

Rscript 09-run-model-liver.r --featurefile "$FEATUREFILE" --model "$MODEL" --outfile "$OUTFILE"
echo "Finished making the Predictions!"

#example usage: sbatch -p long RunModel_Liver.sh -w /well/ludwig/users/cnr137/HCC_samples_lifted -f /well/ludwig/users/cnr137/HCC_samples_lifted/allfeatures.csv -m /well/ludwig/users/cnr137/models/delfi-GBM_risk.rds -o /well/ludwig/users/cnr137/HCC_samples_lifted/results_model_GBM.csv
#example usage: sbatch -p long RunModel_Liver.sh -w /well/ludwig/users/cnr137/post_alignment_hightumour_lifted -f /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/allfeatures_hightumour.csv -m /well/ludwig/users/cnr137/models/delfi-GBM_risk.rds  -o /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/results_model_GBM_hightumour.csv
#sbatch -p long RunModel_Liver.sh -w /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references -f /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references/allfeatures_references.csv -m /well/ludwig/users/cnr137/models/delfi-GBM_risk.rds  -o /well/ludwig/users/cnr137/post_alignment_hightumour_lifted/references/results_model_GBM_references.csv
#sbatch -p long RunModel_Liver.sh -w /well/ludwig/users/cnr137/CD_cohort -f /well/ludwig/users/cnr137/CD_cohort/allfeatures_CD.csv -m /well/ludwig/users/cnr137/models/delfi-GBM_risk.rds -o /well/ludwig/users/cnr137/CD_cohort/resultsCD_model_GBM.csv
#sbatch -p long RunModel_Liver.sh -w /well/ludwig/users/cnr137/HCC_validation_DELFI -f /well/ludwig/users/cnr137/HCC_validation_DELFI/validation_features.csv -m /well/ludwig/users/cnr137/models/delfi-GBM_risk.rds -o /well/ludwig/users/cnr137/HCC_validation_DELFI/results_val_model_GBM.csv