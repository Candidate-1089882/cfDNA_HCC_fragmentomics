#!/bin/bash
#SBATCH --job-name=modeling
#SBATCH --output=/well/ludwig/users/cnr137/methylation_model/logs/val/%A.out
#SBATCH --error=/well/ludwig/users/cnr137/methylation_model/logs/val/%A.err
#SBATCH --time=96:00:00
#SBATCH --mem=200G

###------------------------------------------------- module loading 
#modules

#activate conda env with fetchChromSize, python, numpy, sklearn, xgboost and pandas installed
source /apps/eb/el8/2023a/skylake/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh
conda activate /users/ludwig/cnr137/.conda/envs/epigenetics_env

###------------------------------------------------- flag definition and default definition

while getopts "w:c:n:a:b:t:o:" flag; do
    case "${flag}" in
        w) WORKDIR="${OPTARG}" ;;  
        n) FEATUREDIRNC="${OPTARG}" ;;     
        c) FEATUREDIRC="${OPTARG}" ;;
        a) TRAININGC="${OPTARG}" ;;
        b) TRAININGNC="${OPTARG}" ;;
        t) TRAININGTARG="${OPTARG}" ;;
        o) OUTDIR="${OPTARG}" ;;      
    esac
done


# Set Default variables:

# Make sure feature dirs exists
DEFAULT_FEATUREDIRNC="./featuredir"
if [[ -z "$FEATUREDIRNC" ]]; then
    FEATUREDIRNC="$DEFAULT_FEATUREDIRNC"
fi

DEFAULT_FEATUREDIRC="./featuredir"
if [[ -z "$FEATUREDIRC" ]]; then
    FEATUREDIRC="$DEFAULT_FEATUREDIRC"
fi

#If no workdir provided, make a workdir in the current directory
DEFAULT_WORKDIR="./workdir"
if [[ -z "$WORKDIR" ]]; then
    WORKDIR="$DEFAULT_WORKDIR"
fi
mkdir -p "$WORKDIR"
cd "$WORKDIR"

#If no outdir provided, make an outdir in the current directory
DEFAULT_OUTDIR="./outdir"
if [[ -z "$OUTDIR" ]]; then
    OUTDIR="$DEFAULT_OUTDIR"
fi
mkdir -p "$OUTDIR"


###------------------------------------------------------ processing pipeline

#----Step 1: Create feature matrix for non corrected features --------
echo "Preparing feature matrix Non Corrected..."

# Step 1a: Use first bed file to generate feature (column) names
first_file=$(ls "$FEATUREDIRNC"/*.bed* | head -n 1)
awk '{print $1":"$2":"$3}' "$first_file" | paste -sd',' - > "$OUTDIR/featuresNC.txt"
echo "sample_id,$(cat "$OUTDIR/featuresNC.txt")" > "$OUTDIR/FeatureMatrixNC.csv"

# Step 1b: For each file, extract sample name and feature values (fractions in 4th column)
for file in "$FEATUREDIRNC"/*.bed*; do
    sample_name=$(basename "$file")
    values=$(awk '{print $4}' "$file" | paste -sd',' -)
    echo "$sample_name,$values" >> "$OUTDIR/tmp_valuesNC.csv"
done

# Step 1c: Combine all into final matrix
cat "$OUTDIR/tmp_valuesNC.csv" >> "$OUTDIR/FeatureMatrixNC.csv"

echo "Feature matrix written to $OUTDIR/FeatureMatrixNC.csv"

#----Step 2: Create feature matrix for corrected features --------
echo "Preparing feature matrix Corrected..."

# Step 2a: Use first bed file to generate feature (column) names
first_file=$(ls "$FEATUREDIRC"/*.bed* | head -n 1)
awk '{print $1":"$2":"$3}' "$first_file" | paste -sd',' - > "$OUTDIR/featuresC.txt"
echo "sample_id,$(cat "$OUTDIR/featuresC.txt")" > "$OUTDIR/FeatureMatrixC.csv"

# Step 2b: For each file, extract sample name and feature values (fractions in 4th column)
for file in "$FEATUREDIRC"/*.bed*; do
    sample_name=$(basename "$file")
    values=$(awk '{print $4}' "$file" | paste -sd',' -)
    echo "$sample_name,$values" >> "$OUTDIR/tmp_valuesC.csv"
done

# Step 2c: Combine all into final matrix
cat "$OUTDIR/tmp_valuesC.csv" >> "$OUTDIR/FeatureMatrixC.csv"

echo "Feature matrix written to $OUTDIR/FeatureMatrixC.csv"

#----Step 3: Create target matrix --------
echo "sample_id,tumour" > "$OUTDIR/Target.csv"
for file in "$FEATUREDIRNC"/*.bed*; do
    sample_name=$(basename "$file")
    
    if [[ "$sample_name" == *HCC* ]]; then
        label=1
    elif [[ "$sample_name" == *control* ]]; then
        label=0
    fi
    echo "$sample_name,$label" >> "$OUTDIR/Target.csv"
done

# Step 1e: Cleanup
rm "$OUTDIR/tmp_valuesNC.csv" "$OUTDIR/featuresNC.txt"
rm "$OUTDIR/tmp_valuesC.csv" "$OUTDIR/featuresC.txt"


#----- Step 4: Validate the models for non corrected features ------
python Validation_model_noncorr.py \
    --Featurematrix "$TRAININGNC" \
    --Target "$TRAININGTARG" \
    --ValidationFeatures "$OUTDIR/FeatureMatrixNC.csv" \
    --ValidationTarget "$OUTDIR/Target.csv" \
    --output_dir "$OUTDIR"

#----- Step 5: Validate the models for corrected features ------
python Validation_model_corr.py \
    --Featurematrix "$TRAININGC" \
    --Target "$TRAININGTARG" \
    --ValidationFeatures "$OUTDIR/FeatureMatrixC.csv" \
    --ValidationTarget "$OUTDIR/Target.csv" \
    --output_dir "$OUTDIR"

# sbatch -p long 05-validation-models.sh -w /users/ludwig/cnr137 -n /well/ludwig/users/cnr137/methylation_model/validation_samples/not_corr/features -c /well/ludwig/users/cnr137/methylation_model/validation_samples/corr/features_corr -o /well/ludwig/users/cnr137/methylation_model/validation_samples_feat -a /well/ludwig/users/cnr137/methylation_model/generated_samples/model_eval/FeatureMatrix.csv -b /well/ludwig/users/cnr137/methylation_model/generated_samples/model_eval/corr/FeatureMatrix.csv -t /well/ludwig/users/cnr137/methylation_model/generated_samples/model_eval/Target.csv


