#!/bin/bash
#SBATCH --job-name=modeling
#SBATCH --output=/well/ludwig/users/cnr137/methylation_model/logs/modeling/%A.out
#SBATCH --error=/well/ludwig/users/cnr137/methylation_model/logs/modeling/%A.err
#SBATCH --time=96:00:00
#SBATCH --mem=600G

###------------------------------------------------- module loading 
#modules

#activate conda env with fetchChromSize, python, numpy, sklearn, xgboost and pandas installed
source /apps/eb/el8/2023a/skylake/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh
conda activate /users/ludwig/cnr137/.conda/envs/epigenetics_env

###------------------------------------------------- flag definition and default definition

while getopts "w:f:o:" flag; do
    case "${flag}" in
        w) WORKDIR="${OPTARG}" ;;       
        f) FEATUREDIR="${OPTARG}" ;;
        o) OUTDIR="${OPTARG}" ;;      
    esac
done


# Set Default variables:

# Make sure a feature dir exists
DEFAULT_FEATUREDIR="./featuredir"
if [[ -z "$FEATUREDIR" ]]; then
    FEATUREDIR="$DEFAULT_FEATUREDIR"
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

# Clean up any existing outputs from previous runs
rm -f "$OUTDIR/FeatureMatrix.csv" "$OUTDIR/tmp_values.csv" "$OUTDIR/features.txt" "$OUTDIR/Target.csv"

###------------------------------------------------------ processing pipeline

#----Step 1: Create feature matrix --------
echo "Preparing feature matrix..."

# Step 1a: Use first bed file to generate feature (column) names
first_file=$(ls "$FEATUREDIR"/*.bed* | head -n 1)
awk '{print $1":"$2":"$3}' "$first_file" | paste -sd',' - > "$OUTDIR/features.txt"
echo "sample_id,$(cat "$OUTDIR/features.txt")" > "$OUTDIR/FeatureMatrix.csv"

# Step 1b: For each file, extract sample name and feature values (fractions in 4th column)
for file in "$FEATUREDIR"/*.bed*; do
    sample_name=$(basename "$file")
    values=$(awk '{print $4}' "$file" | paste -sd',' -)
    echo "$sample_name,$values" >> "$OUTDIR/tmp_values.csv"
done

# Step 1c: Combine all into final matrix
cat "$OUTDIR/tmp_values.csv" >> "$OUTDIR/FeatureMatrix.csv"

# Step 1d: Create target matrix
echo "sample_id,tumour" > "$OUTDIR/Target.csv"
for file in "$FEATUREDIR"/*.bed*; do
    sample_name=$(basename "$file")
    
    if [[ "$sample_name" == *tumour* ]]; then
        label=1
    elif [[ "$sample_name" == *healthy* || "$sample_name" == *cirrhosis* ]]; then
        label=0
    fi
    echo "$sample_name,$label" >> "$OUTDIR/Target.csv"
done

# Step 1e: Cleanup
rm "$OUTDIR/tmp_values.csv" "$OUTDIR/features.txt"

echo "Feature matrix written to $OUTDIR/FeatureMatrix.csv"


#----- Step 2: Evaluate models and feature selection methods ------
python Select_model.py \
    --Featurematrix "$OUTDIR/FeatureMatrix.csv" \
    --Target "$OUTDIR/Target.csv" \
    --output_dir "$OUTDIR"

# sbatch -p long 03-modeling.sh -w /users/ludwig/cnr137 -f /well/ludwig/users/cnr137/methylation_model/generated_samples/features -o /well/ludwig/users/cnr137/methylation_model/generated_samples/model_eval
# sbatch -p long 03-modeling.sh -w /users/ludwig/cnr137 -f /well/ludwig/users/cnr137/methylation_model/generated_samples/features_corr/corr/ -o /well/ludwig/users/cnr137/methylation_model/generated_samples/model_eval/corr
# sbatch -p long 03-modeling.sh -w /users/ludwig/cnr137 -f /well/ludwig/users/cnr137/methylation_model/generated_samples/features_corr_2/corr -o /well/ludwig/users/cnr137/methylation_model/generated_samples/model_eval/corr_2
# sbatch -p long 03-modeling.sh -w /users/ludwig/cnr137 -f /well/ludwig/users/cnr137/methylation_model/generated_samples/features_newsamples -o /well/ludwig/users/cnr137/methylation_model/generated_samples/new_model_eval
# sbatch -p long 03-modeling.sh -w /users/ludwig/cnr137 -f /well/ludwig/users/cnr137/methylation_model/generated_samples/features_newsamples/corr -o /well/ludwig/users/cnr137/methylation_model/generated_samples/new_model_eval_corr
