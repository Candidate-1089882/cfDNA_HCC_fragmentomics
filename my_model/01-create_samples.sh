#!/bin/bash
#SBATCH --job-name=create_samples
#SBATCH --output=/well/ludwig/users/cnr137/methylation_model/logs/samples/%A.out
#SBATCH --error=/well/ludwig/users/cnr137/methylation_model/logs/samples/%A.err
#SBATCH --time=96:00:00
#SBATCH --mem=600G


###------------------------------------------------- module loading 
# controls: /well/ludwig/users/ikb109/OAC_Immuno_Trial/TAPS_cfDNA/PerReadCalls
# cfDNA calls (validation): /well/ludwig/users/ikb109/Deliver_calls/1.2/PerReadCalls
# solid tissue samples: /well/ludwig/users/ikb109/Tissue_map/Results/1.2/PerReadCalls

#modules
module load BEDTools/2.31.0-GCC-12.3.0 

#activate conda env with python3.11.13 and numpy2.2.6 installed
source /apps/eb/el8/2023a/skylake/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh
conda activate /users/ludwig/cnr137/.conda/envs/epigenetics_env
###------------------------------------------------- flag definition and default definition

while getopts "w:c:t:o:" flag; do
    case "${flag}" in
        w) WORKDIR="${OPTARG}" ;;       # workdir
        c) CFDNA_DIR="${OPTARG}" ;;     # cfDNA directory
        t) TISSUE_DIR="${OPTARG}" ;;    # tissue directory
        o) OUTDIR="${OPTARG}" ;;        # output directory
    esac
done

# Set Default variables:

#If no workdir provided, make a workdir in the current directory
DEFAULT_WORKDIR="./workdir"
if [[ -z "$WORKDIR" ]]; then
    WORKDIR="$DEFAULT_WORKDIR"
fi
mkdir -p "$WORKDIR"

#If no outdir provided, make an outdir in the current directory
DEFAULT_OUTDIR="./outdir"
if [[ -z "$OUTDIR" ]]; then
    OUTDIR="$DEFAULT_OUTDIR"
fi
mkdir -p "$OUTDIR"

###------------------------------------------------------ processing pipeline

# Filter cfDNA files to only R1
# NR=1: ignore the header, col 9 contains flags
for f in "$CFDNA_DIR"/*.bed.gz; do
    echo "Filtering R1 reads in $f"
    zcat "$f" | awk 'NR==1 || and($9,64)' | gzip > "${f}.tmp" && mv "${f}.tmp" "$f"
done

# Filter tissue files to only R1
# NR=1: ignore the header, col 9 contains flags
for f in "$TISSUE_DIR"/*.bed.gz; do
    echo "Filtering R1 reads in $f"
    zcat "$f" | awk 'NR==1 || and($9,64)' | gzip > "${f}.tmp" && mv "${f}.tmp" "$f"
done

cd "$WORKDIR" || { echo "Error: Cannot change to working directory $WORKDIR"; exit 1; }

python Generate_samples.py \
    --cfDNA_dir "$CFDNA_DIR" \
    --tissue_dir "$TISSUE_DIR" \
    --output_dir "$OUTDIR"


# Example usage:
# sbatch -p long 01-create_samples.sh -w /users/ludwig/cnr137 -c /well/ludwig/users/cnr137/methylation_model/healthy_cfdna_samples -t /well/ludwig/users/cnr137/methylation_model/tissue_samples -o /well/ludwig/users/cnr137/methylation_model/generated_samples
