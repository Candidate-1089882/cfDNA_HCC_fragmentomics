#!/bin/bash
#SBATCH --mem=300G

module load R/4.3.2-gfbf-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2

cd /well/ludwig/users/cnr137/references

# Run the R script
Rscript /users/ludwig/cnr137/generate_cytosine_ref_Hg19.R