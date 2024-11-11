#!/bin/bash

# Set base directory
BASE_DIR="/BASE_PATH"

# Set project directory
PRJ_DIR="${BASE_DIR}/PATH_TO_DATA/HCC-sp-RNAseq/cytospace"

# Set paths for cytospace input and output
INPUT_DIR="${PRJ_DIR}/input"
OUTPUT_DIR_RELATIVE="output" # Relative path

# Ensure the output directory exists
mkdir -p "${OUTPUT_DIR_RELATIVE}"

# Argument passed to the script
sample=$1

# Construct paths for current sample
scRNA_path="${INPUT_DIR}/SC/scRNA_scRNA_data.txt"
cell_type_path="${INPUT_DIR}/SC/scRNA_cell_type_labels.txt"
st_path="${INPUT_DIR}/ST/${sample}_ST_data.txt"
coordinates_path="${INPUT_DIR}/ST/${sample}_Coordinates.txt"

# Activate conda environment (if you use conda for cytospace)
echo "Activating conda environment 'cytospace'..."
source activate cytospace

# Run cytospace
cd "${PRJ_DIR}" # Ensure we are in the right working directory
cytospace \
   --st-path "${st_path}" \
   --coordinates-path "${coordinates_path}" \
   --scRNA-path "${scRNA_path}" \
   --cell-type-path "${cell_type_path}" \
   --output-folder "${OUTPUT_DIR_RELATIVE}" \
   --number-of-processors 96 \
   --output-prefix "${sample}_" || exit 1

echo "Cytospace processing completed for ${sample}."