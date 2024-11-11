#!/bin/bash

# Set base directory
BASE_DIR="/BASE_PATH"

# Set project directory
PRJ_DIR="${BASE_DIR}/PATH_TO_DATA/HCC-sp-RNAseq/cytospace"

# Set paths for cytospace input and output
INPUT_DIR="${PRJ_DIR}/cytospace/input/ST"
OUTPUT_DIR="${PRJ_DIR}/cytospace/output"

# Ensure the output directory exists
mkdir -p "${OUTPUT_DIR}"

# Option 1: Automatically detect samples based on the available ST_data files in the directory
samples=($(ls "${INPUT_DIR}" | grep "_ST_data.txt" | sed 's/_ST_data.txt//'))

# Option 2: Manually specify the samples (Uncomment the line below and provide the sample names)
samples=("cHC-1T" "HCC-1T" "HCC-2T" "HCC-3T" "HCC-4T" "HCC-5A" "HCC-5B" "HCC-5C" "HCC-5D") 

# Run cytospace for each sample using parallel
echo "Running cytospace for each sample..."

# Uncomment the following section if you want to submit jobs sequentially. 
# Each job will wait for the previous one to finish before starting.
for sample in "${samples[@]}"; do
    bash run_cytospace.sh "${sample}" > "${OUTPUT_DIR}/${sample}_$(date +%Y%m%d).log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error encountered while processing sample ${sample}. Check the log file for details."
        exit 1
    fi
done

# Uncomment the following section if you want to submit all jobs in the background.
# Each job will start after a specified delay (in this case, 1000 seconds).
# This approach allows for some degree of parallel processing, but with controlled submission to avoid overloading the system.
# for sample in "${samples[@]}"; do
#     nohup bash run_cytospace.sh "${sample}" > "${OUTPUT_DIR}/${sample}_$(date +%Y%m%d).log" 2>&1 &
#     sleep 1000 # This will add a 1000-second delay between job submissions.
# done

echo "All cytospace jobs submitted."