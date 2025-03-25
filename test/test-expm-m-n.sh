#!/bin/bash

# Output file, append by default
OUTPUT_FILE="test/res/test-expm-m-n.csv"
if [[ $1 == "-t" ]]
then
    > "$OUTPUT_FILE"  # Truncate file
    echo "M,N,ERROR" >> "$OUTPUT_FILE"
else
    echo "" >> "$OUTPUT_FILE"
fi

# Fixed values
SEED=-1
ROW=0
COL=1

# Ranges
BASE_M=65536
EXP_MAX_M=8

BASE_N=256
EXP_MAX_N=8

# Repeat (M, N) this number of times
REPETITIONS=10

# Run REPETITIONS time and write results
run_experiment() {
    local M=$1
    local N=$2
    for ((k=0; k<REPETITIONS; k++)); do
        OUTPUT=$(./main exp ./test/laplacians/laplacian_3d_32768.csv \
            -r ./test/laplacians/laplacian_3d_32768_exp.csv \
            "$M" "$N" "$ROW" "$COL" "$SEED")
        ERROR=$(echo "$OUTPUT" | awk '/Error:/ {print $2}')
        
        echo "$M,$N,$ERROR"
        echo "$M,$N,$ERROR" >> "$OUTPUT_FILE"
    done
}

# Loop over exponentially growing values of M
N=4096
for ((i=0; i<=EXP_MAX_M; i++))
do
    M=$(( BASE_M * 2**i ))
    run_experiment "$M" "$N"
done

# Loop over exponentially growing values of N
M=1048576
for ((i=0; i<=EXP_MAX_N; i++))
do
    N=$(( BASE_N * 2**i ))
    run_experiment "$M" "$N"
done
