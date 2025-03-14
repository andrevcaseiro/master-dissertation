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
BASE_M=263168
EXP_MAX_M=0

BASE_N=16
EXP_MAX_N=10

# Repeat (M, N) this number of times
REPETITIONS=10

# Loop over exponentially growing values of M and N
for ((i=0; i<=EXP_MAX_M; i++))
do
    M=$(( BASE_M * 2**i ))
    
    for ((j=0; j<=EXP_MAX_N; j++))
    do
        N=$(( BASE_N * 2**j ))

        for ((k=0; k<=REPETITIONS; k++))
        do
            OUTPUT=$(./main expm-test-coo ./test/laplacians/laplacian_3d_32768.csv ./test/laplacians/laplacian_3d_32768_exp.csv "$M" "$N" "$SEED" "$ROW" "$COL")
            ERROR=$(echo "$OUTPUT" | awk '{print $1}' | tr -d '%')
            
            echo "$M,$N,$ERROR"
            echo "$M,$N,$ERROR" >> "$OUTPUT_FILE"
        done
    done
done
