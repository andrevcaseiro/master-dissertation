#!/bin/bash

# Output file, truncated at start
OUTPUT_FILE="test-expm-m-n.csv"
> "$OUTPUT_FILE"

# Fixed values
SEED=-1
ROW=1
COL=1

# Ranges
BASE_M=1028
EXP_MAX_M=2

BASE_N=8
EXP_MAX_N=4

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
            OUTPUT=$(./main expm-time-coo ./test/large/laplacian_3d_262144.csv "$M" "$N" "$SEED" "$ROW" "$COL")
            RESULT=$(echo "$OUTPUT" | sed -n '1p')
            TIME=$(echo "$OUTPUT" | sed -n '2p' | awk '{print $3}')  # Extract execution time
            
            echo "$M $N $RESULT $AVG_TIME"
            echo "$M $N $RESULT $AVG_TIME" >> "$OUTPUT_FILE"
        done
    done
done
