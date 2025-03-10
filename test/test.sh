#!/bin/bash

# Define the output file
OUTPUT_FILE="results.txt"
> "$OUTPUT_FILE"  # Clear file before writing

# Define range
BASE_M=1028       # Starting value for M
EXP_MAX_M=2   # M will be BASE_M * 2^i for i in [0, EXP_MAX_M]

BASE_N=8       # Starting value for N
EXP_MAX_N=4   # N will be BASE_N * 2^i for i in [0, EXP_MAX_N]

# Repeat (M, N) if below this time up to this number of repetitions
MIN_TIME=30
MAX_EXECUTIONS=5

# Loop over exponentially growing values of M and N
for ((i=0; i<=EXP_MAX_M; i++)); do
    M=$(( BASE_M * 2**i ))
    
    for ((j=0; j<=EXP_MAX_N; j++)); do
        N=$(( BASE_N * 2**j ))
        
        SUM_TIME=0
        COUNT=0  # Number of runs

        while (( $(echo "$SUM_TIME < $MIN_TIME" | bc -l) && COUNT < MAX_EXECUTIONS )); do
            OUTPUT=$(./main expm-time-coo ./test/large/laplacian_3d_262144.csv "$M" "$N" -1 10 10)
            RESULT=$(echo "$OUTPUT" | sed -n '1p')
            TIME=$(echo "$OUTPUT" | sed -n '2p' | awk '{print $3}')  # Extract execution time
            
            SUM_TIME=$(echo "$SUM_TIME + $TIME" | bc)
            ((COUNT++))
            echo $SUM_TIME
        done
        
        AVG_TIME=$(echo "$SUM_TIME / $COUNT" | bc -l)

        # Format the result with a leading zero and the specified decimal places
        AVG_TIME=$(printf "%0.6f" "$AVG_TIME")
        
        echo "$M $N $RESULT $AVG_TIME"
        echo "$M $N $RESULT $AVG_TIME" >> "$OUTPUT_FILE"
    done
done
