#!/bin/bash

# Output file, append by default
OUTPUT_FILE="time-expm-m-n.csv"
if [[ $1 == "-t" ]]
then
    > "$OUTPUT_FILE"  # Truncate file
    echo "THREADS,M,N,SEED,ROW,COL,RES,TIME" >> "$OUTPUT_FILE"
else
    echo "" >> "$OUTPUT_FILE"
fi

# Fixed values
SEED=1234
ROW=1
COL=1

# Ranges
BASE_M=8224
EXP_MAX_M=1

BASE_N=8
EXP_MAX_N=10

# Repeat (M, N) if below this time up to this number of repetitions
MIN_TIME=30
MAX_EXECUTIONS=5

# Execute on a single thread
THREADS=1
export OMP_NUM_THREADS=$THREADS

# Loop over exponentially growing values of M and N
for ((i=0; i<=EXP_MAX_M; i++))
do
    M=$(( BASE_M * 2**i ))
    
    for ((j=0; j<=EXP_MAX_N; j++))
    do
        N=$(( BASE_N * 2**j ))
        
        SUM_TIME=0
        COUNT=0  # Number of runs

        while (( $(echo "$SUM_TIME < $MIN_TIME" | bc -l) && COUNT < MAX_EXECUTIONS ))
        do
            # Execute test and parse output
            OUTPUT=$(./main expm-time-coo ./test/large/laplacian_3d_262144.csv "$M" "$N" "$SEED" "$ROW" "$COL")
            RESULT=$(echo "$OUTPUT" | sed -n '1p')
            TIME=$(echo "$OUTPUT" | sed -n '2p' | awk '{print $3}')  # Extract execution time
            
            SUM_TIME=$(echo "$SUM_TIME + $TIME" | bc)
            ((COUNT++))
        done
        
        # Calculate average and format number of zeros
        AVG_TIME=$(echo "$SUM_TIME / $COUNT" | bc -l)
        AVG_TIME=$(printf "%0.6f" "$AVG_TIME")
        
        # Output to terminal and file
        echo "$THREADS,$M,$N,$SEED,$ROW,$COL,$RESULT,$AVG_TIME"
        echo "$THREADS,$M,$N,$SEED,$ROW,$COL,$RESULT,$AVG_TIME" >> "$OUTPUT_FILE"
    done
done
