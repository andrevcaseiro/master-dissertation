#!/bin/bash

# Output file, append by default
OUTPUT_FILE="test/res/time-expm-m-n.csv"
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
BASE_M=4096
EXP_MAX_M=16

BASE_N=16
EXP_MAX_N=16

# Repeat (M, N) if below this time up to this number of repetitions
MIN_TIME=30
MAX_EXECUTIONS=5

# Execute on a single thread
THREADS=1
export OMP_NUM_THREADS=$THREADS

# Run and write results
run_experiment() {
    local M=$1
    local N=$2

    SUM_TIME=0
    COUNT=0  # Number of runs

    while (( $(echo "$SUM_TIME < $MIN_TIME" | bc -l) && COUNT < MAX_EXECUTIONS ))
    do
        # Execute test and parse output
        OUTPUT=$(./main exp ./test/large/laplacian_3d_262144.csv \
            "$M" "$N" "$ROW" "$COL" "$SEED")
        RESULT=$(echo "$OUTPUT" | awk '/Result:/ {print $2}')
        TIME=$(echo "$OUTPUT" | awk '/Time:/ {print $2}')
        
        SUM_TIME=$(echo "$SUM_TIME + $TIME" | bc)
        ((COUNT++))
    done
    
    # Calculate average and format number of zeros
    AVG_TIME=$(echo "$SUM_TIME / $COUNT" | bc -l)
    AVG_TIME=$(printf "%0.5f" "$AVG_TIME")
    
    # Output to terminal and file
    echo "$THREADS,$M,$N,$SEED,$ROW,$COL,$RESULT,$AVG_TIME"
    echo "$THREADS,$M,$N,$SEED,$ROW,$COL,$RESULT,$AVG_TIME" >> "$OUTPUT_FILE"
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
