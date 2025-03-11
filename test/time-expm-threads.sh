#!/bin/bash

# Output file, truncated at start
OUTPUT_FILE="results.txt"
> "$OUTPUT_FILE"

# Fixed values
M=16777216
N=512
SEED=1234
ROW=1
COL=1

# Header
echo "THREADS,M,N,SEED,ROW,COL,RESULT,TIME" >> "$OUTPUT_FILE"

# Start threads, doubled each iteration until the max
THREADS=1
MAX_THREADS=$(nproc)
while [ $THREADS -le $MAX_THREADS ]
do
    # Set the OMP_NUM_THREADS environment variable to the current thread count
    export OMP_NUM_THREADS=$THREADS

    # Execute test and parse output
    OUTPUT=$(./main expm-time-coo ./test/large/laplacian_3d_262144.csv "$M" "$N" "$SEED" "$ROW" "$COL")
    RESULT=$(echo "$OUTPUT" | sed -n '1p')
    TIME=$(echo "$OUTPUT" | sed -n '2p' | awk '{print $3}')  # Extract execution time

    # Output to terminal and file
    echo "$THREADS,$M,$N,$SEED,$ROW,$COL,$RESULT,$TIME"
    echo "$THREADS,$M,$N,$SEED,$ROW,$COL,$RESULT,$TIME" >> "$OUTPUT_FILE"

    THREADS=$((THREADS * 2))
done
