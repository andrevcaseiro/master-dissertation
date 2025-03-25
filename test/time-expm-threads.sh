#!/bin/bash

# Output file, append by default
OUTPUT_FILE="test/res/time-expm-threads.csv"
if [[ $1 == "-t" ]]
then
    > "$OUTPUT_FILE"  # Truncate file
    echo "M N RESULT ERROR" >> "$OUTPUT_FILE"
else
    echo "" >> "$OUTPUT_FILE"
fi

# Fixed values
M=1048576
N=4096
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
    # Execute test and parse output
    OUTPUT=$(./main exp ./test/large/laplacian_3d_262144.csv \
        "$M" "$N" "$ROW" "$COL" "$SEED")
    RESULT=$(echo "$OUTPUT" | awk '/Result:/ {print $2}')
    TIME=$(echo "$OUTPUT" | awk '/Time:/ {print $2}')

    # Output to terminal and file
    echo "$THREADS,$M,$N,$SEED,$ROW,$COL,$RESULT,$TIME"
    echo "$THREADS,$M,$N,$SEED,$ROW,$COL,$RESULT,$TIME" >> "$OUTPUT_FILE"

    THREADS=$((THREADS * 2))
done
