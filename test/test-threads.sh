#!/bin/bash

# Define the output file
OUTPUT_FILE="results.txt"
> "$OUTPUT_FILE"  # Clear file before writing

# Fixed values for M and N (modify these as needed)
M=1000000
N=100

# Maximum number of threads to test
MAX_THREADS=$(nproc)

# Loop through increasing number of threads
THREADS=1
while [ $THREADS -le $MAX_THREADS ]
do
    # Set the OMP_NUM_THREADS environment variable to the current thread count
    export OMP_NUM_THREADS=$THREADS

    # Call your program with the fixed values M and N and the current number of threads
    echo "Running with $THREADS threads..."
    OUTPUT=$(./main expm-time-coo ./test/large/laplacian_3d_262144.csv "$M" "$N" -1 1 1)
    TIME=$(echo "$OUTPUT" | sed -n '2p' | awk '{print $3}')  # Extract execution time

    echo "$M $N $TIME"
    echo "$M $N $TIME" >> "$OUTPUT_FILE"

    THREADS=$((THREADS * 2))
done
