#!/bin/bash

# Output CSV header
echo "program,threads,steps,size,time(s)" > results.csv

for file in *.o[0-9]*; do
    base=$(basename "$file")
    IFS='_' read -r program threads steps size_ext <<< "$base"
    size=$(echo "$size_ext" | cut -d'.' -f1)

    # Extract time from inside the file
    time=$(grep -oP 'Execution Time\s*=\s*\K[0-9.]+' "$file")

    echo "$program,$threads,$steps,$size,$time" >> results.csv
done

echo "Results CSV generated @ results.csv"
