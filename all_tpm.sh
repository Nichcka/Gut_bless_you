#!/bin/bash

DIR=$1
OUTPUT_FILE=$2

> "$OUTPUT_FILE"

for dir in "$DIR"/*/; do
    quant_file="${dir}quant.sf"
    folder_name=$(basename "$dir")

    if [[ -f "$quant_file" ]]; then

        tail -n +2 "$quant_file" | awk -v sample="$folder_name" '$4 != 0 {print sample "\t" $1 "\t" $4}' >> "$OUTPUT_FILE"

    fi
done

