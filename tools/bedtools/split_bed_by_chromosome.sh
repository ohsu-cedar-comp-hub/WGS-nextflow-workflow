#!/bin/bash
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <bed_file>"
    exit 1
fi

input_file=$1

# Loop through each line in the file
while IFS= read -r line; do
    # Extract the chromosome name, which is the first field in a BED file
    chr=$(echo "$line" | cut -f 1)
    
    # Check if the chromosome is one of the expected chromosomes
    if [[ "$chr" =~ ^(chr[1-9]|chr1[0-9]|chr2[0-2]|chrX)$ ]]; then
        # Write the line to the corresponding chromosome file
        echo "$line" >> "${chr}.bed"
    else
        echo "Skipping line with unexpected chromosome: $chr"
    fi
done < "$input_file"

