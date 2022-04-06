#!/usr/bin/env bash

# This script grabs the accession ids from the first column of the SRA runtable and writes it to a new .txt file.

input_file="$1"
echo "Processing Input File: $input_file"
output_file="../accessions/accessions_ace2covid86.txt"

field_num=1

echo

echo "Output Filename: $output_file"
echo "Field Number: $field_num"

cat $input_file | cut -f$field_num -d','  | sed '1d' > $output_file

echo "Accessions have been extracted from the SRA runtable and written to: "$output_file
