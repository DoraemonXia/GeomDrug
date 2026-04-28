#!/bin/bash

# Extract sequence IDs from a FASTA file
# Check if a file path was provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <path_to_fasta_file>"
  exit 1
fi

fasta_file=$1

# Use awk to extract sequence IDs from the FASTA file
sequence_ids=$(awk '/^>/ {print substr($0, 2)}' "$fasta_file")

# Store sequence IDs into an array
sequence_id_array=($sequence_ids)

# Output sequence IDs
for id in "${sequence_id_array[@]}"; do
  echo "$id"
done