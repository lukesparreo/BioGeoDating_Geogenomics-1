#!/bin/bash

# Paths to tree files
FILE1="/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/informed_uniform/output_informed_uniform1.tre"
FILE2="/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/informed_uniform/output_informed_uniform2.tre"
OUTPUT="/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/informed_uniform/output_informed_uniform_migration0_combined.tre"

# Temporary files
TMP1=$(mktemp)
TMP2=$(mktemp)

# Get number of lines (excluding header)
LINES1=$(tail -n +2 "$FILE1" | wc -l)
LINES2=$(tail -n +2 "$FILE2" | wc -l)

# Calculate number of lines to keep (after 25% burn-in)
START1=$(( (LINES1 * 25 / 100) + 2 ))
START2=$(( (LINES2 * 25 / 100) + 2 ))

# Write header from first file
head -n 1 "$FILE1" > "$OUTPUT"

# Append post-burnin trees
tail -n +"$START1" "$FILE1" >> "$OUTPUT"
tail -n +"$START2" "$FILE2" >> "$OUTPUT"

# Clean up
rm "$TMP1" "$TMP2"

echo "Combined tree file written to $OUTPUT"

