#!/bin/bash

file="LastManeuverDate.txt"

# Check if the file has at least two lines
if [ $(wc -l < "$file") -lt 2 ]; then
    echo "The file has fewer than two lines. Cannot remove last two lines."
    exit 1
fi

# Remove the last two lines by extracting all lines except the last two
head -n -2 "$file" > "${file}.tmp"

# Replace the original file with the temporary file
mv "${file}.tmp" "$file"

echo "Removed the last two lines from $file."

