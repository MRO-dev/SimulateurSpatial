#!/bin/bash

# Default file name is LastManeuverDate.txt (for "blue1" mode)
file="LastManeuverDate.txt"

# If a mode parameter was provided, set the file name accordingly
if [ "$#" -ge 1 ]; then
    case "$1" in
        blue1)
            file="LastManeuverDate.txt"
            ;;
        blue2)
            file="LastManeuverDate2.txt"
            ;;
        red1)
            file="LastManeuverDate3.txt"
            ;;
        red2)
            file="LastManeuverDate4.txt"
            ;;
        *)
            echo "Unknown mode parameter: '$1'"
            exit 1
            ;;
    esac
fi

# Check if the file exists
if [ ! -f "$file" ]; then
    echo "File '$file' does not exist."
    exit 1
fi

# Count how many lines the file has
line_count=$(wc -l < "$file")

# Ensure there are at least 2 lines
if [ "$line_count" -lt 2 ]; then
    echo "File '$file' has fewer than two lines. Cannot extract last two lines."
    exit 1
fi

# Extract the last two lines into a temporary file
tail -n 2 "$file" > "${file}.tmp"

# Replace the original file with the temporary file
mv "${file}.tmp" "$file"

echo "Extracted the last two lines from '$file'."

