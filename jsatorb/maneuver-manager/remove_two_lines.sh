#!/bin/bash

# Default file name
file="LastManeuverDate.txt"

# Check for a mode parameter
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
            echo "Unknown mode parameter: $1"
            exit 1
            ;;
    esac
fi

# Check if the file exists
if [ ! -f "$file" ]; then
    echo "File '$file' does not exist."
    exit 1
fi

# Count the file’s lines
line_count=$(wc -l < "$file")

# Check if the file has at least two lines
if [ "$line_count" -lt 2 ]; then
    echo "The file '$file' has fewer than two lines. Cannot remove last two lines."
    exit 1
fi

# Remove the last two lines by extracting all lines except the last two
head -n -2 "$file" > "${file}.tmp" || exit 1

# Replace the original file with the temporary file
mv "${file}.tmp" "$file"

echo "Removed the last two lines from $file."

