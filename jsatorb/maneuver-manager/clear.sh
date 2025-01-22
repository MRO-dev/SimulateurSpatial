#!/bin/bash

file="LastManeuverDate.txt"

# Extract last two lines into a temporary file
tail -n 2 "$file" > "${file}.tmp"

# Replace the original file with the temporary file
mv "${file}.tmp" "$file"
