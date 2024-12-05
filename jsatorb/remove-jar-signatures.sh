#!/bin/bash

# Check if a JAR file was provided as an argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <jar_file>"
    exit 1
fi

jar_file="$1"

# Check if the JAR file exists
if [ ! -f "$jar_file" ]; then
    echo "Error: File '$jar_file' not found!"
    exit 1
fi

# Create a temporary directory
temp_dir=$(mktemp -d)
echo "Created temporary directory: $temp_dir"

# Extract the JAR file
echo "Extracting JAR file..."
unzip "$jar_file" -d "$temp_dir"

# Remove signature files from META-INF
echo "Removing signature files..."
find "$temp_dir/META-INF" -type f \( -name "*.SF" -o -name "*.RSA" -o -name "*.DSA" \) -delete

# Create backup of original JAR
backup_file="${jar_file}.backup"
echo "Creating backup of original JAR as: $backup_file"
cp "$jar_file" "$backup_file"

# Repack the JAR file
echo "Repacking JAR file..."
cd "$temp_dir"
jar -cfM "../${jar_file}" ./*

# Clean up
echo "Cleaning up temporary files..."
cd ..
rm -rf "$temp_dir"

echo "Process completed successfully!"
echo "Original JAR backed up as: $backup_file"
echo "Signature-free JAR created as: $jar_file"
