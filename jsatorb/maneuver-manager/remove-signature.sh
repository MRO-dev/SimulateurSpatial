#!/bin/bash

# Print error message and exit
error_exit() {
    echo "Error: $1" >&2
    exit 1
}

# Check if required commands are available
command -v unzip >/dev/null 2>&1 || error_exit "unzip command not found. Please install unzip."
command -v jar >/dev/null 2>&1 || error_exit "jar command not found. Please install Java Development Kit."

# Check if a JAR file was provided as an argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <jar_file>"
    exit 1
fi

jar_file="$1"

# Check if the JAR file exists
if [ ! -f "$jar_file" ]; then
    error_exit "File '$jar_file' not found!"
fi

# Check if the file is actually a JAR file
file_type=$(file "$jar_file")
if ! echo "$file_type" | grep -q "Zip archive data\|Java archive data"; then
    error_exit "'$jar_file' does not appear to be a valid JAR file!"
fi

# Create a temporary directory
temp_dir=$(mktemp -d)
echo "Created temporary directory: $temp_dir"

# Extract the JAR file
echo "Extracting JAR file..."
unzip -q "$jar_file" -d "$temp_dir" || error_exit "Failed to extract JAR file"

# Check if META-INF directory exists
if [ ! -d "$temp_dir/META-INF" ]; then
    error_exit "META-INF directory not found in JAR file"
fi

# Remove signature files from META-INF
echo "Removing signature files..."
signature_files=$(find "$temp_dir/META-INF" -type f \( \
    -name "*.SF" -o \
    -name "*.RSA" -o \
    -name "*.DSA" -o \
    -name "*.EC" -o \
    -name "ECLIPSE_.*" -o \
    -name "SIG-*" -o \
    -name "*.SIGN" \
    \) -print)

if [ -n "$signature_files" ]; then
    echo "Found signature files to remove:"
    echo "$signature_files"
    find "$temp_dir/META-INF" -type f \( \
        -name "*.SF" -o \
        -name "*.RSA" -o \
        -name "*.DSA" -o \
        -name "*.EC" -o \
        -name "ECLIPSE_.*" -o \
        -name "SIG-*" -o \
        -name "*.SIGN" \
        \) -delete
    echo "Signature files removed successfully"
else
    echo "No signature files found"
fi

# Clean up MANIFEST.MF
if [ -f "$temp_dir/META-INF/MANIFEST.MF" ]; then
    echo "Cleaning MANIFEST.MF..."
    # Create a temporary file for the new manifest
    temp_manifest=$(mktemp)
    # Remove signature-related entries from MANIFEST.MF
    grep -v -E "^(Sealed:|Signature-Version:|Created-By:|SHA1-Digest:|SHA-256-Digest:|MD5-Digest:|Sign:|Certificate:|LastModified:|BundleSign:|SHA-Digest:|DSA-Signature:|RSA-Signature:)" "$temp_dir/META-INF/MANIFEST.MF" > "$temp_manifest"
    mv "$temp_manifest" "$temp_dir/META-INF/MANIFEST.MF"
fi

# Create backup of original JAR
backup_file="${jar_file}.backup"
echo "Creating backup of original JAR as: $backup_file"
cp "$jar_file" "$backup_file" || error_exit "Failed to create backup file"

# Repack the JAR file
echo "Repacking JAR file..."
original_dir=$(pwd)
cd "$temp_dir" || error_exit "Failed to change to temporary directory"
jar -cfM "$original_dir/${jar_file}" ./* || error_exit "Failed to create new JAR file"
cd "$original_dir" || error_exit "Failed to return to original directory"

# Clean up
echo "Cleaning up temporary files..."
rm -rf "$temp_dir"

echo "Process completed successfully!"
echo "Original JAR backed up as: $backup_file"
echo "Signature-free JAR created as: $jar_file"
