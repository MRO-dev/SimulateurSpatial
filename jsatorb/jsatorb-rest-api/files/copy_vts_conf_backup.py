import os
import shutil
import sys

def copy_parent_vts(target_dir):
    # Get the parent directory (one level up)
    parent_dir = os.path.dirname(target_dir)
    
    # Get the list of files in the parent directory
    parent_files = os.listdir(parent_dir)
    
    # Copy each .vts file from parent to current directory
    for file in parent_files:
        if file.endswith(".vts"):
            source_path = os.path.join(parent_dir, file)
            target_path = os.path.join(target_dir, file)
            
            if os.path.exists(target_path):
                print("Overwriting existing file: {}".format(target_path))
            shutil.copy2(source_path, target_path)
            print("Copied {} to {}".format(source_path, target_path))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        target_directory = os.getcwd()  # Use the current working directory
    else:
        target_directory = sys.argv[1]
    copy_parent_vts(target_directory)
