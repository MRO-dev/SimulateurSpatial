import os
import shutil
import sys

def backup_files(source_dir):
    parent_dir = os.path.dirname(source_dir)
    files = os.listdir(source_dir)
    
    for file in files:
        if file.endswith(".vts"):
            source_path = os.path.join(source_dir, file)
            backup_path = os.path.join(parent_dir, file)
            print(source_path)
            print(backup_path)
            
            if os.path.exists(backup_path):
                print("Overwriting existing file: {}".format(backup_path))
            shutil.copy2(source_path, backup_path)
            print("Copied {} to {}".format(source_path, backup_path))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        source_directory = os.getcwd()
    else:
        source_directory = sys.argv[1]
    backup_files(source_directory)
