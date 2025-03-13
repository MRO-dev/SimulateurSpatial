# Anthony Le Batteux EAE 130723
import os
import shutil
import sys

def backup_files(source_dir):
    print("TEST")
    # Create the path for the backup directory
    backup_dir = os.path.join(source_dir, "Backup")
    
    # If the backup folder exists, remove it
    if os.path.exists(backup_dir):
        shutil.rmtree(backup_dir)
    
    # Recreate the (empty) backup folder
    os.makedirs(backup_dir)
    
    # List files in the source directory (top level only)
    files = os.listdir(source_dir)

    # Copy each .TXT file from the top level into Backup
    for file in files:
        source_path = os.path.join(source_dir, file)
        backup_path = os.path.join(backup_dir, file)
        # Check if it is a file and ends with .TXT
        if os.path.isfile(source_path) and file.endswith(".TXT"):
            shutil.copy2(source_path, backup_path)
            print("Copied {} to {}".format(source_path, backup_path))

    # We intentionally do NOT copy the "Manoeuver" directory anymore.
    # That ensures .TXT files inside the original Manoeuver folder are not copied.

print("TEST2")

currentWorkingDirectory = os.getcwd()
backup_files(currentWorkingDirectory)
backup_dir = os.path.join(currentWorkingDirectory, "Backup")

# Directory containing the 'Result...' and '...ManeuverDate...' files
result_file_dir = "/app/maneuver-manager"

# Now we copy these special files directly into Backup:
def copy_if_exists(filename):
    """
    Helper function to check if a file exists in result_file_dir
    and copy it to the backup directory.
    """
    result_file = os.path.join(result_file_dir, filename)
    if os.path.exists(result_file):
        shutil.copy2(result_file, backup_dir)
        print("Copied '{}' from {} to {}".format(filename, result_file, backup_dir))


# List of result/date filenames you want to copy to Backup
files_to_copy = [
    "Result.txt",
    "Result2.txt",
    "Result3.txt",
    "Result4.txt",
    "LastManeuverDate.txt",
    "LastManeuverDate2.txt",
    "LastManeuverDate3.txt",
    "LastManeuverDate4.txt",
    "PostManeuverDate.txt",
    "PostManeuverDate2.txt",
    "PostManeuverDate3.txt",
    "PostManeuverDate4.txt",
    "blue1-tle-storage.json",
    "blue1-tle-mongo.json",
    "blue2-tle-storage.json",
    "blue2-tle-mongo.json",
    "red1-tle-storage.json",
    "red1-tle-mongo.json",
    "red2-tle-storage.json",
    "red2-tle-mongo.json",
    "time-persistence.json"
]

for f in files_to_copy:
    copy_if_exists(f)

