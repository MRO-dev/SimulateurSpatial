import os
import shutil

def copy_files_from_manoeuver():
    # Get the current working directory
    current_dir = os.getcwd()

    # Define the source directory path (Manoeuver folder inside current directory)
    source_dir = os.path.join(current_dir, "Manoeuver")

    # Define the destination directory path (current directory)
    destination_dir = current_dir

    # Check if the source directory exists
    if not os.path.exists(source_dir):
        print("Error: Source directory '{}' does not exist.".format(source_dir))
        return

    # Iterate over all items in the source directory
    for item in os.listdir(source_dir):
        source_item = os.path.join(source_dir, item)
        destination_item = os.path.join(destination_dir, item)

        # Check if it's a file (ignore subdirectories)
        if os.path.isfile(source_item):
            # Check if the file already exists in the destination
            if os.path.exists(destination_item):
                print("Skipped (already exists): {}".format(item))
                continue  # Skip copying this file

            try:
                shutil.copy2(source_item, destination_item)
                print("Copied file: {}".format(item))
            except Exception as e:
                print("Failed to copy file {}. Reason: {}".format(item, e))
        else:
            print("Skipped (not a file): {}".format(item))

if __name__ == "__main__":
    copy_files_from_manoeuver()

