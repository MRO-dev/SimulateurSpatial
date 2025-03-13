import os

# Define the list of prefixes and suffixes
PREFIXES = ['BLUE1','BLUE2','RED1','RED2']
SUFFIXES = [
    'AEM_ATTITUDE.TXT',
    'COLOR.TXT',
    'MEM_LLA.TXT',
    'OEM_POSITION.TXT'
]

# Generate the set of target filenames to delete
TARGET_FILES = set()
for prefix in PREFIXES:
    for suffix in SUFFIXES:
        TARGET_FILES.add("{}_{}".format(prefix, suffix).upper())

def delete_specific_files(directory):
    """
    Recursively delete specific files in the given directory and its subdirectories,
    excluding any directories named 'Backup'.

    Parameters:
    - directory (str): The path to the directory to process.
    """
    print("MEMxterminator called from {}".format(directory))

    try:
        for fileName in os.listdir(directory):
            filePath = os.path.join(directory, fileName)

            # Skip the 'Backup' directory
            if fileName == "Backup" and os.path.isdir(filePath):
                print("Skipping Backup directory: {}".format(filePath))
                continue  # Ignore the Backup folder

            # Check if the current item is a file and matches one of the target filenames
            if os.path.isfile(filePath) and fileName.upper() in TARGET_FILES:
                try:
                    os.remove(filePath)
                    print("{} was deleted successfully.".format(fileName))
                except PermissionError:
                    print("Permission denied: Could not delete {}.".format(fileName))
                except FileNotFoundError:
                    print("File not found: {} does not exist.".format(fileName))
                except Exception as e:
                    print("Error deleting {}: {}".format(fileName, e))

            # If the current item is a directory (and not 'Backup'), recurse into it
            elif os.path.isdir(filePath):
                try:
                    delete_specific_files(filePath)
                except Exception as e:
                    print("Error accessing directory {}: {}".format(filePath, e))

    except PermissionError:
        print("Permission denied: Cannot access directory {}.".format(directory))
    except FileNotFoundError:
        print("Directory not found: {} does not exist.".format(directory))
    except Exception as e:
        print("Error accessing directory {}: {}".format(directory, e))

if __name__ == "__main__":
    currentWorkingDirectory = os.getcwd()
    delete_specific_files(currentWorkingDirectory)

