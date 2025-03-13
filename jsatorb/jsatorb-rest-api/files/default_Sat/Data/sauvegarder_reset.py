#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import time
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def ensure_directory_removed(directory):
    """
    Ensures a directory is removed, with multiple retries if needed.
    Returns True if successful, False otherwise.
    """
    max_retries = 3
    retry_delay = 1.0  # seconds
    
    for attempt in range(max_retries):
        try:
            if os.path.exists(directory):
                shutil.rmtree(directory)
            return True
        except Exception as e:
            logger.warning("Attempt {} failed to remove {}: {}".format(attempt + 1, directory, e))
            if attempt < max_retries - 1:
                time.sleep(retry_delay)
    return False

def safe_copy(src, dst):
    """
    Safely copy a file with error handling and logging.
    """
    try:
        shutil.copy2(src, dst)
        logger.info("Copied {} => {}".format(src, dst))
        return True
    except Exception as e:
        logger.error("Failed to copy {} to {}: {}".format(src, dst, e))
        return False

def backup_files(source_dir):
    """
    Copies .TXT files from `source_dir` into `source_dir/Reset/`,
    and copies the entire `Manoeuver/` folder (if it exists) inside `source_dir`
    into `source_dir/Reset/Manoeuver/`.
    """
    backup_dir = os.path.join(source_dir, "Reset")
    
    # 1. Remove existing backup directory with retry mechanism
    if not ensure_directory_removed(backup_dir):
        raise RuntimeError("Failed to remove existing backup directory: {}".format(backup_dir))
    
    # 2. Create new backup directory
    try:
        os.makedirs(backup_dir)
        logger.info("Created new backup folder: {}".format(backup_dir))
    except Exception as e:
        raise RuntimeError("Failed to create backup directory {}: {}".format(backup_dir, e))

    # 3. Copy all .TXT files
    for file_name in os.listdir(source_dir):
        source_path = os.path.join(source_dir, file_name)
        if os.path.isfile(source_path) and file_name.endswith(".TXT"):
            destination_path = os.path.join(backup_dir, file_name)
            safe_copy(source_path, destination_path)

    # 4. Copy Maneuver directory if it exists
    manoeuver_dir = os.path.join(source_dir, "Manoeuver")
    if os.path.isdir(manoeuver_dir):
        backup_manoeuver_dir = os.path.join(backup_dir, "Manoeuver")
        try:
            shutil.copytree(manoeuver_dir, backup_manoeuver_dir)
            logger.info("Copied entire 'Manoeuver' directory => {}".format(backup_manoeuver_dir))
        except Exception as e:
            raise RuntimeError("Failed to copy Manoeuver directory: {}".format(e))

def modify_last_maneuver_date(file_path):
    """
    Removes the last two lines from the file at the given path, 
    ensuring no residual newline characters are left.
    """
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        
        # Remove last two lines if they exist
        if len(lines) >= 2:
            modified_lines = lines[:-2]
        else:
            modified_lines = lines
            
        # Write modified content
        with open(file_path, 'w') as file:
            file.write(''.join(modified_lines))
        logger.info("Removed last two lines from {}".format(file_path))
    except Exception as e:
        raise RuntimeError("Error modifying {}: {}".format(file_path, e))
        
def remove_last_couple_of_lines(file_path):
    """
    Removes the last 'couple' of lines from the file.
    In a scenario where each 'couple' is effectively 3 lines
    (including a blank line), this removes the last 3 lines.
    """
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Remove the last 3 lines if the file has at least 3 lines
        if len(lines) >= 3:
            lines = lines[:-3]

        with open(file_path, 'w') as file:
            file.write(''.join(lines))

        logger.info("Removed the last couple of lines from {}".format(file_path))
    except Exception as e:
        raise RuntimeError("Error modifying {}: {}".format(file_path, e))

def keep_only_last_two_lines(file_path):
    """
    Keeps only the last two lines in the given file, removing all other lines.
    """
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Keep only the last two lines (if the file has at least two)
        if len(lines) > 2:
            lines = lines[-2:]

        # Overwrite the file with the reduced content
        with open(file_path, 'w') as file:
            file.write(''.join(lines))
        
        logger.info("Kept only the last two lines in {}".format(file_path))
    except Exception as e:
        raise RuntimeError("Error modifying {}: {}".format(file_path, e))


def main():
    try:
        current_working_dir = os.getcwd()
        backup_files(current_working_dir)

        # Copy special maneuver-manager files
        backup_dir = os.path.join(current_working_dir, "Reset")
        backup_manoeuver_dir = os.path.join(backup_dir, "Manoeuver")
        result_file_dir = "/app/maneuver-manager"
        
        special_files = [
            "Result.txt",
            "LastManeuverDate.txt",
            "PostManeuverDate.txt",
            "Result2.txt",
            "LastManeuverDate2.txt",
            "PostManeuverDate2.txt",
            "Result3.txt",
            "LastManeuverDate3.txt",
            "PostManeuverDate3.txt",
            "Result4.txt",
            "LastManeuverDate4.txt",
            "PostManeuverDate4.txt",
            "time-persistence.json"
        ]
        for special_file in special_files:
            src_path = os.path.join(result_file_dir, special_file)
            if os.path.exists(src_path):
                dst_path = os.path.join(backup_manoeuver_dir, special_file)
                safe_copy(src_path, dst_path)
                if special_file in ["LastManeuverDate.txt", "LastManeuverDate2.txt", "LastManeuverDate3.txt","LastManeuverDate4.txt"]:
                    remove_last_couple_of_lines(dst_path)
                    keep_only_last_two_lines(dst_path)

    except Exception as e:
        logger.error("Backup process failed: {}".format(e))
        raise

if __name__ == "__main__":
    main()
