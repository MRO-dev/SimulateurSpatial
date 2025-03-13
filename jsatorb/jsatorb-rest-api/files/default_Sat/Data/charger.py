# -*- coding: utf-8 -*-
# Anthony Le Batteux EAE 130723
import os
import shutil
import sys

def backup_files(source_dir):
    # Chemin du répertoire de sauvegarde
    backup_dir = os.path.join(source_dir, "Backup")
    
    # Supprime le dossier Backup s'il existe, puis le recrée
    if os.path.exists(backup_dir):
        shutil.rmtree(backup_dir)
    os.makedirs(backup_dir)
    
    # Liste des fichiers dans le répertoire source
    files = os.listdir(source_dir)

    # Copier chaque fichier .TXT dans le répertoire de sauvegarde
    for file in files:
        source_path = os.path.join(source_dir, file)
        backup_path = os.path.join(backup_dir, file)
        # Vérifie si c'est un fichier et s'il se termine par .TXT
        if os.path.isfile(source_path) and file.endswith(".TXT"):
            shutil.copy2(source_path, backup_path)
            print("Copied {} to {}".format(source_path, backup_path))

    # Copie tout le dossier "Manoeuver" dans "Backup"
    manoeuver_dir = os.path.join(source_dir, "Manoeuver")
    if os.path.exists(manoeuver_dir) and os.path.isdir(manoeuver_dir):
        backup_manoeuver_dir = os.path.join(backup_dir, "Manoeuver")
        shutil.copytree(manoeuver_dir, backup_manoeuver_dir)
        print("Copied entire 'Manoeuver' directory to {}".format(backup_manoeuver_dir))

def delete_txt_files(source_dir):
    # Supprime chaque fichier .TXT dans le répertoire courant
    for file in os.listdir(source_dir):
        file_path = os.path.join(source_dir, file)
        if os.path.isfile(file_path) and file.endswith(".TXT"):
            os.remove(file_path)
            print("Deleted {}".format(file_path))

def delete_manoeuver_dir(source_dir):
    # Supprime le dossier "Manoeuver" dans le répertoire courant
    manoeuver_dir = os.path.join(source_dir, "Manoeuver")
    if os.path.exists(manoeuver_dir) and os.path.isdir(manoeuver_dir):
        shutil.rmtree(manoeuver_dir)
        print("Deleted 'Manoeuver' directory from {}".format(source_dir))

def copy_backup_to_source(backup_dir, source_dir):
    """
    Copies files from the Backup folder to source_dir:
      - Regular .txt files go to the top-level of source_dir,
      - EXCEPT the special .txt files (listed below) which, along with all non-.txt files, go to source_dir/Manoeuver.
    """
    # List (set) of special .txt files that must NOT be copied to the top-level:
    special_txt_files = {
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
        "PostManeuverDate4.txt"
    }
    
    local_manoeuver = os.path.join(source_dir, "Manoeuver")
    if not os.path.exists(local_manoeuver):
        os.makedirs(local_manoeuver)

    # Walk through the entire Backup folder structure
    for root, dirs, files in os.walk(backup_dir):
        # rel_path is relative to the Backup folder
        rel_path = os.path.relpath(root, backup_dir)
        for file in files:
            backup_file_path = os.path.join(root, file)
            # Determine destination:
            # If the file is a .txt file AND NOT one of the special files, copy it to the top level.
            # Otherwise, copy it into the local Manoeuver folder.
            if file.lower().endswith(".txt") and file not in special_txt_files:
                dest_dir = source_dir
            else:
                dest_dir = os.path.join(local_manoeuver, rel_path)
                if not os.path.exists(dest_dir):
                    os.makedirs(dest_dir)
            dest_path = os.path.join(dest_dir, file)
            shutil.copy2(backup_file_path, dest_path)
            print("Copied {} to {}".format(backup_file_path, dest_path))

# ----------------------------------------------------------------
# Main "Loading" logic
# ----------------------------------------------------------------
print("TEST2")
currentWorkingDirectory = os.getcwd()

# 1. Supprime les fichiers .TXT du répertoire courant
delete_txt_files(currentWorkingDirectory)

# 2. Supprime le dossier "Manoeuver" du répertoire courant
delete_manoeuver_dir(currentWorkingDirectory)

# 3. Copie tout le contenu du dossier "Backup" vers le répertoire courant,
#    en mettant les .txt réguliers dans le top-level et les fichiers spéciaux dans Manoeuver
backup_directory = os.path.join(currentWorkingDirectory, "Backup")
copy_backup_to_source(backup_directory, currentWorkingDirectory)

# 4. Ensuite, copie certains fichiers depuis le dossier local "Manoeuver" vers /app/maneuver-manager
local_manoeuver = os.path.join(currentWorkingDirectory, "Manoeuver")
result_file_dir = "/app/maneuver-manager"

def copy_manoeuver_file_to_result(filename):
    source_file = os.path.join(local_manoeuver, filename)
    if os.path.exists(source_file):
        dest_file = os.path.join(result_file_dir, filename)
        shutil.copy2(source_file, dest_file)
        print("Copied '{}' from {} to {}".format(filename, source_file, dest_file))

files_to_copy_back = [
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
    "blue1-tle-mongo.json",
    "blue2-tle-mongo.json",
    "red1-tle-mongo.json",
    "red2-tle-mongo.json",
    "time-persistence.json"
]

for f in files_to_copy_back:
    copy_manoeuver_file_to_result(f)

# Fin du script.

