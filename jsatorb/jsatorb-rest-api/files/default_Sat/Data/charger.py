# Anthony Le Batteux EAE 130723
import os
import shutil
import sys

def backup_files(source_dir):
    # Creer le chemin du repertoire de sauvegarde
    backup_dir = os.path.join(source_dir, "Backup")
    
    # Verifie si le dossier de sauvegarde existe et le supprime s'il existe deja
    if os.path.exists(backup_dir):
        shutil.rmtree(backup_dir)  # Supprime entierement le dossier Backup
    
    # Recreer le dossier Backup apres l'avoir vide
    os.makedirs(backup_dir)
    
    # Liste des fichiers dans le repertoire source
    files = os.listdir(source_dir)

    # Copier chaque fichier .TXT dans le repertoire de sauvegarde
    for file in files:
        source_path = os.path.join(source_dir, file)
        backup_path = os.path.join(backup_dir, file)
        
        # Verifie si c'est un fichier et s'il se termine par .TXT
        if os.path.isfile(source_path) and file.endswith(".TXT"):
            # Copie le fichier (remplace s'il existe deja dans le dossier de backup)
            shutil.copy2(source_path, backup_path)
            print("Copied {} to {}".format(source_path, backup_path))

    # Chemin du dossier "Manoeuver"
    manoeuver_dir = os.path.join(source_dir, "Manoeuver")
    if os.path.exists(manoeuver_dir) and os.path.isdir(manoeuver_dir):
        # Chemin de destination du dossier "Manoeuver" dans le dossier "Backup"
        backup_manoeuver_dir = os.path.join(backup_dir, "Manoeuver")

        # Copie tout le contenu du dossier "Manoeuver" vers le dossier "Backup"
        shutil.copytree(manoeuver_dir, backup_manoeuver_dir)
        print("Copied entire 'Manoeuver' directory to {}".format(backup_manoeuver_dir))

def delete_txt_files(source_dir):
    # Liste des fichiers dans le repertoire source
    files = os.listdir(source_dir)

    # Supprime chaque fichier .TXT dans le repertoire courant
    for file in files:
        file_path = os.path.join(source_dir, file)
        if os.path.isfile(file_path) and file.endswith(".TXT"):
            os.remove(file_path)
            print("Deleted {}".format(file_path))

def delete_manoeuver_dir(source_dir):
    # Chemin du dossier "Manoeuver"
    manoeuver_dir = os.path.join(source_dir, "Manoeuver")
    # Verifie si le dossier existe, et le supprime s'il existe
    if os.path.exists(manoeuver_dir) and os.path.isdir(manoeuver_dir):
        shutil.rmtree(manoeuver_dir)  # Supprime entierement le dossier "Manoeuver"
        print("Deleted 'Manoeuver' directory from {}".format(source_dir))

def copy_backup_to_source(backup_dir, source_dir):
    # Liste des fichiers dans le dossier Backup
    for root, dirs, files in os.walk(backup_dir):
        # Obtenir le chemin relatif a partir du repertoire Backup
        rel_path = os.path.relpath(root, backup_dir)
        # Calculer le chemin de destination dans le repertoire source
        dest_dir = os.path.join(source_dir, rel_path)

        # Creez le dossier de destination s'il n'existe pas
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)

        # Copier chaque fichier dans le repertoire de destination
        for file in files:
            source_file = os.path.join(root, file)
            dest_file = os.path.join(dest_dir, file)
            shutil.copy2(source_file, dest_file)
            print("Copied {} to {}".format(source_file, dest_file))

print("TEST2")
currentWorkingDirectory = os.getcwd()
# Dossier "Manoeuver"
manoeuver_directory = os.path.join(currentWorkingDirectory, "Manoeuver")
# Supprimer les fichiers .TXT du repertoire courant
delete_txt_files(currentWorkingDirectory)
# Supprimer le dossier "Manoeuver" du repertoire courant
delete_manoeuver_dir(currentWorkingDirectory)
# Copier tout le contenu du dossier "Backup" vers le repertoire courant
backup_directory = os.path.join(currentWorkingDirectory, "Backup")
copy_backup_to_source(backup_directory, currentWorkingDirectory)

manoeuver_dir = os.path.join(currentWorkingDirectory , "Manoeuver")

backup_dir = os.path.join(currentWorkingDirectory, "Backup")
backup_manoeuver_dir = os.path.join(backup_dir, "Manoeuver")
# Chemin du repertoire contenant Result.txt
result_file_dir = "/app/maneuver-manager"

# Copie du fichier Result.txt a partir de son emplacement specifique
result_file = os.path.join(backup_manoeuver_dir, "Result.txt")
if os.path.exists(result_file):
    backup_result_path = os.path.join(result_file_dir, "Result.txt")
    shutil.copy2(result_file, backup_result_path)
    print("Copied 'Result.txt' from {} to {}".format(result_file, backup_result_path))
    
    # Copie du fichier Result.txt a partir de son emplacement specifique
result_file = os.path.join(backup_manoeuver_dir, "LastManeuverDate.txt")
if os.path.exists(result_file):
    backup_result_path = os.path.join(result_file_dir, "LastManeuverDate.txt")
    shutil.copy2(result_file, backup_result_path)
    print("Copied 'LastManeuverDate.txt' from {} to {}".format(result_file, backup_result_path))
    
# Copie du fichier Result.txt a partir de son emplacement specifique
result_file = os.path.join(backup_manoeuver_dir, "blue1-tle-mongo.json")
if os.path.exists(result_file):
    backup_result_path = os.path.join(result_file_dir, "blue1-tle-mongo.json")
    shutil.copy2(result_file, backup_result_path)
    print("Copied 'Result.txt' from {} to {}".format(result_file, backup_result_path))    

# Copie du fichier PostManeuverDate.txt a partir de son emplacement specifique
result_file = os.path.join(backup_manoeuver_dir, "PostManeuverDate.txt")
if os.path.exists(result_file):
    backup_result_path = os.path.join(result_file_dir, "PostManeuverDate.txt")
    shutil.copy2(result_file, backup_result_path)
    print("Copied 'PostManeuverDate.txt' from {} to {}".format(result_file, backup_result_path))   


