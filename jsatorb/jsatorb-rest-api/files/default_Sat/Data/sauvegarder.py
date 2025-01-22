# Anthony Le Batteux EAE 130723
import os
import shutil
import sys

def backup_files(source_dir):
    print("TEST")
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

print("TEST2")
currentWorkingDirectory = os.getcwd()
backup_files(currentWorkingDirectory)
backup_dir = os.path.join(currentWorkingDirectory, "Backup")
backup_manoeuver_dir = os.path.join(backup_dir, "Manoeuver")

# Chemin du repertoire contenant Result.txt
result_file_dir = "/app/maneuver-manager"

manoeuver_dir = os.path.join(currentWorkingDirectory , "Manoeuver")

# Copie du fichier Result.txt a partir de son emplacement specifique
result_file = os.path.join(result_file_dir, "Result.txt")
if os.path.exists(result_file):
    shutil.copy2(result_file, backup_manoeuver_dir)
    print("Copied 'Result.txt' from {} to {}".format(result_file, backup_manoeuver_dir))
    
    # Copie du fichier LastManeuverDate.txt a partir de son emplacement specifique
result_file = os.path.join(result_file_dir, "LastManeuverDate.txt")
if os.path.exists(result_file):
    shutil.copy2(result_file, backup_manoeuver_dir)
    print("Copied 'LastManeuverDate.txt' from {} to {}".format(result_file, backup_manoeuver_dir))


# Copie du fichier PostManeuverDate.txt a partir de son emplacement specifique
result_file = os.path.join(result_file_dir, "PostManeuverDate.txt")
if os.path.exists(result_file):
    shutil.copy2(result_file, backup_manoeuver_dir)
    print("Copied 'PostManeuverDate.txt' from {} to {}".format(result_file, backup_manoeuver_dir))   


