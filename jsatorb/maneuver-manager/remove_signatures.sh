#!/bin/bash

# Vérifier si un fichier JAR est fourni en argument
if [ -z "$1" ]; then
  echo "Usage: $0 fichier.jar"
  exit 1
fi

# Vérifier si le fichier JAR existe
if [ ! -f "$1" ]; then
  echo "Erreur: le fichier $1 n'existe pas."
  exit 1
fi

# Nom du fichier JAR
JAR_FILE=$1
NEW_JAR_FILE="${JAR_FILE%.jar}-unsigned.jar"

# Créer un répertoire temporaire pour extraire le contenu du JAR
TMP_DIR=$(mktemp -d)

# Extraire le contenu du fichier JAR dans le répertoire temporaire
echo "Extraction du fichier JAR..."
unzip -q "$JAR_FILE" -d "$TMP_DIR"
if [ $? -ne 0 ]; then
  echo "Erreur lors de l'extraction du fichier JAR."
  rm -rf "$TMP_DIR"
  exit 1
fi

# Vérifier si le répertoire META-INF existe
if [ ! -d "$TMP_DIR/META-INF" ]; then
  echo "Le répertoire META-INF n'existe pas dans ce JAR."
  rm -rf "$TMP_DIR"
  exit 1
fi

# Supprimer les fichiers de signature dans META-INF
echo "Suppression des fichiers de signature dans META-INF..."
find "$TMP_DIR/META-INF" -type f \( -name "*.SF" -o -name "*.RSA" -o -name "*.DSA" \) -exec rm -f {} +

# Recréer le fichier JAR sans les fichiers de signatures
echo "Création d'un nouveau fichier JAR sans signatures..."
cd "$TMP_DIR"
zip -r -q "$NEW_JAR_FILE" . 
if [ $? -ne 0 ]; then
  echo "Erreur lors de la création du nouveau fichier JAR."
  cd -
  rm -rf "$TMP_DIR"
  exit 1
fi

# Déplacer le nouveau fichier JAR à l'emplacement d'origine
mv "$NEW_JAR_FILE" "$OLDPWD"
cd -

# Nettoyage du répertoire temporaire
rm -rf "$TMP_DIR"

echo "Le fichier signé a été recréé sous le nom $NEW_JAR_FILE sans les fichiers de signature."
