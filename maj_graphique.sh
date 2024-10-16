#!/bin/bash

# Fonction pour mettre à jour les pilotes graphiques
update_graphics_drivers() {
    echo "Mise à jour des pilotes graphiques et de MESA..."
    sudo apt update && sudo apt upgrade -y
    
    # Installer ou mettre à jour les paquets MESA
    sudo apt install --only-upgrade mesa-utils libgl1-mesa-dri libglx-mesa0 libegl-mesa0 libgles2-mesa libosmesa6 -y

    # Vérifier si une carte NVIDIA est présente, et mettre à jour les pilotes propriétaires
    if lspci | grep -i nvidia; then
        echo "Carte NVIDIA détectée. Mise à jour des pilotes NVIDIA..."
        sudo apt install --only-upgrade nvidia-driver-470 -y
    fi

    # Vérifier si une carte AMD est présente, et mettre à jour les pilotes
    if lspci | grep -i amd; then
        echo "Carte AMD détectée. Mise à jour des pilotes AMD..."
        sudo apt install --only-upgrade mesa-vulkan-drivers mesa-vdpau-drivers -y
    fi

    echo "Mise à jour terminée. Un redémarrage peut être nécessaire pour appliquer les changements."
}

# Fonction pour désactiver l'accélération matérielle
disable_hardware_acceleration() {
    echo "Désactivation de l'accélération matérielle dans Qt..."
    
    # Créer un fichier de configuration Qt si nécessaire
    QT_CONFIG_DIR="$HOME/.config/QtProject"
    mkdir -p "$QT_CONFIG_DIR"

    # Ajouter la désactivation de l'accélération matérielle
    echo "[Qt]" > "$QT_CONFIG_DIR/qt.conf"
    echo "Platform=xcb" >> "$QT_CONFIG_DIR/qt.conf"
    echo "Désactivation de l'accélération matérielle ajoutée dans $QT_CONFIG_DIR/qt.conf."
}

# Fonction principale
main() {
    echo "Script de gestion des pilotes graphiques et de désactivation de l'accélération matérielle"
    
    # Mettre à jour les pilotes graphiques
    read -p "Voulez-vous mettre à jour les pilotes graphiques et MESA ? (y/n) " update_choice
    if [[ $update_choice == "y" || $update_choice == "Y" ]]; then
        update_graphics_drivers
    fi

    # Désactiver l'accélération matérielle
    read -p "Voulez-vous désactiver l'accélération matérielle pour Qt ? (y/n) " disable_choice
    if [[ $disable_choice == "y" || $disable_choice == "Y" ]]; then
        disable_hardware_acceleration
    fi

    echo "Script terminé."
}

# Lancer le script principal
main
