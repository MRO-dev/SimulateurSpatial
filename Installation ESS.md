**Guide d'installation du simulateur spatial (ESS)**

L’ESS (Educational Space Simulator) est un prototype permettant de réaliser des manœuvres orbitales sur des satellites définis par l’utilisateur. Actuellement, ce simulateur est déployé dans une topologie sous Hynesim (serveur 172.16.38.128 et topologie MCO, dont le nom peut évoluer).

Ce guide présente les étapes nécessaires pour installer le simulateur ESS sur une nouvelle machine, en cas de besoin d'un déploiement à partir de zéro. Le simulateur peut être installé sous Linux, mais également sous Windows moyennant certaines adaptations.

Les sources du projet se trouvent dans le répertoire git suivant :\
[MRO-dev/SimulateurSpatial: Projet de simulateur spatial](https://github.com/MRO-dev/SimulateurSpatial)

Le projet initial est issu du travail de M. Anthony LEBATTEUX :\
[AnthonyLB2140/SIMU](https://github.com/AnthonyLB2140/SIMU/tree/master)

Les instructions ci-dessous concernent spécifiquement l’installation sous Linux :

**1. Installation préalable des outils :**

- **Docker** avec docker compose (et non docker-compose).
- **Outil de visualisation des trajectoires VTS** : disponible via Timeloop (une adresse e-mail est requise pour obtenir le lien de téléchargement).
- Une version compatible de **Node.js** et **npm** (indispensables pour Node-Red).
- **Node-Red** : disponible sur Node-RED.
- Les librairies précisées dans le fichier **DockerFile** (cf. NodeRed-Dockerfile).

**2. Installation des images Docker :**

Exécutez les commandes suivantes pour récupérer les images Docker nécessaires :

docker pull anthonyeae/reposetoibien:monnode

docker pull anthonyeae/reposetoibien:monjsat

docker pull anthonyeae/reposetoibien:moncelestrak

docker pull anthonyeae/reposetoibien:monmongo

**3. Récupération du code et des utilitaires :**

- Copier les dossiers **jsatorb** et **mosquitto** vers votre répertoire $HOME.
- Copier les scripts présents dans .local/share/applications vers le dossier d'applications utilisateur (typiquement .local/share/applications pour une distribution Debian classique).
- Copier les icônes associées dans .local/share/icons.
- Copier le fichier **flows.json** du répertoire .nodered du git vers votre répertoire local .node-red.

**4. Lancement du simulateur :**

- Lancez les conteneurs Docker suivants :

  - Docker **jsatorb** (icône EAE) 

    ![Une image contenant texte, Police, logo, graphisme

Le contenu généré par l’IA peut être incorrect.](Aspose.Words.e359b00d-41a8-4161-a519-9a4f8e625838.001.png)


  - Docker **Node-Red** (icône par défaut de Node-Red). 

![Une image contenant symbole, Panneau de stop, Panneau de signalisation, rouge

Le contenu généré par l’IA peut être incorrect.](Aspose.Words.e359b00d-41a8-4161-a519-9a4f8e625838.002.png)



- L’icône **VTS** sert à lancer l’interface de visualisation des satellites avec notamment une modélisation 3D ainsi qu’un planisphère. 

![](Aspose.Words.e359b00d-41a8-4161-a519-9a4f8e625838.003.png)



![](Aspose.Words.e359b00d-41a8-4161-a519-9a4f8e625838.004.png)

Une fois ces services démarrés, une page web contenant le tableau de bord du simulateur spatial s’ouvrira automatiquement.

