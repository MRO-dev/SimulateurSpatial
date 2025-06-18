**Guide d’utilisation du simulateur spatial (ESS)**

L'ESS (Educational Space Simulator) est un prototype conçu pour réaliser des manœuvres orbitales sur des satellites définis par l'utilisateur.

Ce simulateur est actuellement déployé au sein d'une topologie Hynesim sur le serveur **172.16.38.128** et dans la topologie **MCO** (ce nom pouvant être modifié ultérieurement).

Ce guide d'utilisation vise à faciliter une prise en main rapide du prototype du simulateur spatial.

Les sources du projet se trouvent dans le répertoire git suivant :\
[MRO-dev/SimulateurSpatial: Projet de simulateur spatial](https://github.com/MRO-dev/SimulateurSpatial)

Le projet initial est issu du travail de M. Anthony LEBATTEUX :\
[AnthonyLB2140/SIMU](https://github.com/AnthonyLB2140/SIMU/tree/master)






















1. **Tableau de bord**

![Une image contenant texte, capture d’écran, Logiciel multimédia, logiciel

Le contenu généré par l’IA peut être incorrect.](Aspose.Words.3bae6ac2-a381-4f90-870e-74b47082fe6b.001.png)

**I.1. 	Présentation générale du tableau de bord** 

Le tableau de bord comprend une page principale dotée de plusieurs boutons, dont les fonctionnalités sont détaillées ci-dessous.

Cette page principale est divisée en plusieurs onglets accessibles par navigation. Les onglets utilisés dans le cadre de ce projet sont les suivants :

- **« MCS Home » :** permet de revenir au tableau de bord principal.
- **« MCS Maneuver Manager – BLUE1 » :** conduit à la page spécifique au satellite BLUE1, permettant de gérer les manœuvres orbitales de ce satellite.
- **« MCS Maneuver Manager – BLUE2 » :** conduit à la page spécifique au satellite BLUE2, permettant de gérer les manœuvres orbitales de ce satellite.
- **« MCS Maneuver Manager – RED1 » :** conduit à la page spécifique au satellite RED1, permettant de gérer les manœuvres orbitales de ce satellite.
- **« MCS Maneuver Manager – RED2 » :** conduit à la page spécifique au satellite RED2, permettant de gérer les manœuvres orbitales de ce satellite.

Les autres onglets ne sont pas utilisés dans le cadre de ce projet.


**II.2.	Fonctionnalités du tableau de bord**


**Définition de l’horizon de simulation**

Cette fonctionnalité permet de fixer une période de simulation en indiquant une date de début et une date de fin. Toutes les manœuvres orbitales effectuées doivent obligatoirement s’inscrire dans cet intervalle temporel.


![Une image contenant texte, capture d’écran, logiciel, Logiciel multimédia

Le contenu généré par l’IA peut être incorrect.](Aspose.Words.3bae6ac2-a381-4f90-870e-74b47082fe6b.002.png)











**Description des boutons du DashBoard**

![Une image contenant texte, capture d’écran, logiciel, Système d’exploitation

Le contenu généré par l’IA peut être incorrect.](Aspose.Words.3bae6ac2-a381-4f90-870e-74b47082fe6b.003.png)
























- **« DELETE MEMS » :**\
  Supprime les fichiers MEMs contenant les données des trajectoires des satellites utilisateurs (BLUE/RED).
- **« DELETE TLE STORAGE » :**\
  Supprime l’ensemble des données TLE ainsi que les masses des satellites utilisateurs. Utile particulièrement pour démarrer une toute nouvelle simulation.
- **« SAVE SITUATION » :**\
  Permet de sauvegarder la configuration en cours.
- **« TLE AND VISI WINDOWS UPDATE » :**\
  Réinitialise les TLE des satellites selon la configuration par défaut du logiciel. Peu utile dans ce cadre, car les satellites seront reconfigurés par l’utilisateur.
- **« Simulation Time Reset » :**\
  Réinitialise l’horizon de simulation aux dates actuelles.

**Remarque :** En cas de problème de sauvegarde, notamment lors de la création d'un nouveau dossier, il est parfois nécessaire de soumettre à nouveau la sauvegarde. Si les difficultés persistent, rafraîchissez la page et recommencez l'opération. À la fin, les données doivent être correctement enregistrées dans le dossier choisi.

- **« LOAD SITUATION » :**\
  Permet de charger une configuration existante. Veillez à spécifier un dossier valide contenant des données réelles.
- **« RESCUE TLE DATABASE » :**\
  Actualise la base de données des satellites prédéfinis (ex. PLEIADES). Intéressant lors d'une nouvelle simulation.
- **« REFRESH VTS » :**\
  Génère à nouveau les fichiers MEMs des satellites utilisateurs pour la période définie.
- **« RESCUE TLE AIGLONS » :**\
  Met à jour la base des TLE des satellites Aiglons. **Inutile pour cette étude.**
- **« REFRESH VTS AIGLONS » :**\
  Recrée les fichiers MEMs des satellites Aiglons. **Inutile pour cette étude.**
- **« DELETE ALL MEMS » :**\
  Efface l'ensemble des fichiers MEMs existants. Fonction indispensable pour démarrer un nouvel horizon de simulation.

**Important :** Après avoir exécuté cette opération, il est nécessaire d'effectuer dans l’ordre :

1. un **RESCUE TLE DATABASE**,
1. puis un **REFRESH ALL VTS SATS**.
- **« REFRESH ALL VTS SATS » :**\
  Régénère les fichiers MEMs de tous les satellites. Essentiel dans le cadre d'une simulation entièrement nouvelle.












**II.	Réalisation d’une manœuvre sur un satellite**

![Une image contenant texte, capture d’écran, logiciel, Logiciel multimédia

Le contenu généré par l’IA peut être incorrect.](Aspose.Words.3bae6ac2-a381-4f90-870e-74b47082fe6b.004.jpeg)

**II.1.	Initialisation des paramètres**

Initialisation des paramètres selon le formulaire suivant :

![Une image contenant texte, logiciel, Logiciel multimédia, Icône d’ordinateur

Le contenu généré par l’IA peut être incorrect.](Aspose.Words.3bae6ac2-a381-4f90-870e-74b47082fe6b.005.jpeg)

- Vérifiez et validez les paramètres suivants : inclinaison, RAAN, excentricité, AoP, anomalie moyenne, SMA, masse sèche (Dry Mass) et masse initiale d’ergol (Initial Ergol Mass).
- Sélectionnez ensuite une date initiale, puis cliquez sur **« B1 TLE UPDATE »** et validez afin d'initialiser le satellite.

**II.2.	Réalisation d’une manœuvre manuelle**


![Une image contenant texte, capture d’écran, logiciel, Logiciel multimédia

Le contenu généré par l’IA peut être incorrect.](Aspose.Words.3bae6ac2-a381-4f90-870e-74b47082fe6b.006.png)



















**Procédure pour paramétrer une manœuvre :**

- Définir les angles d’inclinaison des satellites Alpha et Beta.
- Sélectionner une date initiale pour la manœuvre
- Choisir le type de manœuvre (**Maneuver Type**) : « Impulse » (seul choix disponible).
- Définir les paramètres de poussée :
  - ΔV (en m/s)
  - ISP (en secondes)
- Cliquer sur le bouton **« SUBMIT »** pour valider la saisie.


**II.3. Réalisation d’une manœuvre prédéfinie**

![Une image contenant texte, Appareils électroniques, capture d’écran, logiciel

Le contenu généré par l’IA peut être incorrect.](Aspose.Words.3bae6ac2-a381-4f90-870e-74b47082fe6b.007.png)

































Pour effectuer une manœuvre prédéfinie, suivez les étapes suivantes :

- Sélectionnez une manœuvre parmi les choix disponibles :\
  (**Hohmann**, **Incli**, **Phasage**, **Bi-elliptic**).
- Complétez les paramètres spécifiques requis selon la manœuvre sélectionnée.
- Choisissez la date prévue pour réaliser la manœuvre.
- Définissez la valeur de l’ISP (en secondes).
- Cliquez sur le bouton **« SUBMIT »** pour valider la manœuvre.


**II.3. 	Reset**

Le bouton **RESET**, situé dans l’interface de configuration du satellite, permet de restaurer l’état du satellite tel qu’il était juste avant l’exécution de la dernière manœuvre.\
Cette fonctionnalité est utile pour tester différentes manœuvres en ayant la possibilité de revenir rapidement à l'état antérieur.

