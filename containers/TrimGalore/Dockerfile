# Spécification de l'image Ubuntu comme base
FROM ubuntu:22.04 

# Mise à jour du système et installation des dépendances nécessaires
RUN apt-get update && apt-get install -y wget

# Téléchargement de Trim Galore version 0.6.4 depuis le site officiel
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.4.tar.gz

# Décompression de l'archive
RUN tar xzf 0.6.4.tar.gz

# Ajout du répertoire contenant Trim Galore à la variable PATH
ENV PATH="/TrimGalore-0.6.4:${PATH}"

# Suppresion des fichiers temporaires et inutiles
RUN rm -rf 0.6.4.tar.gz

# Commande par défaut à exécuter lorsque le conteneur démarre
CMD ["trim_galore"]