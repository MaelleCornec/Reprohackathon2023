#Spécification de l'image de base qu'on souhaite utiliser 
FROM ubuntu:22.04 

# Mettez à jour le système et installez les dépendances nécessaires
RUN apt-get update && apt-get install -y wget

# Téléchargez featureCounts version 1.4.6-p3 depuis le site officiel
RUN wget https://downloads.sourceforge.net/project/subread/subread-1.4.6-p3/subread-1.4.6-p3-Linux-x86_64.tar.gz

# Décompressez l'archive
RUN tar xzf subread-1.4.6-p3-Linux-x86_64.tar.gz

# Ajoutez le répertoire contenant featureCounts à la variable PATH
ENV PATH="/subread-1.4.6-p3-Linux-x86_64/bin:${PATH}"

# Supprimez les fichiers temporaires et inutiles
RUN rm -rf subread-1.4.6-p3-Linux-x86_64.tar.gz

# Commande par défaut à exécuter lorsque le conteneur démarre
CMD ["featureCounts"]

