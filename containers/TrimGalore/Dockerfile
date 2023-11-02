# Utilisez une image Ubuntu comme base
FROM ubuntu:22.04 

# Mettez à jour le système et installez les dépendances nécessaires
RUN apt-get update && apt-get install -y wget

# Télécharge Trim Galore! version 0.6.4 depuis le site officiel
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.4.tar.gz

# Décompressez l'archive
RUN tar xzf 0.6.4.tar.gz

# Ajoutez le répertoire contenant Trim Galore! à la variable PATH
ENV PATH="/TrimGalore-0.6.4:${PATH}"

# Supprimez les fichiers temporaires et inutiles
RUN rm -rf 0.6.4.tar.gz

# Commande par défaut à exécuter lorsque le conteneur démarre
CMD ["trim_galore"]


# Installation des dépendances et des logiciels (wget et unzip)
#RUN apt-get update && apt-get install -y \
#    wget \
#    unzip \
#    && rm -rf /var/lib/apt/lists/*

# Téléchargement et installation de TrimGalore version 0.6.7
#RUN wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.7.zip -O TrimGalore.zip && \
#    unzip TrimGalore.zip && \
#   mv TrimGalore-0.6.7/trim_galore /usr/local/bin/ && \
#   rm -r TrimGalore-0.6.7 TrimGalore.zip

# Commande par défaut exécutée lorsque le conteneur démarre
#CMD ["trim_galore", "--version"]
