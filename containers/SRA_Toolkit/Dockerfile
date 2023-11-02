# Spécification de l'image Ubuntu 22.04 comme base
FROM ubuntu:22.04

# Mettez à jour le système et installez les dépendances nécessaires
RUN apt-get update && apt-get install -y wget

# Téléchargez le SRA Toolkit version 2.10.0 depuis le site officiel
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.0/sratoolkit.2.10.0-ubuntu64.tar.gz

# Décompressez l'archive
RUN tar xzf sratoolkit.2.10.0-ubuntu64.tar.gz

# Ajoutez le répertoire contenant les outils du SRA Toolkit à la variable PATH
ENV PATH="/sratoolkit.2.10.0-ubuntu64/bin:${PATH}"

# Supprimez les fichiers temporaires et inutiles
RUN rm -rf sratoolkit.2.10.0-ubuntu64.tar.gz

# Commande par défaut à exécuter lorsque le conteneur démarre
CMD ["fastq-dump"]
