# Utilisez une image Ubuntu comme base
FROM ubuntu:22.04 

# Installation des dépendances et des logiciels (wget et unzip)
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Téléchargement et installation de TrimGalore
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.7.zip -O TrimGalore.zip && \
    unzip TrimGalore.zip && \
    mv TrimGalore-0.6.7/trim_galore /usr/local/bin/ && \
    rm -r TrimGalore-0.6.7 TrimGalore.zip

# Commande par défaut exécutée lorsque le conteneur démarre
CMD ["trim_galore", "--version"]
