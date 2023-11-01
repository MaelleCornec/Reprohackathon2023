#Spécification de l'image de base qu'on souhaite utiliser 
FROM ubuntu:22.04 

#Installation des dépendances et des logiciels
RUN apt-get update && apt-get install -y build-essential

####### Téléchargement de Bowtie Version 0.12.7 depuis la source officielle pour l'alignement des séquences
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip/download
RUN unzip download
RUN touch bowtie-build.sh
RUN echo bowtie-build 

#RUN ./bowtie-0.12.7/bowtie-build ../reference.fasta ../reference.gff

# Par exemple, pour vérifier que Bowtie est installé correctement
CMD ["bowtie", "--version"]