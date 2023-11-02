# Déplacement du fichier fasta du génome de référence dans le dossier avec le Dockerfile et l'exécutable
# cp /home/ubuntu/Reprohackathon2023/reference.fasta /home/ubuntu/Reprohackathon2023/containers/Bowtie/

# Construction de l'image Docker
docker build -t bowtie:0.12.7 /home/ubuntu/Reprohackathon2023/containers/Bowtie

# Exécution de l'image
docker run -it bowtie:0.12.7

# Indexation du génome de référence 
#bowtie-build reference.fasta genome_ref_index
