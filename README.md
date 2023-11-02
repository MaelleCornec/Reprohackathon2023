# Reprohackathon2023
Repôt pour le projet dans le cadre de l'UC Reprohackathon - Formation IODAA APT/U. P-S

Idée pour créer le workflow :   
    * Créer les conteneurs suivants (nom de l'image : description):   
        * cutadapt version 1.11 ??? Seuls les séquences de plus de 25 nucléotides ont été considérés pour la suite.   
        Est-ce que c'est la même chose que TrimGalore ?    
        * trimGalore:v1 : un conteneur pour créer un environnement où l'outil TrimDalore sera utilisé.    
        * bowtie:v1 : un conteneur pour créer un environnement où la version 0.12.7 avec ses paramètres par défaut de l'outil Bowtie sera utilisé.   
        * featureCount:v1 : un conteneur pour créer un environnement où la version 1.4.6-p3 avec ses paramètres par défaut de l'outil Bowtie sera utilisé (from Subreads package (parameters: -t gene -g ID -s 1) Je ne sais paq ce que ça veut dire ???).    
        * deSeq2:v1 : un conteneur pour créer un environnement où la version 1.16 avec ses paramètres par défaut de l'outil DESeq2 sera utilisé.    
    * Créer les scripts (python ou bash ou etc.) pour chaque étape du workflow :    
        * trimming.py : le script pour éliminer les séquences de moins de 25 nucléotides les nucléotides pas assez précis (ce n'est pas obligé d'être en pyhthon mais c'est pour l'idée).    
        * indexing.py : le script pour créer les index vis-à-vis du génôme de référence.   
        * mapping.py : le script pour poser des index sur les séquences téléchargées.   
        * counting.py : le script pour compter les motifs (je ne suis pas sure de comprendre ce qu'il se passe ici).   
        * analysis.py (ou peut-être analysis.r avec R version 3.4.1 ?) : le script pour l'analyse des données (enfin !)    
    * Créer le worflow prenant en compte tous ces conteneurs et les scripts d'analyse des fichiers fastq initiaux.   
   
Information pour plus tard (extrait de l'article) :   
"For over-representation analysis, S. aureus KEGG gene-sets were
downloaded thanks to the EnrichmentBrowser R package version 2.14.3 (organism
code sao). All the 106 KEGG sets were then tested for the over-representation in
differentially expressed genes using the Fisher statistical test. Only gene-sets with a
FDR lower than 0.05 were considered significantly enriched."