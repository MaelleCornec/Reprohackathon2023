# Téléchargement du fichier kegg de tous les gènes du staphylocoque doré
wget "https://www.genome.jp/kegg-bin/download_htext?htext=sao00001.keg&format=htext&filedir=" -O sao00001.keg

# Suppression des doublons
awk '!seen[$0]++' sao00001.keg > sao00002.keg

# Sélection dans le document téléchargé des lignes qui nous intéressent
awk '
/D      SAOUHSC_01493 rpsA; 30S ribosomal protein S1	K02945 RP-S1; small subunit ribosomal protein S1/,
/D      SAOUHSC_R00011 5S Ribosomal RNA	K01985 5SrRNA; 5S ribosomal RNA/ {if (NR > 1) print}
/D      SAOUHSC_T0001 tRNA-Ala	K14218 tRNA-Ala; tRNA Ala/,
/D      SAOUHSC_T00060 tRNA-Val	K14237 tRNA-Val; tRNA Val/ {if (NR > 1) print}
/D      SAOUHSC_00509 gltX; glutamyl-tRNA synthetase	K09698 gltX; nondiscriminating glutamyl-tRNA synthetase/,
/D      SAOUHSC_00933 tryptophanyl-tRNA synthetase	K01867 WARS; tryptophanyl-tRNA synthetase/ {if (NR > 1) print}
/D      SAOUHSC_01246 infB; translation initiation factor IF-2	K02519 infB; translation initiation factor IF-2/,
/D      SAOUHSC_00956 prfC; peptide chain release factor 3	K02837 prfC; peptide chain release factor 3/ {if (NR > 1) print}
/D      SAOUHSC_00892 hypothetical protein	K07570 GSP13; general stress protein 13/,
/D      SAOUHSC_00483 hypothetical protein	K07571 K07571; S1 RNA binding domain protein/ {if (NR > 1) print}' sao00002.keg > sao00003.keg

# Supression de tout ce qu'il y a avant SAOUHSC
sed 's/.*SAOUHSC/SAOUHSC/' sao00003.keg > sao00004.keg

# Suppression de tout ce qui se trouve après le premier espace pour n'avoir que le nom du gène
awk '{sub(/ .*/, ""); print}' sao00004.keg > scripts/genes/translation_genes.csv

# Suppression des fichiers intermédiaires
rm sao00001.keg sao00002.keg sao00003.keg sao00004.keg