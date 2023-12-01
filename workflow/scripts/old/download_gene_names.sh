#!/bin/bash
wget "https://www.genome.jp/kegg-bin/download_htext?htext=sao00001.keg&format=htext&filedir=" -O sao00001.keg

awk '/D      SAOUHSC_01493 rpsA; 30S ribosomal protein S1	K02945 RP-S1; small subunit ribosomal protein S1/,/D      SAOUHSC_00933 tryptophanyl-tRNA synthetase	K01867 WARS; tryptophanyl-tRNA synthetase [EC:6.1.1.2]/ {if (NR>1) print}' sao00001.keg> sao00002.keg

#Supression de tout ce qu'il y a avant SAOUHSC
sed 's/.*SAOUHSC/SAOUHSC/' sao00002.keg > sao00003.keg

#Suppression des derniers caractères pour n'avoir que le nom du gène
#awk '{print substr($0, 8, 13)}' sao00002.keg > sao00003.keg
awk '{
    if (NR >= 1 && NR <= 54) {
        print substr($0, 1, 13)
    } else if (NR >= 55 && NR <= 133) {
        print substr($0, 1, 14)
    } else {
        print substr($0, 1, 13)
    }
}' sao00003.keg > sao00004.keg

awk '/C    03013 Nu/{exit} {print}' sao00004.keg > sao00005.keg

grep "^SAOUHSC" sao00005.keg > sao00006.keg

# Noms des lignes à ajouter
lignes_a_ajouter=(
  "SAOUHSC_01236"
  "SAOUHSC_02489"
  "SAOUHSC_00475"
  "SAOUHSC_01234"
  "SAOUHSC_01786"
  "SAOUHSC_01246"
)

# Chemin du fichier où ajouter les lignes
fichier_destination="sao00006.keg"  # Remplacez par le nom de votre fichier

# Ajouter chaque ligne au fichier
for ligne in "${lignes_a_ajouter[@]}"; do
  echo "$ligne" >> "$fichier_destination"
done

#suppression des espaces en trop
sed 's/ //g' sao00006.keg > Gene_Names_1col_true.csv