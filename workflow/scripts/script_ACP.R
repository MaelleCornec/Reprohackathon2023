#############################################
# Téléchargement des librairies nécessaires #
#############################################

library(ggplot2)
library(dplyr)
library(factoextra)

###################################################################
# Chargement des données nécessaires et premières transformations #
###################################################################

# Données de comptage #############################################

# Chargement fichier counts.txt dans la variable full_data
full_data <- read.table("scripts/genes/counts.txt", header = TRUE, skip = 1)
# Suppression des 5 premiers caractères de la colonne Geneid
full_data$Geneid <- substring(full_data$Geneid, first = 6)
# Renommage des colonnes 7 à 12
colnames(full_data)[7:12] <- c("Treatment_1",
                               "Treatment_2",
                               "Treatment_3",
                               "Control_1",
                               "Control_2",
                               "Control_3")

# Gènes obtenus par la documentation ##############################

# Chargement fichier Gene_Names dans la variable genes_names
genes_names <- read.csv("scripts/genes/translation_genes.csv",
                        header = TRUE)
# Sélection des gènes dans full_data correspondant aux gènes de genes_names
indices <- which(full_data$Geneid %in% genes_names$Name)
translation_genes <- full_data[indices, c(1:12)]

# Création des dataframes utiles à l'analyse ######################

# Création du dataframe de comptage des échantillons dans la variable count_data
count_data <- full_data[c(7:12)]
count_data <- na.omit(count_data)

##################################################
# Analyse des données totales à l'aide de DESeq2 #
##################################################

# Construction des figures d'ACP #################

# Calcul de l'ACP avec prcomp()
# Utilisation de scale. = TRUE pour centrer et réduire les variables
res_acp <- prcomp(count_data, scale. = TRUE)
# Résumé des résultats de l'ACP et affichage des composantes principales
summary(res_acp)
print(res_acp)
# Explication fournie par la variance de chaque composante principale
table_acp <- factoextra::get_eig(res_acp)
# Visualisation des cercles des variables dans le plan factoriel de l'ACP
png("reports/PCA_all_genes.png", width = 800, height = 800)
fviz_pca_var(res_acp,
             col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
dev.off()

#########################################################
# Analyse des données de translation à l'aide de DESeq2 #
#########################################################

# Construction des figures d'ACP ########################

# Calcul de l'ACP avec prcomp()
# Utilisation de scale. = TRUE pour centrer et réduire les variables
res_acp_trans <- prcomp(translation_genes[,c(7,8,9,10,11,12)], scale. = TRUE)
# Résumé des résultats de l'ACP et affichage des composantes principales
summary(res_acp_trans)
print(res_acp_trans)
# Explication fournie par la variance de chaque composante principale
table_acp_trans <- factoextra::get_eig(res_acp_trans)
# Visualisation des cercles des variables dans le plan factoriel de l'ACP
png("reports/PCA_translation_genes.png", width = 800, height = 800)
fviz_pca_var(res_acp_trans,
             col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
dev.off()
