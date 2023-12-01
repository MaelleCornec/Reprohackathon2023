#############################################
# Téléchargement des librairies nécessaires #
#############################################

library(ggplot2)
library(dplyr)
library(BiocManager)
library(DESeq2)
#library(factoextra)

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
genes_names <- read.csv("scripts/genes/Gene_Names_1col_true.csv",
                        header = TRUE,
                        sep = ";")
# Sélection des gènes dans full_data correspondant aux gènes de genes_names
indices <- which(full_data$Geneid %in% genes_names$Name)
translation_genes <- full_data[indices, c(1:12)]

# Création des dataframes utiles à l'analyse ######################

# Création du dataframe de comptage des échantillons dans la variable count_data
count_data <- full_data[c(7:12)]
count_data <- na.omit(count_data)
# Création du dataframe des conditions dans la variable col_data
col_data <- data.frame(condition = c("Treatment",
                                     "Treatment",
                                     "Treatment",
                                     "Control",
                                     "Control",
                                     "Control"))

##################################################
# Analyse des données totales à l'aide de DESeq2 #
##################################################

# Construction du MA-plot ########################

# Construction d'un objet DESeqDataSet
dds_full <- DESeqDataSetFromMatrix(countData = count_data,
                                   colData = col_data,
                                   design = ~ condition)
# Exécution de DESeq2
dds_full <- DESeq(dds_full)
# Nombre de gènes significatifs up-régulés et down-régulés
res <- results(dds_full, alpha = 0.05)
head(res)
summary(res)
# Enregistrement au format png du MA-plot fait à partir de res
png("reports/MA-plot_all_genes.png", width = 800, height = 600)
plotMA(res,
       main = "MA-plot, all genes",
       colSig = "red",
       cex = 0.9)
# Légende pour les points significatifs et non-significatifs
legend("bottomleft",
       legend = c("Significatif", "Non-significatif"),
       pch = 16,
       col = c("red", "black"))
dev.off()

# Construction des figures d'ACP #################

# Calcul de l'ACP avec prcomp()
# Utilisation de scale. = TRUE pour centrer et réduire les variables
res_acp <- prcomp(count_data, scale. = TRUE)
# Résumé des résultats de l'ACP et affichage des composantes principales
summary(res_acp)
print(res_acp)
# Explication fournie par la variance de chaque composante principale
#table_acp <- factoextra::get_eig(res_acp)
# Visualisation des cercles des variables dans le plan factoriel de l'ACP
#png("reports/PCA_all_genes.png", width = 800, height = 600)
#fviz_pca_var(res_acp,
#             col.var = "cos2",
#             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#             repel = TRUE)
#dev.off()

# Construction du Volcano plot ###################

# Création de l'image du Volcano plot
png("reports/Volcano_plot_all_genes.png", width = 800, height = 600)
# Volcano plot: Log2FoldChange en fonction de - ln de la p-value ajustée
plot(x = res$log2FoldChange,
     y = -log10(res$padj),
     main = "Volcano-plot, all genes",
     ylim = c(0, 20),
     cex = 0.5,
     xlab = bquote(~log[2] ~ fold ~ change),
     ylab = bquote(~-log[10] ~ Q ~ value))
abline(h = 0, col = "black")
dev.off()

#########################################################
# Analyse des données de translation à l'aide de DESeq2 #
#########################################################

# Sélection des lignes des gènes de translation dans la variable res précédente
resultat_final <- res[indices,
                      c("log2FoldChange", "baseMean", "pvalue", "padj")]
#Définition des 6 gènes particuliers
genes_to_label <- c("SAOUHSC_01236", "SAOUHSC_02489", "SAOUHSC_00475",
                    "SAOUHSC_01234", "SAOUHSC_01786", "SAOUHSC_01246")
genes_code <- c("frr", "infA", "pth", "tsf", "infC", "infB")
# Création de l'index pour les 6 gènes particuliers
index <- c()
for (i in seq_along(genes_to_label)) {
  gene <- genes_to_label[i]
  index <- append(index,
                  which(translation_genes$Geneid == gene))
}

# Création de l'index pour les points à encadrer en bleu
# Liés à l'aminoacyl ARNt synthétase
# Trouver l'indice où Geneid == "SAOUHSC_T0001"
start_index <- which(genes_names$Name == "SAOUHSC_T0001")
# Liste des gènes liés à la synthèse de l'aminoacyl ARNt
noms_genes_arnt <- genes_names$Name[start_index:162]
# Liste des index dans GENES_TRANSLATION_DATA
index_to_encircle <- c()
for (i in seq_along(noms_genes_arnt)) {
  index_to_encircle <- append(index_to_encircle,
                       which(translation_genes$Geneid == noms_genes_arnt[i]))
}

# Construction du MA-plot ###############################

#Plot classique
png("reports/MA-plot_translation_genes.png", width = 800, height = 600)
plotMA(resultat_final,
       main = "MA-plot, translation genes",
       colSig = "red",
       cex = 0.9)
#text(x = resultat_final$baseMean[index],
#     y = resultat_final$log2FoldChange[index],
#     pos = 4,
#     labels = genes_code,
#     cex = 1)
# Légende pour les points significatifs et non-significatifs
legend("bottomleft",
       legend = c("Significatif", "Non-significatif"),
       pch = 16,
       col = c("red", "black"))
dev.off()

# Construction des figures d'ACP #################

# Calcul de l'ACP avec prcomp()
# Utilisation de scale. = TRUE pour centrer et réduire les variables
res_acp_trans <- prcomp(translation_genes[,c(7,8,9,10,11,12)], scale. = TRUE)
# Résumé des résultats de l'ACP et affichage des composantes principales
summary(res_acp_trans)
print(res_acp_trans)
# Explication fournie par la variance de chaque composante principale
#table_acp_trans <- factoextra::get_eig(res_acp_trans)
# Visualisation des cercles des variables dans le plan factoriel de l'ACP
#png("reports/PCA_translation_genes.png", width = 800, height = 600)
#fviz_pca_var(res_acp_trans,
#             col.var = "cos2",
#             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#             repel = TRUE)
#dev.off()

# Construction du Volcano plot ###################

# Création de l'image du Volcano plot
png("reports/Volcano_plot_translation_genes.png", width = 800, height = 600)
# Volcano plot: Log2FoldChange en fonction de - ln de la p-value ajustée
plot(x = resultat_final$log2FoldChange,
     y = -log10(resultat_final$padj),
     main = "Volcano-plot, translation genes",
     ylim = c(0, 20),
     cex = 0.5,
     xlab = bquote(~log[2] ~ fold ~ change),
     ylab = bquote(~-log[10] ~ Q ~ value))
abline(h = 0, col = "black")
dev.off()
