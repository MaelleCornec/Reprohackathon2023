
#Imports nécessaires
library(ggplot2)
library(dplyr)
library(BiocManager)
library(DESeq2)


#Import des données
Full_DATA <- read.table("counts.txt", header = TRUE, skip=1)
Full_DATA$Geneid <- substring(Full_DATA$Geneid, first = 6)

#Sélection des colonnes d'expression des gènes et filtration des valeurs NA
colonnes_selectionnees_Full_DATA <- Full_DATA[, c(7,8,9,10,11,12)]
colonnes_selectionnees_Full_DATA<-na.omit(colonnes_selectionnees_Full_DATA)

#Création d'un objet dds avec DESeq2
dds_full <- DESeqDataSetFromMatrix(countData = Full_DATA[c(7,8,9,10,11,12)],
                                   colData = DataFrame(condition =c("Treatment","Treatment","Treatment",
                                                                    "Control","Control","Control")),
                                   design = ~ condition)

# Normalisation DESeq2
dds_full <- DESeq(dds_full)

res1 <- results(dds_full)

#Passage de BaseMean au log2
res1$baseMean <- log2(res1$baseMean)

#Affichage de la MA-plot: Moyenne d'expression en fonction de Log2FoldChange
plotMA(res1, main = "MA-plot with DESeq2", ylim = c(-4, 4), x="Log2(mean of normalized counts)")

#Plot classique du MA-plot sans utiliser la fonction MA-plot
plot(x=res1$baseMean, y=res1$log2FoldChange, main="MA-plot, all genes", ylim = c(-6, 5))
points(res1$baseMean, cex=0.30, pch=10, res1$log2FoldChange, col = ifelse(res1$pvalue < 0.05, "red", "grey"))

abline(h = 0, col = "black")

#Volcano plot: Log2FoldChange en fonction de la p-value
plot(x=res1$log2FoldChange, y=-log10(res1$pvalue), main="Volcano-plot, all genes", ylim = c(0, 1), cex=0.01)

abline(h = 0, col = "black")

#***************************************************************************************
#***************************************************************************************
#*#***************************************************************************************
#*#***************************************************************************************

#Avec les gènes liés à la traduction
#Lecture des noms des gènes liés à la traduction
noms_genes <- read.csv("Gene_Names_1col_true.csv", header = TRUE, sep = ";")

#Extraction des lignes de res1 qui correspondent aux gènes de la traduction
indices <- which(Full_DATA$Geneid %in% noms_genes$Name)
resultat_final <- res1[indices, c("log2FoldChange", "baseMean", "pvalue")]

#Plot classique du MA-plot sans utiliser la fonction MA-plot
plot(x=resultat_final$baseMean, y=resultat_final$log2FoldChange, main="MA-plot, translation genes", ylim = c(-6, 5))
points(resultat_final$baseMean, cex=0.30, pch=10, resultat_final$log2FoldChange, col = ifelse(resultat_final$pvalue < 0.05, "red", "grey"))

abline(h = 0, col = "black")

#Volcano plot: Log2FoldChange en fonction de la p-value
plot(x=resultat_final$log2FoldChange, y=resultat_final$pvalue, main="Volcano-plot, translation genes", ylim = c(0, 1), cex=0.5)
abline(h = 0, col = "black")




