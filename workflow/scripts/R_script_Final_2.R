
library(ggplot2)
library(dplyr)
library(BiocManager)
library(DESeq2)
library(factoextra)

Full_DATA <- read.table("genes/counts.txt", header = TRUE, skip=1)
Full_DATA$Geneid <- substring(Full_DATA$Geneid, first = 6)

colonnes_selectionnees_Full_DATA <- Full_DATA[, c(7,8,9,10,11,12)]

colonnes_selectionnees_Full_DATA<-na.omit(colonnes_selectionnees_Full_DATA)

dds_full <- DESeqDataSetFromMatrix(countData = Full_DATA[c(7,8,9,10,11,12)],
                                   colData = DataFrame(condition =c("Treatment","Treatment","Treatment",
                                                                    "Control","Control","Control")),
                                   design = ~ condition)

# Normalisation DESeq2
dds_full <- DESeq(dds_full)

res1 <- results(dds_full)

#On passe la moyenne des expressions au log2
res1$baseMean <- log2(res1$baseMean)

#Plot classique du MA-plot sans utiliser la fonction MA-plot
png("MA-plot_all_genes.png", width = 800, height = 600)

plot(x=res1$baseMean, y=res1$log2FoldChange, main="MA-plot, all genes", ylim = c(-6, 5), cex=0.9)
points(res1$baseMean, cex=0.50, pch=16, res1$log2FoldChange, col = ifelse(res1$pvalue < 0.05, "red", "black"))

abline(h = 0, col = "black")
legend("bottomleft", legend = c("Significatif", "Non-significatif"), pch = 16, col = c("red", "black"))

dev.off()

png("PCA_all_genes.png", width = 800, height = 600)
# Calcul de l'ACP avec prcomp()
resultat_acp <- prcomp(colonnes_selectionnees_Full_DATA, scale. = TRUE)  # Utilisation de scale. = TRUE pour centrer et réduire les variables

# Résumé des résultats de l'ACP
summary(resultat_acp)

# Affichage des composantes principales
print(resultat_acp)

# Visualisation des cercles des variables dans le plan factoriel de l'ACP
fviz_pca_var(resultat_acp, col.var = "black") 
dev.off()

#Volcano plot: Log2FoldChange en fonction de la p-value
png("Volcano_plot_all_genes.png", width = 800, height = 600)

plot(x=res1$log2FoldChange, y=-log10(res1$pvalue), main="Volcano-plot, all genes", ylim = c(0, 1), cex=0.01)

abline(h = 0, col = "black")
dev.off()

noms_genes <- read.csv("Gene_Names_1col_true.csv", header = TRUE, sep = ";")
indices <- which(Full_DATA$Geneid %in% noms_genes$Name)
resultat_final <- res1[indices, c("log2FoldChange", "baseMean", "pvalue")]

GENES_TRANSLATION_DATA <- Full_DATA[indices, c(1,2,3,4,5,6,7,8,9,10,11,12)]

#Plot classique
png("MA-plot_translation_genes.png", width = 800, height = 600)

#Définition des 6 gènes particuliers
genes_to_label <- c("SAOUHSC_01236", "SAOUHSC_02489", "SAOUHSC_00475", "SAOUHSC_01234", "SAOUHSC_01786", "SAOUHSC_01246")
genes_code <- c("frr", "infA", "pth", "tsf", "infC", "infB")

index <- c()

# Création de l'index pour les 6 gènes particuliers
for (i in seq_along(genes_to_label)) {
  gene <- genes_to_label[i]
  index <- append(index, which(GENES_TRANSLATION_DATA$Geneid == gene)) 
}

# Création de l'index pour les points à encadrer en bleu (liés à l'aminoacyl ARNt synthétase)
start_index <- which(noms_genes$Name == "SAOUHSC_T0001")  # Trouver l'indice où Geneid == "SAOUHSC_T0001"
#liste des gènes liés à la synthèse de l'aminoacyl ARNt 
noms_genes_ARNt<-noms_genes$Name[start_index:162]
#liste des index dans GENES_TRANSLATION_DATA
index_to_encircle<- c()
for (i in seq_along(noms_genes_ARNt)) {
  gene <- noms_genes_ARNt[i]
  index_to_encircle <- append(index_to_encircle, which(GENES_TRANSLATION_DATA$Geneid == gene)) 
}

# Le plot
plot(
  x = resultat_final$baseMean,
  y = resultat_final$log2FoldChange,
  main = "MA-plot, translation genes",
  ylim = c(-6, 5),
  cex=0.9
)

#Coloration des points en fonction de la valeur de la p-value
points(
  resultat_final$baseMean,
  cex = 0.50,
  pch = 16,
  resultat_final$log2FoldChange,
  col = ifelse(resultat_final$pvalue < 0.05, "red", "black"),
)
# Encadrer les points correspondants aux gènes liés à la synthèse de l'aminoacyl ARNt synthétase en bleu
x_vals <- resultat_final$baseMean[index_to_encircle]
y_vals <- resultat_final$log2FoldChange[index_to_encircle]

# Calculer les dimensions du carré
square_size <- 0.2  # Choisir la taille du carré

# Calculer les coins du carré pour chaque point
x1 <- x_vals - square_size/2
x2 <- x_vals + square_size/2
y1 <- y_vals - square_size/2
y2 <- y_vals + square_size/2

for (i in seq_along(index_to_encircle)) {
  rect(
    xleft = x1[i],
    ybottom = y1[i],
    xright = x2[i],
    ytop = y2[i],
    border = "green",
    lwd = 0.5  # Choisir l'épaisseur de la ligne du carré
  )
}

#tracé de la ligne noire en y=0
abline(h = 0, col = "black")

#Ecriture des textes pour les 6 gènes particuliers
text(
  x = resultat_final$baseMean[index],
  y = resultat_final$log2FoldChange[index],
  pos = 4,
  labels = genes_code,
  cex = 1
)

# Légende pour les points significatifs et non-significatifs
legend("bottomleft", legend = c("Significatif", "Non-significatif", "AA-tRNA synthetases"), pch = 16, col = c("red", "black","green"))

dev.off()

# Calcul de l'ACP avec prcomp()
resultat_acp_translation <- prcomp(GENES_TRANSLATION_DATA[,c(7,8,9,10,11,12)], scale. = TRUE)  # Utilisation de scale. = TRUE pour centrer et réduire les variables

# Résumé des résultats de l'ACP
summary(resultat_acp_translation)

# Affichage des composantes principales
print(resultat_acp_translation)

png("PCA_translation_genes.png", width = 800, height = 600)

# Visualisation des cercles des variables dans le plan factoriel de l'ACP
fviz_pca_var(resultat_acp_translation, col.var = "black") 
dev.off()

#Volcano plot: Log2FoldChange en fonction de la p-value
png("Volcano_plot_translation_genes.png", width = 800, height = 600)

plot(x=resultat_final$log2FoldChange, y=-log10(resultat_final$pvalue), main="Volcano-plot, translation genes", ylim = c(0, 1), cex=0.5)

abline(h = 0, col = "black")

dev.off()
