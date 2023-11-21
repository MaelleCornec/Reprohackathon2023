#Importation des librairies
library(DESeq2) #pour l'analyse différentielle de l'expression génique
library(dplyr) #pour la manipulation des données
library(ggplot2) #pour l'affichage

#Chargement des données de comptage
##Stockage des chemins vers les fichiers de comptage dans vecteur path. 
paths <- c(snakemake@input[[1]],
           snakemake@input[[2]],
           snakemake@input[[3]],
           snakemake@input[[4]],
           snakemake@input[[5]],
           snakemake@input[[6]])

#Lecture des fichiers(read.table) et stockage des données de comptage dans une liste ensembleCounts pour chaque échantillon
ensembleCounts = c()
for(i in 1:6){
  path = paths[i]
  ensembleCounts[[i]] <- read.table(path, header = T, skip = 1, sep = '\t')
}

#Création d'une table de comptage à partir des données lues.
##en colonne : les 6 échantillons. 
##en ligne : les gènes. Les noms des lignes sont les noms des gènes extraits de l'échantillon 1.
tableCounts <- data.frame(sample1 = ensembleCounts[[1]]$mapped.sort_SRR10379721.bam,
                          sample2 = ensembleCounts[[2]]$mapped.sort_SRR10379722.bam, 
                          sample3 = ensembleCounts[[3]]$mapped.sort_SRR10379723.bam,
                          sample4 = ensembleCounts[[4]]$mapped.sort_SRR10379724.bam, 
                          sample5 = ensembleCounts[[5]]$mapped.sort_SRR10379725.bam,
                          sample6 = ensembleCounts[[6]]$mapped.sort_SRR10379726.bam)
rownames(tableCounts) <- ensembleCounts[[1]]$Geneid

#Création d'un objet DESeq2 à partir des données de comptage
samples_column <- c(colnames(tableCounts))
type_column <- c("persister","persister","persister","control","control","control") #les 6 échantillons
type_tab <- data.frame(Type = type_column)
rownames(type_tab) <- samples_column

dds <- DESeqDataSetFromMatrix(countData = tableCounts,
                              colData = type_tab,
                              design = ~ Type)
dds$Type <- relevel(dds$Type, ref = "control") 

#Analyse différentielle
dds_post <- DESeq(dds)
res <- results(dds_post, alpha = 0.05) #stockage des résultats 


# MA-plot pour l'ensemble du dataset
dev.size()
diff.df %>% ggplot(aes(x = baseMean, y = log2FoldChange, col = padj < 0.1)) + #coloration des points selon la significativité stats
  geom_point() + 
  theme_classic() +
  scale_color_manual(values = c("grey50", "red")) +
  geom_hline(yintercept = 0, linetype = "dashed") + #ligne en pointillés indiquant la ligne de référence 
  ggtitle("MA-plot of complete RNA-seq dataset") +
  xlab("Mean of normalized counts")
ggsave(snakemake@output[[2]]) #sauvegarde des résultats
dev.off()




