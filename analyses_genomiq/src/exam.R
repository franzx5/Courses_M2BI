#!/usr/bin/Rscript
#Cours: Analyses génomiques
#Prof: Karine Audouze 
#Université Paris Diderot ~ Master Biologie Informatique M2

#Author: AKE Franz-Arnold ~ Pastelle
#Date: 28/09/2018

#Set current working directory
setwd(dir= "~/Master/Courses_M2BI/analyses_genomiq/src/")

# import libraries
source("https://bioconductor.org/biocLite.R")
require(oligo)
require(annotate)
require(genefilter)
#library(affy)
library(pd.mogene.2.0.st)
library(mogene20sttranscriptcluster.db)
library(hgu133a.db)
library(FactoMineR)
library(factoextra)
library(rmarkdown)
library(ensembldb)
library("RDAVIDWebService")
library(gplots)


#Reading Data CeL (needed mogene library)
data.Cel = read.celfiles(filenames = sort(list.celfiles("../data/CEL-20180928/Archive/", full = TRUE)))
#normalisation des données
ei = exprs(rma(data.Cel))
#Annotate Genes names
data.Cel.ID <- rownames(ei)
GeneNames <- unlist(contents(mogene20sttranscriptclusterGENENAME)[data.Cel.ID])
GeneSymb <- unlist(contents(mogene20sttranscriptclusterSYMBOL)[data.Cel.ID])
GeneEnsembl <- unlist(contents(mogene20sttranscriptclusterENSEMBL)[data.Cel.ID])

#Visualisation
par(mfrow = c(1,2))
plot(density(ei[,1], na.rm = TRUE), xlim = c(0,16), ylim = c(0,0.3), main
     = "Normalized Data", xlab = "Intensity")
for(i in 2:ncol( ei ) ) {
  points( density( ei[,i], na.rm = TRUE ), type = "l", col = rainbow(20)[i] )
}
plot(density(ei[,1:9], na.rm = TRUE), xlim = c(0,16), ylim = c(0,0.3), main
     = "Normalized Data without sample I", xlab = "Intensity")
for(i in 2:ncol( ei )-1) {
  points( density( ei[,i], na.rm = TRUE ), type = "l", col = rainbow(20)[i] )
}
dev.off()

cat("
    Q1:
    We can observe (plot 1) that all the sample are normalized and centred except sample I which have a bad fit
    with others sample (plot 2). It must be an outlier.
")


cat("Q2 : Combien	de	gènes	uniques	sont	présents? ~ Pourquoi	certains	sont	ils	répliqués ? ")
nb.unique.genes =  length(GeneNames[!(duplicated(GeneNames)|duplicated(GeneNames, fromLast=TRUE))])
cat(paste("we have",nb.unique.genes,"unique genes"))
cat("
    we have multiples same genes because the probes match probably
    multiples same transcript of the same gene
")


cat("Q3:  Réaliser	une	PCA 	Y	a	t’il	un	échantillon	‘outlier’ ?	si	oui	lequel.	Dans	ce	cas	ne	pas	
 l’utiliser	pour	la	suite	des	analyses.")
design <- c(rep(c("heal","sick"), each=1, 3), rep(c("heal"),2), rep(c("sick"),2))
design <- data.frame(design)
data.Cel.design <- cbind(design, t(ei))  #concatenate design with ei
rownames(data.Cel.design) = c("A", "J", "B", "C", "D", "E", "F", "G", "H", "I")  #rename samples rows
colnames(ei) = rownames(data.Cel.design)
data.Cel.PCA = PCA(data.Cel.design, ncp=5, quali.sup = 1, graph = F)  #compute PCA

#Visualisation PCA ~ With the Outlier I
#Percentage of cumulative Variance
#barplot(data.Cel.PCA$eig[,1], main = "Eboulis des Eigen Values of ACP")
fviz_screeplot(data.Cel.PCA, addlabels = TRUE, ylim = c(0, 50),
               main = "Evolution de la variance expliquée ~ Dimensions")
# plot.PCA(data.Cel.PCA, axes = c(1,4), choix = "ind", habillage = 1)
fviz_pca_ind(data.Cel.PCA, habillage = data.Cel.design[,1], addEllipses = F)

#Echantillon_Outlier
cat("Oui nous avons un echantillon outlier, en l'occurence l'echantillon I
      Elimination de l'echantillon outlier...")
ei = ei[,1:9]
design <- c(rep(c("heal","sick"), each=1, 3), rep(c("heal"),2), rep(c("sick"),1))
design <- data.frame(design)
data.Cel.design <- cbind(design, t(ei))  #concatenate design with ei
rownames(data.Cel.design) = c("A", "J", "B", "C", "D", "E", "F", "G", "H")  #rename samples rows


cat("Q4: Réaliser	un	volcano	plot.	Sélectionner	un	seuil.	Pourquoi	avez-vous	choisi ce	seuil ?")
#T-tests
t.pval = rowttests(ei,design$design)$p.value
names(t.pval) = rownames(ei)

#Fold_change
fc_cel = rowMeans(ei[, design$design == "heal"]) - rowMeans(ei[, design$design == "sick"])

#Random
random.labels = sample(design$design)
perm_t.pval = rowttests(ei, random.labels)$p.value
perm_fc_cel = rowMeans(ei[, random.labels == "heal"]) - rowMeans(ei[, random.labels == "sick"])

#Visualisation
plot(fc_cel, t.pval, main = "Volcano Plot", log = "y", xlab = "M (log2 fold change)", 
     ylab = "p-value", pch = 20, col = "blue", xlim = c(-2.5,2.5))
points(perm_fc_cel, perm_t.pval, type = "p", pch = 20, col = "red")
legend("topleft", col=c("red","blue"), legend=c("perm", "real"), pch=5, bg="white", cex = 0.5)

abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=0.008, col="black", lty=4, lwd=2.0) 
cat("threshold défini à 0.008")
#why ?

#Selection of significatives genes
ei_top = ei[names(t.pval[t.pval<0.001]),] #we have 37 genes
rownames(ei_top) = GeneSymb[rownames(ei_top)]
cat(paste("we have", dim(ei_top)[1], "genes"))


#Q6 - Realiser un enrichissement biologique (GSEA) ~ dbase KEGG ~ tool DAVID



#Q7 - Compute HeatMap
#visualisation
scaleyellowred <- colorRampPalette(c("white", "green", "blue", "orange", "red"), space = "rgb")(256)
heatmap.2(ei_top, col=scaleyellowred, trace = "none", density.info = "none")









