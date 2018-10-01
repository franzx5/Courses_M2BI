#!/usr/bin/Rscript
#Cours: Analyses génomiques
#Prof: Karine Audouze 
#Université Paris Diderot ~ Master Biologie Informatique M2

#Author: AKE Franz-Arnold ~ Pastelle
#Date: 28/09/2018

# import libraries
source("https://bioconductor.org/biocLite.R")
biocLite("pd.mogene.2.0.st")
biocLite("mogene20sttranscriptcluster.db")
biocLite("RDAVIDWebService")
require(oligo)
require(annotate)
require(genefilter)
#library(affy)
library(pd.mogene.2.0.st)
library(mogene20sttranscriptcluster.db)
library(hgu133a.db)
library(FactoMineR)




#Reading Data CeL (needed mogene library)
data.Cel = read.celfiles(filenames = sort(list.celfiles("../data/CEL-20180928/Archive/", full = TRUE)))

#normalisation des données
ei = exprs(rma(data.Cel))

#Annotate Genes names
data.Cel.ID <- rownames( ei )
GeneNames <- getSYMBOL( data.Cel.ID , "mogene20sttranscriptcluster") 
#How to get "gene_symbol" & "ensembl_id" ?


# 2 /Combien	de	gènes	uniques	sont	présents? ~ Pourquoi	certains	sont	ils	répliqués ?
nb.unique.genes =  length(GeneNames[!(duplicated(GeneNames)|duplicated(GeneNames, fromLast=TRUE))])
#we have multiples same genes because the probes match probably
#multiples same transcript of the genes


# # 3/ - Réaliser	une	PCA 	Y	a	t’il	un	échantillon	‘outlier’ ?	si	oui	lequel.	Dans	ce	cas	ne	pas	
# l’utiliser	pour	la	suite	des	analyses.
design <- c(rep(c("heal","sick"), each=1, 3), rep(c("heal"),2), rep(c("sick"),2))
design <- data.frame(design)

data.Cel.design <- cbind(design, t(ei))  #concatenate design with ei
rownames(data.Cel.design) = c("A", "J", "B", "C", "D", "E", "F", "G", "H", "I")  #rename samples rows
colnames(ei) = rownames(data.Cel.design)
data.Cel.PCA = PCA(data.Cel.design, ncp=5, quali.sup = 1, graph = F)  #compute PCA

#Visualisation PCA
plot.PCA(data.Cel.PCA, axes = c(1,2), choix = "ind", habillage = 1)



# 4- Réaliser	un	volcano	plot.	Sélectionner	un	seuil.	Pourquoi	avez-vous	choisi ce	seuil ?
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
legend("topleft", col=c("red","blue"), legend=c("perm", "real"), pch=10, bg="white")

abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=0.001, col="black", lty=4, lwd=2.0) #threshold défini à 0.001

#Selection of significatives genes
ei_top = ei[names(t.pval[t.pval<0.001]),] #we have 37 genes


#Q6 - Realiser un enrichissement biologique (GSEA) ~ dbase KEGG ~ tool DAVID
#?


#Q7 - Compute HeatMap
#visualisation
heatmap(ei_top, col = cm.colors(256))









