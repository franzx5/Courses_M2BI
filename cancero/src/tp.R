#!/usr/bin/Rscript

#Cours: Methodologie Bioinfo ~
# Analyse genomique intégrée des adénomes hépatocellulaires
#Prof: Stephane Caruso / Eric Letouzé
#Université Paris Diderot ~ Master Biologie Informatique M2

#Author: AKE Franz-Arnold ~ Pastelle
#Date: 05/19/2018
#Last update : ?

#Set Working directory
setwd("~/Master/Courses_M2BI/cancero/src")


#imported libraries
library(ConsensusClusterPlus)
library(FactoMineR)


#**********************************************************************
#1 - Classification non supervisée des AHC sur la base du transcriptome
#**********************************************************************

#Q1: Reading data
load("../data/AHC_data/expression_matrix_31T_3N_samples.RData")
load("../data/AHC_data/Clinical_annotations.RData")
load("../data/AHC_data/expression_array_annotations.RData")

cat("matrice de taille", dim(exp))
cat("\nla puce contient", dim(exp)[1], "sondes")

#Q2 : selection des 500 sondes les plus variables
top_probes <- tail(sort(apply(exp, 1, sd)),500)


#Q3/Substract Mean values
mat <- exp[names(top_probes),]
cat("table mat de dimension", dim(mat))
mat <- mat - rowMeans(mat) #substract the mean value from the data

#Q4 ~ Hierarchical clustering of the datas
#Computing distance
dist_mat_euclidean <- dist(t(mat), "euclidean")
dist_mat_manhattan <- dist(t(mat), "manhattan")
dist_mat_maximum <- dist(t(mat), "maximum")

#compute clustering
hc_euclidean <- hclust(dist_mat_euclidean, method = "ward.D")
hc_manhattan <- hclust(dist_mat_manhattan, method = "ward.D")
hc_maximum <- hclust(dist_mat_maximum, method = "ward.D") 

#Visualisation
plot(hc_euclidean)
plot(hc_manhattan)
plot(hc_maximum)

#Q5 : Consensus Clustering
pdf("../results/consensus_cluster_figure.pdf")
consclust <- ConsensusClusterPlus(mat, clusterAlg = "hc", distance = "euclidean",
                     innerLinkage = "ward.D", finalLinkage = "ward.D", maxK = 10,
                     pFeature = 0.8, pItem = 0.8)
dev.off()

#Q6 Annotations works
# annot = cbind(annot, rep(NA, dim(annot)[1]))
# names(annot[,13]) = c("expGroup2")
# #group_table <- consclust[[6]]$consensusClass
#annot[match(names(group_table), annot[,1]), "expGroup"]#

#Q7 ~ 
table(annot[c("expGroup","Molecular.group")])

#Q8~ Principal Component Analysis
mat.PCA = PCA(t(mat), graph = FALSE)

#Visualisation
plot.PCA(mat.PCA, axes = c(1,2), choix = "ind", habillage = annot$expGroup)




#**********************************************************
#2 / identification du/des genes drivers dans chaque groupe
#**********************************************************

#Q1 ~ Reading exome
exom <- read.delim("../data/AHC_data/WES_mutations.txt")

#Q2

#Q3
exom$groupes <- annot[match(exom$Sample.ID, annot$Sample.ID), "Molecular.group"]

#Q4 
tail(sort(table(exom$Gene.symbol)))

#Q5
exom_gene_drivers <- exom$groupes[which(exom$Gene.symbol == c("CTNNB1", "HNF1A"))]
table(exom_gene_drivers)



#**************************************************************
#3 - Recherche du mécanismem oncogénique impliqué dans les UHCA
#***************************************************************

#1/
exp <- exp[which(rowMeans(exp) >= 3),]

#3/
source("../data/AHC_data/Functions_for_supervised_analysis.R") #loading some add_functions

#4/ 
refsamp.input <- c("CHC934N", "CHC926N", "CHC469N")
gpsamp.input <- c("CHC603T", "CHC605T", "CHC950T", "CHC1317T")
#Diff expression analysis
ExpCompare(exp, refsamp = refsam.input, gpsamp = gpsamp.input, featureTab = gene_table)










