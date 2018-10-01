#!/usr/bin/Rscript
#Cours: Analyses génomiques
#Prof: Karine Audouze 
#Université Paris Diderot ~ Master Biologie Informatique M2

#Author: AKE Franz-Arnold
#Date: 25/09/2018

#Imported libraries
#source("https://bioconductor.org/biocLite.R")
require(genefilter)
require(annotate)
library(affy)
library(hgu133a.db)
require(annotate)
#Exo 3: Normalization	ans	statistical	analysis	of	microarray	data	exercise
    
#Introduction
  #Question
#1 - Which	genes	are	expressed	differently	in	HIV	infected	patients ?- the	disease	effect	

# 2. Do	the	cell	types	respond	in	different	ways	to	the	infection,		i.e.	there	is	any
# difference	between	which	genes	that	we	find	differentially	expressed	in	the	two	cell
# types ?	(healthy	versus	infected	T-helper	cells	and	healthy	versus	infected	CTLs)	-
# the	interaction	effect

# 3. Are	they	any	significant	differential	gene	expression	between	the	two	cell	types ?	(Thelper
# cells	versus	CTLs)- the	cell	type	effect



#Loading and Normalizing Data
#****************************
HIV.affyb <- ReadAffy(filenames = sort(list.celfiles("../data/cellfiles", full = TRUE))) #load data
# rma ~ Robust multi array method giving porbes values into expression index values
ei = exprs(rma(HIV.affyb))
#index-lines ~ probe set
#index-cols ~ sample

#...describing data structures (...Can also be done by reading an input design file)
design <- data.frame(HIV=rep(c("ctrl","HIV", "HIV","ctrl"),each=5), cell=
  rep(c("CTL","Th"),each=10), row.names=colnames(ei))


#splitting plot windows
par(mfrow = c(1,2))
#density plot of normalized data
plot ( density ( ei[,1], na.rm=TRUE ), xlim = c(0,16), ylim = c(0,0.3), main
   = "Normalized", xlab = "Intensity" )

#..For all the sample ( make a loop on all the sample with a different color)
for(i in 2:ncol( ei ) ) {
  points( density( ei[,i], na.rm = TRUE ), type = "l", col = rainbow(20)[i] )
}
#Normalisation n'est pas optimale, l'on peut passer par d'autres algo de normalisation

#Exemple ~ for the non normalized dataset
hist(HIV.affyb, main="Un-normalized")

#Make an pdf file output
dev.copy2pdf(file = "../doc/DensityPlots.pdf")
dev.off()

#Q1 - What	do	you	see	in	this	plot	and	does	the	data	look	normalized?
# the courbes are better superposed compared to the non normalized plot dataset,
# the means values are very close

#Q2 - Could the data could be more comparable
# yes, the density values of each plot are not very equals, we can change the normalization method,
# (we can test the QSpline normalization method) to improve that



#Singular	Value	Decomposition	and	outlier	assessment
#***************************************************
# An	SVD	analysis	kind	of	breaks	your	data	down	into	X	
# components	(where	X	=	number	of	samples).	The	first	component	contains	the	most	
# variation,	and	the	last	the	least	(none).

m <- ncol(ei)
X <- ei - rowMeans(ei) #(val_actuelles - val_moyennes_samples)
svd <- svd (X)
V <- svd$v
S <- matrix(0 , m, m)
diag(S)<-svd$d

pnames <- paste(design$HIV, design$cell, rep(1:5, 4), sep="-")
cat_col <- rep(c("blue","darkgreen","red","orange"), each=5)

par(mfrow = c(1,2))
#How the fig image inform on the 
image( S%*%t(V), main="Singular Value Decomposition",
       ylab = "Samples ( 1-20 )", xlab = "The 20 Components" )

plot(V[,1], V[,2], col = "white", main = "Singular Value Decomposition",
     xlab = "First Component", ylab = "Second Component")
text(V[,1], V[,2], pnames, col = cat_col )
plot(V[,1], V[,3], col = "white", main = "Singular Value Decomposition",
     xlab = "First Component", ylab = "Third Component")
text(V[,1], V[,3], pnames, col = cat_col )
dev.print(device = pdf, file = "../results/SVD.pdf" )
dev.off()


#Q3:	Do	you	see	any	systematic	patterns?	And	what	do	you	think	this	means?
# we can observe a particular structure of two clusters in our dataset, based on the observation of the SVD.
# an cluster is constituted by control sample et the second by non-control sample.
# this observation means that we have a close expression of genes of certains sample,
# (HIV-Th, ctrl-Th, ctrl-CTRL)  compared to (HIV-CTL).


#Q4: Which	two	components	contain	the	most	information	in	regard	to	your	categories?
#the first and third component are most informative ...


#The T-test
#**********
#t-test for each cell type
t.pval.CTL = rowttests(ei[,design$cell=="CTL"],design$HIV[design$cell=="CTL"])$p.value
t.pval.Th = rowttests(ei[,design$cell=="Th"],design$HIV[design$cell=="Th"])$p.value
names(t.pval.Th) = names(t.pval.CTL) = rownames(ei) 

#... alternative way to group the data
rowttests(ei[,1:10],	factor(c(1,1,1,1,1,2,2,2,2,2)))

#Get the genes name
AffyID <- rownames( ei )
GeneNames <- getSYMBOL( AffyID , "hgu133a" )

#names of the most significatives exprimed genes
GeneNames[names(sort(t.pval.CTL)[1:10])]
GeneNames[names(sort(t.pval.Th)[1:10])]

#for knowing sens of expression (neg or positive?) --> log2 ~ fold change
fc.CTL<- rowMeans(ei[,design$cell=="CTL" & design$HIV=="HIV"]) -
  rowMeans(ei[,design$cell=="CTL" & design$HIV=="ctrl"])
fc.Th <- rowMeans(ei[,design$cell=="Th" & design$HIV=="HIV"]) -
  rowMeans(ei[,design$cell=="Th" & design$HIV=="ctrl"])

#Take	a	look	at	the	fold	changes	of	the	30	most	significant	genes:
fc.CTL[ order(t.pval.CTL) ][1:30]
fc.Th[ order(t.pval.Th) ][1:30]

# Visualisation in a heatmap
Th_Top <- ei[ order(t.pval.Th)[1:30], design$cell=="Th"]
rownames(Th_Top)<-GeneNames[order(t.pval.Th)[1:30]]
heatmap(Th_Top , main = "Th-cells", scale="none", mar=c(8,5))

CTL_Top<-ei[ order(t.pval.CTL)[1:30], design$cell=="CTL"]
rownames(CTL_Top)<-GeneNames[order(t.pval.CTL)[1:30]]
heatmap(CTL_Top, main = "CTL-cells", scale= "none", mar=c(8,5))


#Q5 -	What	are	we	testing	in	these	experiments?
x <- rownames(Th_Top)
y <- rownames(CTL_Top)
  
# Q6:	How	many	genes	have	we	tested?	hint:	length(ei[,1])

# Q7:	At	which	p-value	would	you	apply	the	cutoff	if	you	were	to	use	the	Bonferroni	
# correction	at	a	0.05	level?

# Q8:	How	many	genes	would	we	trust	in	each	of	the	two	experiments	using	this	cutoff?	hint:	
#   min(t.pval.Th <= the bonferoni cutoff 
      
# Q10:	If	you	choose	to	use	the	Benjamini-Hochberg	cutoff	in	the	"Th"	test,	how	many	genes	
# will	you	then	find	to	be	significant?




#Volcano Plots
#*************
random_lables<-sample(design$HIV[design$cell=="CTL"])

# log2	foldchange	and the	p-value
perm_t.pval.CTL <- rowttests(ei[,design$cell=="CTL"], random_lables)$p.value
perm_fc.CTL<- rowMeans(ei[,design$cell=="CTL" & random_lables=="HIV"]) -
  rowMeans(ei[,design$cell=="CTL" & random_lables=="ctrl"])

#Volcano_Plots visualisation
plot(fc.CTL, t.pval.CTL, main = "Volcano Plot\nCTL", log = "y",
     xlab = "M (log2 fold change)", ylab = "p-value", pch = 20, col = "blue")
points(perm_fc.CTL, perm_t.pval.CTL, type = "p", pch = 20, col = "red")
legend("topleft", col=c("red","blue"), legend=c("perm", "real"), pch=20, bg="white")

# Q11:	Do	you	think	that	your	cutoff	from	Q5	and	Q7	were	reasonable?	
# Normally	one	should	make	a	couple	of	Volcano	plots	using	balanced	randomization's	in	
# order	to decide	where	to	put	your	cutoff/	or	estimate	FDR,	but	for	now	one	is	fine.	

# Q12:	Where	do	you	think	the	best	cutoff	would	be?

# Q13:	Do	you	think	that	these	plots/calculations	could	somehow	be	used	to	estimate	the	
# False	Discovery	Rate?	How	would	you	do	that?

# Q14:	Are	we	now	able	to	answer	the	first	biological	question?	
# 1.	Which	genes	are	expressed	differently	between	HIV-infected	and	control	cells?


#The	Two-Way	ANOVA:
#*******************

anova.2F <- function(vector, F1, F2, design.df){
  y<-anova(aov(vector ~ design.df[,F1]*design.df[,F2]))$Pr[1:3]
  return(y)
}

#ANOVA.pval <- apply(ei, 1, anova.2F, F1="HIV", F2="cell", design)

#Extract	the	three	lists	of	p-values:
HIV.pval <- ANOVA.pval[1,]
CellType.pval <- ANOVA.pval[2,]
Interaction.pval <- ANOVA.pval[3,]






