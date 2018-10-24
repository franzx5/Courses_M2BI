
#Polymorphisme & Human Genom
#Date : 12/10/2018
#Prof : J. Noirel
#Source tp (http://jnoirel.fr/teaching/M2BIB/)

#Setting working directories
setwd("~/Master/Courses_M2BI/polymorphisme_human_genom/src/")


#importing packages
source("https://www.bioconductor.org/biocLite.R")
biocLite("SNPlocs.Hsapiens.dbSNP141.GRCh38")
biocLite("biomaRt")

#libraries
library("SNPlocs.Hsapiens.dbSNP141.GRCh38")
library(biomaRt)
library(tidyverse)
library(dplyr)

#readind data
snpmart = useMart(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")

#research
listFilters(snpmart) %>% filter(str_detect(description, "Chr"))

#Loading on biomart
snps = getBM(
  attributes = c('refsnp_id', "minor_allele_freq"),
  filters = c('chr_name', 'start', "end"),
  values = list(8, 100000, 1000000),
  mart = snpmart
)


#graphs
snps = snps %>% filter(! is.na(minor_allele_freq) & minor_allele_freq > 0.05)

#hist
hist(snps$minor_allele_freq)

set.seed(1)
size=20
m = 10
genotypes = matrix(rbinom(size * m, size * 2, prob =0.5), nrow = size)


# missing = matrix(sample(c(NA, 0), size = 

###############################################
#Quality control
#logistic simulation
x = rnorm(1000, mean = 0, sd = 0)
s = 3 + (x/5)
p = 1/(1+ exp(-S))



phneotypes = sapply(P, FUN = function(p)
  sample(c("Case", "Control"), size = 1, prob=c(p, 1-p), replace = TRUE))

phenotypes = phenotypes == "Case"

glm(phenotypes ~ x, family = "binomial")
 




dat = read_tsv("starified.dat")

dat$Phenotype = dat$Phenotype == "Case"
#dat$Phentotype = as.nuemric(dat$Phenotype == "Case")


model = glm(Phenotype ~ Genotype, data = dat, family = "binomial")
model


summary


#Stratification
library("snpStats")
library(FactoMineR)

bed = 'data/hapmap1/small.bed'
bim = 'data.hapmap1/small.bim'
fam = 'data/hapmap1/small.fam'
dat = read.plink(bed, bim, fam)

genotypes = as.numeric(dat$genotypes)
genotypes[genotypes == 0] = NA
genotypes = matrix(genotypes, nrow = nrow(dat$fam))
genotypes = genotypes -1


#PCA

#Genetic Association Study

dat$fam$PC1 = pca$x[, 1]


association = snp.rhs.estimates(affected ~ PC1 + PC2, #SNP is implicity given
                                family = "binomial",
                                data = dat$fam,
                                snp.data = dat$genotypes)

#Zscore_calcul
zscores = map_dbl(association, ~.$beta/sqrt(.$Var.beta))
pvalues = pnorm(abs(zscores), lower = FALSE)
p.adjusted = p.adjust(pvalues, method = "fdr")
sort(padjusted)[1:10]

























