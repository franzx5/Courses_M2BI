#!/usr/bin/Rscript

#Author: AKE Franz-Arnold


#test de Fisher


#Step 1:

table.1 = matrix(NA, nrow =2, ncol = 2)
colnames(table.1) = c("-GO","+GO")
rownames(table.1) = c("Non-sig", "Signif")

d = 13
c = 100 - 13
b = 125 - 13
a = 6000 - (b + c + d)

table.1[1,1] = a
table.1[1,2] = b
table.1[2,1] = c
table.1[2,2] = d

#Test de Fisher
fisher.test(table.1)$p.value
