#!/usr/bin/Rscript

#Author: AKE Franz-Arnold
#Cours: Simulation en Biologie

#Exercice 2 (moodle M2BI  ~ Simulation en biologie)
#Comparaison de l'efficacite de methode d'integration Monte Carlo

#Setting working directory
setwd("~/Master/Courses_M2BI/process_stochastic/src/")

#import library
require(MASS)

#********************
#function g à integrer
G <- function(x){
  return((exp(x)-1)/(exp(1)-1))
}
#Graphe
valX <- seq(from = 0, to = 2, len = 1000)
valY <- G(valX)
plot(valX, valY, xlab = "x", ylab = "g(x)", type = "l", col = "red") # plotting ...

#integrate I
I <- (exp(2)-3) / (exp(1)-1)
print(I)
print(integrate(G,0,2))



#**********************************
Mc_bw <- function(n, a , b){
  #tirage monte carlo blanc ou noir
  # n ~ nombre de tirages
  # a , b ~ intervalle de selection pour le tirages des valsX
  
  u <- runif(n, min = a, max = b)
  m <- max(G(seq(a, b, len = 1000)))
  v <- runif(n, min = 0, max = m)
  nsucces <- sum(v < G(u))
  return(m * (b-a) * (nsucces/n))
}


#***********************************
Mc_simple <- function(n, a ,b){
  #tirage monte carlo simple
  # n ~ nombre de tirages
  # a , b ~ intervalle de selection pour le tirages des valsX
  
  valX <- runif(n, min = a, max = b)
  h <- (b-a) * G(valX)
  return(mean(h))
}

#***********************************

#***********************************
Mc_importance_beta = function(n, a){
  # n ~ nombre de tirages
  # a ~ coefficient multiplicateur
  #Monte_carlo_importance_beta
  #la loi beta generant des nombres aleatoire sur le support
  #[0-1], necessaire de multiplier par un entier a(~2) > 0 pour avoir des nombres
  #sur le support [0-2]
  valX = a * rbeta(n, 2, 1)  #param alfa = 2 & beta = 1
  fX = (1/a) * dbeta(valX/a, 2, 1)
  return(mean(G(valX)/fX))
}


Compute_MSE <- function(liste.values){
  return(sum(liste.values - 2.55) / length(liste.values))
}

#************************************************************
#Compute Mse Monte Carlo performances
output = c()
for (i in seq(10, 10000, 700)){
  res <- Compute_MSE(replicate(100, Mc_importance_beta(i, 2)))
  print(abs(res))
  output = append(output, abs(res))
}
plot(output, ty = "l", col = "red", main = "Evolution du MSE ~ n tirages", xlab = "n",
     ylab = "MSE_value")


#Comparaison
#***********
n = 100
nrep = 1000000
mse.MC.BW = Compute_MSE(replicate(nrep, Mc_bw(n, 0, 2)))
mse.MC.simple = Compute_MSE(replicate(nrep, Mc_simple(n, 0, 2)))
mse.MC.importance = Compute_MSE(replicate(nrep, Mc_importance_beta(n, 2)))

print(abs(c(mse.MC.BW, mse.MC.simple, mse.MC.importance)))
barplot(abs(c(mse.MC.BW, mse.MC.simple, mse.MC.importance)), col = cm.colors(256))



























