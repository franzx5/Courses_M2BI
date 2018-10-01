#!/usr/bin/Rscript

#Author: AKE Franz-Arnold
#Cours: Simulation en Biologie

#Exercice 2 (moodle M2BI  ~ Simulation en biologie)
#Comparaison de l'efficacite de methode d'integration Monte Carlo

#import library
require(MASS)

#********************
#function g à integrer
G <- function(x){
  return((exp(x)-1)/(exp(1)-1))
}
#testing ...!
valX <- seq(from = 0, to = 2, len = 1000)
valY <- G(valX)
plot(valX, valY, xlab = "x", ylab = "g(y)", type = "l", col = "red") # plotting ...
I <- (exp(2)-3) / (exp(1)-1)
print(I)
print(integrate(G,0,2))


#**********************************
Mc_bw <- function(n, a , b){
  #tirage monte carlo "blanc ou noir"
  valX <- runif(n, min = a, max = b)
  m <- G(b)
  valY <- runif(n, min = 0, max = m)
  nsucces <- sum(valY < G(valX))
  return(m * (b-a) * (nsucces/n))
}
#testing ...!
n = 100
Mc_bw(n, 0, 2)


#******************
Mc_simple <- function(n, a ,b){
  #tirage monte carlo simple
  valX <- runif(n, min = a, max = b)
  h <- (b-a) * G(valX)
  return(mean(h))
}
#testing ...!
n <- 100
Mc_simple(n, 0, 2)


#***************************
Mc_importance_beta = function(n, a){
  #Monte_carlo_importance_beta
  #la loi beta generant des nombres aleatoire sur le support
  #[0-1], necessaire de multiplier par un entier a(~2) pour avoir des nombres
  #sur le support [0-2]
  valX = a * rbeta(n, 2, 1)
  fX = (1/a) * dbeta(valX/a, 2, 1)
  return(mean(G(valX)/fX))
}
#testing ...!
n = 10000
Mc_importance_beta(n, 2)


Compute_MSE <- function(liste.values){
  return(sum(liste.values - mean(liste.values)) / length(liste.values))
}


#*****************************
#Compare Mse
k.n <- 100
k.nrep <- 10
mse.MC.BW <- Compute_MSE(replicate(nrep, Mc_bw(n, 0,2)))
mse.MC.simple <- Compute_MSE(replicate(nrep, Mc_simple(n, 0, 2)))
mse.MC.importance <- Compute_MSE(replicate(nrep, Mc_importance_beta(n, 2)))




























