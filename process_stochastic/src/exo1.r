#!/usr/bin/Rscript

#Author: AKE Franz-Arnold
#Cours: Simulation en Biologie

#Exercice 1 (moodle M2BI  ~ Simulation en biologie)
#implementer une simulation de l'aiguille de Buffon

#Library imported
require(MASS)

#**********
#Functions
#**********
Estime_pi <- function(n, r) {
  #Function to estimate pi value...
  x1 <- runif(n, min = 0, max=1) #generer aleatoire un nombre unif entre 0 et 1
  theta <- runif(n, min = -pi/2, max = pi/2)
  x2 <- x1 + sin(theta) * r
  nsucces <- sum(x2 > 1 | x2 < 0) #decompte du nombre de succès
  return ((2 * n / nsucces) * r)
}

Mse_pi <- function(n, r, nrep){
  #function to compute Mean square error on the number of n throws...
  res <- rep(NA, nrep)
  res <- replicate(nrep, Estime_pi(n, r))
  return((sum(res - mean(res))) / length(res))
}

Mse_pi_variation <- function(n, nrep, r.values){
  #function to compute Mse_variation giving r param variation (vector)
  mse.var <- rep(NA, length(r.values))
  i <- 1
  for (r in r.values){
    mse.var[i] <- Mse_pi(n, r, nrep)
    i <- i + 1
  }
  return(mse.var)
}


#__Main__#
#*********
n <- 100
r.values <- seq(0.1,1,0.05)
nrep <- 10000
mse.variations <- Mse_pi_variation(n, nrep, r.values)
plot(mse.variations~r.values, ty = "l", col = "red", main = "Mse_values depending R parameter") 

