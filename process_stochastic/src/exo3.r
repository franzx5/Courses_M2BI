#!/usr/bin/Rscript

#Author: AKE Franz-Arnold
#Cours: Simulation en Biologie

#Exercice 3 (moodle M2BI  ~ Simulation en biologie)

#Q1/

my_rpois_one <- function(lambda){
  u = runif(1)
  x = 0
  while(u > ppois(x, lambda)){
    x = x + 1
  }
  return(x)
}

my_rpois <- function(n, lambda){
  ret = rep(NA, n)
  for (i in 1:n){
    ret[i] = my_rpois_one(lambda) #or use replicate(n, my_rpois_one)
  }
  return(ret)
}

lambda = 2.5
n = 10000
res = my_rpois(n, lambda)
plot(table(res)/n)

#visualisation graph
valX = seq(0,12)
valY = dpois(valX, lambda)
points(valX, valY, col = "red")


#Q2/
#**********
my_rdiscret_one <- function() {
  u = runif(1)
  x = c(-3,1.3,7,15.2)
  Fx = cumsum(c(0.1,0.4,0.3,0.2))
  i = 1
  while(u > Fx[i]){
    i = i + 1
  }
  return(x[i])
}


#Q3/
#**********
My_dlaplace <- function(x){
  return(1/2 * exp(- abs(x)))
}

My_plaplace <- function(x){
  if (x<=0)
    return(1/2 * exp(x))
  else
    return(1 - (1/2 * exp(-x)))
}
My_qlaplace <- function(p){
  if (p <= 1/2){
    return(log(2 * p))
  } else if (p > 1/2) {
    return(-log(2*(1-p)))
  }
}
My_rlaplace <- function(n){
  u = runif(n)
  res = c()
  for (i in 1:length(u)){
    res[i] = My_qlaplace(u[i])
  }
  return(res)
}
#visualisation 
truehist(My_rlaplace(1000))
x_vect = seq(-6,8, length=10000)
lines(x_vect, My_dlaplace(x_vect), col = "red")



#Q4 ~ Methode de Rejet ~ Loi normale
#**********************
Meth_rejet <- function(n){
  m = sqrt(2*(exp(1)/pi))
  x1 = My_rlaplace(n)
  x2 = runif(n)
  filtre = which(x2 <= dnorm(x1)/(m*My_dlaplace(x1)))
  return(x1[filtre])
}


#Q5 ~ Algorithme MCMC
#********************
D_target <- function(x){
  output = 0.2 * dnorm(x, -3, 2) + 0.5 * dnorm(x, 0, 1) + 0.3 * dnorm(x, 5, 3)
  return(output)
}

R_prop <- function(x, delta){
  return(runif(1, min = x-delta, max = x+delta))
}

My_rtarget <- function(x0, nstep){
  vecX = rep(NA, length=nstep)
  vecX[1] = x0
  for (i in 1:nstep){
    x_i = vecX[i]
    y = R_prop(x_i, 2)
    p = min(1, (D_target(y)/D_target(x_i)))
    if (runif(1, min=0, max=1) < p){
      x_i = y
    }
    vecX[i+1] = x_i
  }
  return(vecX)
}

#visualisation
x.points = seq(-10,10, length = 10000)
lines(x.points, D_target(x.points), col = "red")











