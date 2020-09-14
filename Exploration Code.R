#24th July 2020

#--Load Packages--------------------------------------------------------
library(tidyverse)
library(mvtnorm)

#--Exploration----------------------------------------------------------

#--Variables------------------------------------------------------------
m= 1000000
alpha=0.05
p = 0.5

t = 0.1
mp1_1 = round((m*t))
mp1_2 = round(((1-t)*m)/6)
mp2 = round(((1-t)*m))

xi1=pow
xi2=pow*sqrt(1-t)

pow = 2.16254   #Gives power  of 0.7

#--Simulate Data--------------------------------------------------------
Z1_1 <- rnorm(mp1_1, xi1*sqrt(t), 1)
Z1_2 <- rnorm(mp1_2, xi1*sqrt(1-t), 1)
Z1 <- sqrt(t)*Z1_1 + sqrt(1-t)*Z1_2
Z2 <- rnorm(mp2, (xi2+(Z1_2 - xi1*sqrt(1-t))/2), sqrt(1-t))

#--Conditional error principle------------------------------------------
AZ1_1 <- 1 - pnorm((qnorm(1-alpha)-((sqrt(t))*Z1_1))/(sqrt(1-t)))

#--Indicator variables for local tests of H01 and H02-------------------
I1 <- Z1 > qnorm(1-alpha)
I2 <- Z2 > qnorm(1-alpha)
I1_2 <- Z1_2 > qnorm(1-AZ1_1)

#--Intersection Hypothesis----------------------------------------------
Z12 <- Z1_2*(Z1_2>=Z2) + Z2*(Z1_2<Z2)

n <- 2
corr <- diag(n)
corr[lower.tri(corr)] <- p 
corr[upper.tri(corr)] <- p 

p12 <- rep(-1, length(Z12))
for(i in 1:length(Z12))
{
  p12[i] <- 1-pmvnorm(lower=c(-Inf,-Inf), upper=c(Z12[i],Z12[i]), mean=c(0,0), corr=corr)
}

#--Indicator variable of H0,12---------------------------
I12 <- p12 < AZ1_1 

#--Local level test-------------------------------
mean(I1)
mean(I2)
mean(I12)

#--Global level probabilities of rejecting------------------------------
GI1 <- I1*I12
mean(GI1)
GI2 <- I2*I12
mean(GI2)
GI12 <- GI1*GI2
mean(GI12)

GI1Only <- I1*I12*(1-I2)
mean(GI1Only)
GI2Only <- I2*I12*(1-I1)
mean(GI2Only)
