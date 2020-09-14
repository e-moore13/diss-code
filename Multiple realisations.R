#27th July 2020

library(tidyverse)
library(mvtnorm)

m= 1000000
alpha=0.05
p = 0.5
t = 0.8

mp1_1 = round((m*t))
mp1_2 = round(((1-t)*m)/9)
mp2 = round(((1-t)*m))

pow = 2.16254
xi1=pow
xi2=pow*sqrt(1-t)

#----------------------------------------------------------------------------------------
probs <- function(mp1_1, mp1_2, mp2, t, xi1, xi2,p,alpha){
  Z1_1 <- rnorm(mp1_1, xi1*sqrt(t), 1)
  Z1_2 <- rnorm(mp1_2, xi1*sqrt(1-t), 1)
  Z1 <- sqrt(t)*Z1_1 + sqrt(1-t)*Z1_2
  Z2 <- rnorm(mp2, (xi2+(Z1_2 - xi1*sqrt(1-t))/2), sqrt(1-t))
  
  AZ1_1 <- 1 - pnorm((qnorm(1-alpha)-((sqrt(t))*Z1_1))/(sqrt(1-t)))
  
  I1 <- Z1 > qnorm(1-alpha)
  I2 <- Z2 > qnorm(1-alpha)
  I1_2 <- Z1_2 > qnorm(1-AZ1_1)
  
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
  
  I12 <- p12 < AZ1_1 
  
  L1 <- mean(I1)
  L2 <- mean(I2)
  Lboth <- mean(I12)
  
  GI1 <- I1*I12
  G1 <- mean(GI1)
  GI2 <- I2*I12
  G2 <- mean(GI2)
  GI12 <- GI1*GI2
  G12 <- mean(GI12)
  
  GI1Only <- I1*I12*(1-I2)
  G1O <- mean(GI1Only)
  GI2Only <- I2*I12*(1-I1)
  G2O <- mean(GI2Only)
  
 #msg1 <- paste("The local hypotheses for rejecting H01, H02 and H0,12 are",L1, L2, "and", Lboth)
 msg2 <- paste("The global hypotheses for rejecting H01, H02 and for both are",G1, G2, "and", G12)
 #msg3 <- paste("The global hypotheses for rejecting H01, H02 only are", G1O, "and", G2O)
  
  #print(msg1)
  print(msg2)
  #print(msg3)
}


for (i in c(2,3,4,5,6,7,8,9,10)){
probs(mp1_1=round((m*t)/i), mp1_2, mp2, t= 0.1, xi1, xi2,p,alpha)
}

mp1_2 = round(((1-t)*m)/2)
for (i in c(2,3,4,5,6,7,8,9,10)){
  probs(mp1_1=round((m*t)/i), mp1_2, mp2, t= 0.5, xi1, xi2,p,alpha)
}
 