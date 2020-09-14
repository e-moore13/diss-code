# Bayes Optimisation

#-Load Packages-------------------------------------
library(tidyverse)
library(mvtnorm)

#--Variables------------------------------------------
m= 1000
alpha=0.05
p = 0.5
t = 0.5

mp1_1 = round((m*t))
mp1_2 = round(((1-t)*m)/5)
mp2 = round(((1-t)*m))

#-priors-treatment effects--------------------------
n1 =100
n2 = 100
sd1 = 1.5
sd2 = 1.5

theta1 <- rnorm(n1, 1.446169333, sd1)
theta2 <- rnorm(n2, 1.446169333, sd2)

#-Calculations---------------------------------------

bayes <- function(th1, th2, mp1_1, mp1_2, mp2, t, xi1, xi2,p,alpha){
xi1 <- th1*sd1
xi2 <- th2*sd2*sqrt(1-t)

#Test Statistics
Z1_1 <- rnorm(mp1_1, xi1*sqrt(t), 1)
Z1_2 <- rnorm(mp1_2, xi1*sqrt(1-t), 1)
Z1 <- sqrt(t)*Z1_1 + sqrt(1-t)*Z1_2
Z2 <- rnorm(mp2, (xi2+(Z1_2 - xi1*sqrt(1-t))/2), sqrt(1-t))

#Conditional Error
AZ1_1 <- 1 - pnorm((qnorm(1-alpha)-((sqrt(t))*Z1_1))/(sqrt(1-t)))

#Dunnett's Test
Z12 <- Z1_2*(Z1_2>=Z2) + Z2*(Z1_2<Z2)

nd <- 2
corr <- diag(nd)
corr[lower.tri(corr)] <- p 
corr[upper.tri(corr)] <- p 

p12 <- rep(-1, length(Z12))
for(i in 1:length(Z12))
{
  p12[i] <- 1-pmvnorm(lower=c(-Inf,-Inf), upper=c(Z12[i],Z12[i]), mean=c(0,0), corr=corr)
}

#Indicators of Events 
R1 <- Z1 > qnorm(1-alpha)
R2 <- Z2 > qnorm(1-alpha)
R1_2 <- Z1_2 > qnorm(1-AZ1_1)
R12 <- p12 < AZ1_1 

#Global Rejections
G1 <- mean(R1*R12)
G2 <- mean(R2*R12)
GBoth <- mean(G1*G2)

gain <- mean(th1*R1 + th2*R2)
maxgain <- max(gain)

#msg1 <- paste("The expected gain is", gain)
print(gain)
}

# Loop 
for (j in 1:length(theta1)){
  a[j] <- bayes(theta1[j], theta2[j], mp1_1, mp1_2, mp2, t, xi1, xi2, p ,alpha)
}

max(a) # maximum expected gain 
 
