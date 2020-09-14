# 8th August 2020

#-Load Packages-------------------------------------
library(tidyverse)
library(mvtnorm)
library(reshape)
library(RColorBrewer)
library(wesanderson)

#-variables-----------------------------------------
sd = 1
m = 1000
alpha = 0.05
t = 0.5
#stage 1 randomisation proportions
r0_1 = 1 
r1_1 = 1
#stage 2 randomisation proportions
r0_2 = 1 
r1_2 = 1/10
r2_2 = 1
#Sample sizes for stage 1 
s0_1 = round((r0_1*m*t))
s1_1 = round((r1_1*m*t))
#Sample sizes for stage 2
s0_2 = round(r0_2*(1-t)*m)
s1_2 = round(r1_2*(1-t)*m)
s2_2 = round(r2_2*(1-t)*m)
ybar = 2
xbar1 = 4
xbar2 = 6

#-Treatments----------------------------------------
#-Control-------------------------------------------

#Control treatment simulation for stage 1
Y = rnorm(s0_1, ybar, sd)    
Ybar = mean(Y) #Gives mean of the simulated Ys

#Control treatment simulation for stage 2
Y2 = rnorm(s0_2, ybar, sd) 
Ybar2 = mean(Y2)

#-Treatment 1----------------------------------------

#Treatment 1 simulation for stage 1
X1 = rnorm(s1_1, xbar1, sd)    
X1bar1 = mean(X1) #Gives mean of the simulated Xs

#Treatment 1 simulation for stage 2
X12 = rnorm(s1_2, xbar1, sd) 
X1bar2 = mean(X12)

#-Treatment 2----------------------------------------

#Treatment 2 
X2 <- rnorm(s2_2, xbar2, sd) # randomise new treatment 
Xbar2 = mean(X2)

#-Test Statistics------------------------------------

#Test statistic for treatment 1 vs control
Z11 <- (X1bar1 - Ybar)/(sqrt(4*sd)) 
Z12 <- (X1bar2 - Ybar2)/(sqrt(4*sd))
Z1 = sqrt(t)*Z11 + sqrt(1-t)*Z12

#Test statistic for treatment 2 vs control
Z2 <- (Xbar2 - Ybar2)/(sqrt(4*sd)) 

#-Dunnett Test --------------------------------------

ZD <- Z12*(Z12>=Z2) + Z2*(Z12<Z2)
ZD 

#Distribution
p = 0.5 # correlation
mean = rep(0, 2)
corr <- diag(2)
corr[lower.tri(corr)] <- p
corr[upper.tri(corr)] <- p
lower = -Inf
upper = ZD

# p value for intersection
pd <- 1 - pmvnorm(lower, upper, mean, corr)

#- Conditional Error Function ------------------------

AZ11 <- 1 - pnorm((qnorm(1-alpha)-((sqrt(t))*Z11))/(sqrt(1-t)))  

#- Closed Testing Procedure --------------------------

#Reject H01 at local alpha level when TRUE:
Z12 > qnorm(1-AZ11)
p12 <- 1 - pnorm(Z12)

#Reject H02 at local alpha level when TRUE: 
Z2>qnorm(1-alpha)
p2 <- 1 - pnorm(Z2)

#Reject H012 when TRUE:
pd < AZ11

p12
p2
pd
#------------------------------------------------------

tau <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
r1_2 <- c(1/2,1/3,1/4,1/5,1/6,1/7,1/8,1/9,1/10)

rpvals1 <- c(0.143045, 0.1576106, 0.1365907, 0.1529637, 0.1561877, 0.1654823, 0.1430024, 0.1510821,  0.1898487)
rpvals2 <- c(0.02333906, 0.0227003, 0.02030876, 0.025696, 0.02197633, 0.02177183, 0.0216713, 0.02371666, 0.02582957)
rpvals12 <- c(0.04247825, 0.04136002, 0.03715723, 0.04658928, 0.04009046, 0.03973143, 0.03955487, 0.04313844, 0.04682156)

tpvals1 <- c(0.1504392, 0.1702923, 0.1579101, 0.1412704, 0.1635834, 0.1423823, 0.1882808, 0.1484693, 0.1479845)
tpvals2 <- c(0.02314694, 0.02230831, 0.02218743, 0.02233311, 0.02113309, 0.02225047, 0.02350154, 0.01934923, 0.02248004)
tpvals12 <- c(0.0421421, 0.04067292, 0.04046089, 0.0407164, 0.0386088, 0.04057146, 0.0427624, 0.03546357, 0.04097401)

# Plot - Vary r1_2
df1 <- data.frame(A=r1_2, H01=rpvals1, H02=rpvals2, H012=rpvals12)
d1 <- melt(df1, id.vars="A")
p1 <- ggplot(d1, aes(A, value, col=variable)) + geom_point(size=3) + geom_hline(yintercept=0.05, col="red", size=1.1)
p11 <- p1 + labs(color="p-values", x="r'_1", y="p-values")       
p11 + scale_color_manual(values=brewer.pal(n=3, name="Accent"))  

# Plot - Vary tau
df2 <- data.frame(A=tau, H01=tpvals1, H02=tpvals2, H012=tpvals12)
d2 <- melt(df2, id.vars="A")
p2 <- ggplot(d2, aes(A, value, col=variable)) + geom_point(size=3) + geom_hline(yintercept=0.05, col="red", size=1.1)
p22 <- p2 + labs(color="p-values", x="tau", y="p-values")       
p22 + scale_color_manual(values=brewer.pal(n=3, name="Accent"))  

