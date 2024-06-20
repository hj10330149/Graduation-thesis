# Bivariate skew normal and bivariate gamma distribution

#### BSN mean var corr and Ratio of BSN mean var skewness ####
#Using the "rmsn" function from the "sn" package, various parameters were employed to generate 
#Bivariate Skew Normal (BSN) data and the ratio of Bivariate Skew Normal (RBSN) data.
#Also calculate the mean, variance, correlation, and skewness of these data.

library(sn)
library(moments)
Mean_BSN=c()
Var_BSN=c()
Corr_BSN=c()
Mean_R_BSN=c()
Var_R_BSN=c()
p0_R_BSN=c()
skew_R_BSN=c()
for (k in c(0.1,1,2)) {
  for (rho in c(0.65,0.80,0.95)) {
    xi=c(10,10)
    alpha=c(5,5)
    w1=w2=k
    Omega=matrix(c(w1^2,w1*w2*rho,w1*w2*rho,w2^2),ncol=2)
    data=rmsn(10000000, xi, Omega, alpha, tau = 0, dp = NULL)
    Mean_BSN=c(Mean_BSN, mean(data[,1]))
    Var_BSN=c(Var_BSN, var(data[,1]))
    Corr_BSN=c(Corr_BSN, cor(data[,1],data[,2]))
    
    R=data[,1]/data[,2]
    Mean_R_BSN=c(Mean_R_BSN, mean(R))
    Var_R_BSN=c(Var_R_BSN, var(R))
    
    Y=(R[1:5000000]-R[5000001:10000000])^2/2
    p0_R_BSN=c(p0_R_BSN, round(sum(Y>var(R))/length(Y),4) )
    skew_R_BSN=c(skew_R_BSN, skewness(R))
    
  }
}
xi=rep(10,9)
alpha=rep(5,9)
w=c(rep(0.1,3), rep(1,3), rep(2,3))
rho=c(0.65,0.80,0.95, 0.65,0.80,0.95, 0.65,0.80,0.95)
df1=data.frame(xi,alpha,w,rho,Mean_BSN,Var_BSN,Corr_BSN,Mean_R_BSN,Var_R_BSN,skew_R_BSN)
#write.csv(df1,file="C:/Users/WEIYU/Desktop/政治大學/品管研究計畫/Meeting/bivariate gamma/excel_BG/IC_table_BSN.csv")



#### BG mean var corr and Ratio of BG mean var skewness ####
#Need to run the BSN mean var corr code above
#Use the BSN mean var corr to find the BG parameter which have the same mean var corr

library(moments)
generate_UV=function(total,a,b,c,d){ 
  X=rgamma(total,a,scale=d) #rate=lambda, scale=beta        
  Y=rgamma(total,b,scale=d) 
  Z=rgamma(total,c,scale=d)
  U=(X+Y)
  V=(X+Z)
  UV=matrix(c(U,V),ncol=2,byrow=F)
  return(UV)
}
par=function(mean,var,cor){
  beta=var/mean
  a=cor*mean/beta
  b=mean/beta-a
  c=mean/beta-a
  return(c(a,b,c,beta))
}
a=c()
b=c()
c=c()
d=c()
for (i in 1:length(Mean_BSN)) {
  parameters=par(Mean_BSN[i],Var_BSN[i],Corr_BSN[i])
  a=c(a,parameters[1])
  b=c(b,parameters[2])
  c=c(c,parameters[3])
  d=c(d,parameters[4]) #d=beta
}

Mean_BG=c() #The mean var corr of Bivariate Gamma
Var_BG=c()
Corr_BG=c()
Mean_R_BG=c()
Var_R_BG=c()
p0_R_BG=c()
skew_R_BG=c()
for (j in 1:length(a)) {
  data=generate_UV(total=10000000,a[j],b[j],c[j],d[j])
  Mean_BG=c(Mean_BG, mean(data[,1]))
  Var_BG=c(Var_BG, var(data[,1]))
  Corr_BG=c(Corr_BG, cor(data[,1],data[,2]))
  
  R=data[,1]/data[,2]
  Mean_R_BG=c(Mean_R_BG, mean(R))
  Var_R_BG=c(Var_R_BG, var(R))
  
  Y=(R[1:5000000]-R[5000001:10000000])^2/2
  p0_R_BG=c(p0_R_BG, round(sum(Y>var(R))/length(Y),4) )
  skew_R_BG=c(skew_R_BG, skewness(R))
}

df2=data.frame(a,b,c,d,Mean_BG,Var_BG,Corr_BG,Mean_R_BG,Var_R_BG,skew_R_BG)
#write.csv(df2,file="C:/Users/WEIYU/Desktop/政治大學/品管研究計畫/Meeting/bivariate gamma/excel_BG/IC_BG_W3.csv")






#### Out of control bivariate skew normal ####
#Find the out of control omega(k) under delta=0.5, 1.5, 2.0, 2.5, 3.0

library(sn)
Find_oc_parameter_BSN=function(i,j,k,rho,var0,delta){
  xi=c(i,i)
  alpha=c(j,j)
  w1=w2=k
  Omega=matrix(c(w1^2,w1*w2*rho,w1*w2*rho,w2^2),ncol=2)
  data_oc=rmsn(10000000, xi, Omega, alpha, tau = 0, dp = NULL)
  R_oc=data_oc[,1]/data_oc[,2]
  delta_hat=sd(R_oc)/sqrt(var0) #delta 是標準差相除
  return(delta_hat-delta) #when delta=0 then return delta_hat
}

#Find OC Omega parameter
i=10;j=5;k=1;rho=0.80;var0=0.00351234446412953;delta=0
delta_hat=c()
for (k in c(3.288)) {
  del=Find_oc_parameter_BSN(i,j,k,rho,var0,delta)
  delta_hat=c(delta_hat,del)
}
print(delta_hat)
#Under k=0.1, var0=0.0000394
#OC delta=0.5,    Find  k=0.050, delta_hat=0.5044
#OC delta=1.5,    Find  k=0.150, delta_hat=1.5021
#OC delta=2.0,    Find  k=0.200, delta_hat=1.9956
#OC delta=2.5,    Find  k=0.251, delta_hat=2.4965
#OC delta=3.0,    Find  k=0.303, delta_hat=3.0033
#Under k=1, var0=0.00351234
#OC delta=0.5,    Find  k=0.485, delta_hat=0.500506
#OC delta=1.5,    Find  k=1.543, delta_hat=1.500364 
#OC delta=2.0,    Find  k=2.110, delta_hat=2.002631
#OC delta=2.5,    Find  k=2.692, delta_hat=2.499536
#OC delta=3.0,    Find  k=3.288, delta_hat=3.000275
#Under k=2, var0=0.01275085
#OC delta=0.5,    Find  k=0.950, delta_hat=0.5001
#OC delta=1.5,    Find  k=3.120, delta_hat=1.5008
#OC delta=2.0,    Find  k=4.265, delta_hat=2.0001
#OC delta=2.5,    Find  k=5.385, delta_hat=2.5006
#OC delta=3.0,    Find  k=6.300, delta_hat=2.9763



### Find BSN parameters
#w1=w2=c(0.1, 1, 2)
#var0=c(0.0000394, 0.00351234, 0.01275085)
#p0=c(0.317, 0.316, 0.310)
#w1=0.1 k=c(0.05, 0.15, 0.20, 0.251, 0.303)
#w1=1   k=c(0.485, 1.543, 2.110, 2.692, 3.288)
#w1=2   k=c(0.95, 3.120, 4.265, 5.385, 6.300)
var0=0.00351234
p0=0.316

library(moments)
delta_hat=c()
Mean=c()
Var=c()
Corr=c()
Mean_R=c()
Var_R=c()
Skew_R=c()
P1=c()
P1P0=c()
for (k in c(0.485, 1.543, 2.110, 2.692, 3.288) ) {
  w1=w2=k
  xi=c(10,10)
  alpha=c(5,5)
  rho=0.80
  Omega=matrix(c(w1^2,w1*w2*rho,w1*w2*rho,w2^2),ncol=2)
  
  data_oc=rmsn(10000000, xi=xi, Omega=Omega, alpha=alpha, tau = 0, dp = NULL)
  R_oc=data_oc[,1]/data_oc[,2]
  Y_oc=(R_oc[1:5000000]-R_oc[5000001:10000000])^2/2
  
  delta_hat=c(delta_hat, sd(R_oc)/sqrt(var0)) #delta 是標準差相除
  Mean=c(Mean, mean(data_oc[,1]))
  Var=c(Var, var(data_oc[,1]))
  Corr=c(Corr, cor(data_oc[,1],data_oc[,2]))
  Mean_R=c(Mean_R, mean(R_oc))
  Var_R=c(Var_R, var(R_oc))
  Skew_R=c(Skew_R, skewness(R_oc))
  p1=round(sum(Y_oc>var0)/length(Y_oc),4)
  P1=c(P1, p1)
  P1P0=c(P1P0, p1/p0)
}
Var_R=round(Var_R,7)
OC_BSN=data.frame(delta_hat,Mean,Var,Corr,Mean_R,Var_R,Skew_R)
#write.csv(OC_BSN,file="C:/Users/WEIYU/Desktop/政治大學/品管研究計畫/Meeting/EWMA Sign chart/OC data/OC_BSN_omega_1.csv")




#### Out of control bivariate gamma ####
#Use the OC BSN mean var corr to find the OC BG parameter which have the same mean var corr

generate_UV=function(total,a,b,c,beta){ 
  X=rgamma(total,a,scale=beta) #rate=lambda, scale=beta        
  Y=rgamma(total,b,scale=beta) 
  Z=rgamma(total,c,scale=beta)
  U=(X+Y)
  V=(X+Z)
  UV=matrix(c(U,V),ncol=2,byrow=F)
  return(UV)
}
par=function(mean,var,cor){
  beta=var/mean
  a=cor*mean/beta
  b=mean/beta-a
  c=mean/beta-a
  return(c(a,b,c,beta))
}

# Find the parameters of BG 
Mean=c(10.03765, 10.11290, 10.15049, 10.18898, 10.22806, 10.36506, 11.16091, 11.58849, 12.02638, 12.47476, 10.71508, 12.34917, 13.21085, 14.05455, 14.74312)
Var=c(0.001082905, 0.009743170, 0.017337344, 0.027320447, 0.039726442, 0.1018997, 1.0316968, 1.9294569, 3.1405259, 4.6857097, 0.3912914,  4.2204975,  7.8840184, 12.5689357, 17.2087858)
Corr=c(0.5381808, 0.5385451, 0.5383685, 0.5383924, 0.5380886, 0.5387229, 0.5386399, 0.5384643, 0.5384145, 0.5386356, 0.5385072, 0.5387657, 0.5386846, 0.5387195, 0.5381384)
a=c()
b=c()
c=c()
d=c()
for (i in 1:length(Mean)) {
  parameters=par(mean=Mean[i],var=Var[i],cor=Corr[i])
  a=c(a,parameters[1])
  b=c(b,parameters[2])
  c=c(c,parameters[3])
  d=c(d,parameters[4]) #d=beta
}

library(moments)
var0=c( rep(0.000039,5), rep(0.0035,5), rep(0.0124,5) )
p0=c( rep(0.318,5), rep(0.316,5), rep(0.312,5) )
delta_hat=c()
Mean_BG=c() #The mean var corr of Bivariate Gamma
Var_BG=c()
Cor=c()
Mean_R=c()
Var_R=c()
Skew_R=c()
P1=c()
for (j in 1:length(Mean)) {
  data_oc=generate_UV(total=10000000,a[j],b[j],c[j],d[j])
  Mean_BG=c(Mean_BG, mean(data_oc[,1]))
  Var_BG=c(Var_BG, var(data_oc[,1]))
  Cor=c(Cor, cor(data_oc[,1], data_oc[,2]))
  
  R_oc=data_oc[,1]/data_oc[,2]
  Y_oc=(R_oc[1:5000000]-R_oc[5000001:10000000])^2/2
  delta_hat=c(delta_hat, sd(R_oc)/sqrt(var0[j])) #delta 是標準差相除
  Mean_R=c(Mean_R, mean(R_oc))
  Var_R=c(Var_R, var(R_oc))
  Skew_R=c(Skew_R, skewness(R_oc))
  p1=round(sum(Y_oc>var0[j])/length(Y_oc),4)
  P1=c(P1, p1)
}
Var_R=round(Var_R,6)
OC_BG=data.frame(delta_hat,a,b,c,d,Mean,Var,Cor,Mean_R,Var_R,Skew_R)

