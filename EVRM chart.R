# Mood based EWMA variance of ratio control chart (EVRM control chart)

#### Check the simulated Mean Var are same as the real Mean Var ####

generate_R_gamma=function(total,a,b,c,beta){ 
  X=rgamma(total,a,scale=beta)  #rate=lambda, scale=beta     
  Y=rgamma(total,b,scale=beta) 
  Z=rgamma(total,c,scale=beta)
  R=(X+Y)/(X+Z)
  return(R)
}
N=50;n=4
a_ic=a_oc=41.12; 
b_ic=b_oc=35.22; 
c_ic=c_oc=35.22; 
beta_ic=beta_oc=0.15

start=Sys.time()
N_pool=N+n
Mean=n*(N_pool^2-1)/12
Var=N*n*(N_pool+1)*(N_pool^2-4)/180
M=c()
for (num in 1:100000) {
  R1=generate_R_gamma(N,a_ic,b_ic,c_ic,beta_ic)
  R2=generate_R_gamma(n,a_oc,b_oc,c_oc,beta_oc)
  R_pool=sort(c(R1,R2))
  b1=c()
  for (j in 1:N_pool) {  
    if (R_pool[j] %in% R2) {
      b1=c(b1,(j-(N_pool+1)/2)^2) 
    }
  }
  M=c(M, sum(b1))
}
mean(M);var(M)
Mean;Var
end= Sys.time()
end-start

#### EVRM chart function ####
#The EVRM function.
#When set=0 it output ARL0_hat-ARL*2. This is used for "uniroot" to find L1 and L2
#When set!=0 it output ARL0_hat

library(sn)
Find_L1_Mood_EWMA_BSN=function(N,n,i,j,k_ic,k_oc,rho,L1,lambda,ARL0,set){
  #start= Sys.time()
  M=c()
  EWMA=c()
  Z=c()
  RL=c()
  
  N_pool=N+n
  mean=n*(N_pool^2-1)/12
  var=N*n*(N_pool+1)*(N_pool^2-4)/180
  
  xi=c(i,i)
  alpha=c(j,j)
  w1_ic=w2_ic=k_ic
  Omega_ic=matrix(c(w1_ic^2,w1_ic*w2_ic*rho,w1_ic*w2_ic*rho,w2_ic^2),ncol=2)
  w1_oc=w2_oc=k_oc
  Omega_oc=matrix(c(w1_oc^2,w1_oc*w2_oc*rho,w1_oc*w2_oc*rho,w2_oc^2),ncol=2)
  
  for (num in 1:10000) {
    #R1=R_ic (reference data)
    data_ic=rmsn(N, xi, Omega_ic, alpha, tau = 0, dp = NULL)
    R1=data_ic[,1]/data_ic[,2]
    
    #R2=R_oc (new data)
    data_oc=rmsn(n, xi, Omega_oc, alpha, tau = 0, dp = NULL)
    R2=data_oc[,1]/data_oc[,2]
    
    R_pool=sort(c(R1,R2))
    b1=c()
    for (count in 1:N_pool) {  
      if (R_pool[count] %in% R2) {
        b1=c(b1,(count-(N_pool+1)/2)^2) 
      }
    }
    i=1
    M[i]=sum(b1)
    EWMA[i]=lambda*M[i]+(1-lambda)*mean 
    Z[i]=(EWMA[i]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^(2*i)))
    while (Z[i]<L1) { 
      i=i+1
      
      #R1=R_ic (reference data)
      data_ic=rmsn(N, xi, Omega_ic, alpha, tau = 0, dp = NULL)
      R1=data_ic[,1]/data_ic[,2]
      
      #R2=R_oc (new data)
      data_oc=rmsn(n, xi, Omega_oc, alpha, tau = 0, dp = NULL)
      R2=data_oc[,1]/data_oc[,2]
      
      R_pool=sort(c(R1,R2))
      b1=c()
      for (count in 1:N_pool) {  
        if (R_pool[count] %in% R2) {
          b1=c(b1,(count-(N_pool+1)/2)^2) 
        }
      }
      M[i]=sum(b1)
      EWMA[i]=lambda*M[i]+(1-lambda)*EWMA[i-1]
      Z[i]=(EWMA[i]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^(2*i)))
      if (i>=5000) {break}
    }
    RL=c(RL,i)
  }
  # sum(RL==5000);mean(RL[-which(RL==5000)]);mean(RL)
  # end= Sys.time()
  # end-start 
  
  count=sum(RL==5000)
  if (count>0) {
    ARL0_hat=mean(RL[-which(RL==5000)])
  }else{
    ARL0_hat=mean(RL)
  }
  if (set==0) {
    return(ARL0_hat-2*ARL0)
  }else{
    return(ARL0_hat)
  }
}
Find_L2_Mood_EWMA_BSN=function(N,n,i,j,k_ic,k_oc,rho,L1,L2,lambda,ARL0,set=0){
  #start= Sys.time()
  M=c()
  EWMA=c()
  Z=c()
  RL=c()
  
  N_pool=N+n
  mean=n*(N_pool^2-1)/12
  var=N*n*(N_pool+1)*(N_pool^2-4)/180
  
  xi=c(i,i)
  alpha=c(j,j)
  w1_ic=w2_ic=k_ic
  Omega_ic=matrix(c(w1_ic^2,w1_ic*w2_ic*rho,w1_ic*w2_ic*rho,w2_ic^2),ncol=2)
  w1_oc=w2_oc=k_oc
  Omega_oc=matrix(c(w1_oc^2,w1_oc*w2_oc*rho,w1_oc*w2_oc*rho,w2_oc^2),ncol=2)
  
  for (num in 1:10000) {
    #R1=R_ic (reference data)
    data_ic=rmsn(N, xi, Omega_ic, alpha, tau = 0, dp = NULL)
    R1=data_ic[,1]/data_ic[,2]
    
    #R2=R_oc (new data)
    data_oc=rmsn(n, xi, Omega_oc, alpha, tau = 0, dp = NULL)
    R2=data_oc[,1]/data_oc[,2]
    
    R_pool=sort(c(R1,R2))
    b1=c()
    for (count in 1:N_pool) {  
      if (R_pool[count] %in% R2) {
        b1=c(b1,(count-(N_pool+1)/2)^2) 
      }
    }
    i=1
    M[i]=sum(b1)
    EWMA[i]=lambda*M[i]+(1-lambda)*mean 
    Z[i]=(EWMA[i]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^(2*i)))
    while (Z[i]<L1 && Z[i]>-L2) { 
      i=i+1
      
      #R1=R_ic (reference data)
      data_ic=rmsn(N, xi, Omega_ic, alpha, tau = 0, dp = NULL)
      R1=data_ic[,1]/data_ic[,2]
      
      #R2=R_oc (new data)
      data_oc=rmsn(n, xi, Omega_oc, alpha, tau = 0, dp = NULL)
      R2=data_oc[,1]/data_oc[,2]
      
      R_pool=sort(c(R1,R2))
      b1=c()
      for (count in 1:N_pool) {  
        if (R_pool[count] %in% R2) {
          b1=c(b1,(count-(N_pool+1)/2)^2) 
        }
      }
      M[i]=sum(b1)
      EWMA[i]=lambda*M[i]+(1-lambda)*EWMA[i-1]
      Z[i]=(EWMA[i]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^(2*i)))
      if (i>=5000) {break}
    }
    RL=c(RL,i)
  }
  # sum(RL==5000);mean(RL[-which(RL==5000)]);mean(RL)
  # end= Sys.time()
  # end-start 
  
  count=sum(RL==5000)
  if (count>0) {
    ARL0_hat=mean(RL[-which(RL==5000)])
    MRL=median(RL[-which(RL==5000)])
    SDRL=sd(RL[-which(RL==5000)])
  }else{
    ARL0_hat=mean(RL)
    MRL=median(RL)
    SDRL=sd(RL)
  }
  if (set==0) {
    return(ARL0_hat-ARL0)
  }else{
    return(c(ARL0_hat,MRL,SDRL))
  }
}


#### EVRM chart L1、L2、ARL0 ####
#Use the EVRM function and "uniroot" function to find the L1, L2 and estimate ARL0

start=Sys.time()
i=10;j=5;k_ic=2;k_oc=2;rho=0.80
N=50
Ans=c()
for (n in c(4,8,16)) {
  ARL0_hat=0 
  while (abs(ARL0_hat-740.8)>1) {  
    L1=uniroot(Find_L1_Mood_EWMA_BSN,interval=c(2.2, 2.6),N=N,n=n,i=i,j=j,k_ic=k_ic,k_oc=k_oc,
               rho=rho,lambda=0.05,ARL0=370.4,set=0,extendInt = "yes")$root
    ARL0_hat=Find_L1_Mood_EWMA_BSN(N,n,i,j,k_ic,k_oc,rho,L1,lambda=0.05,ARL0=370.4,set=1)
  }
  while (abs(ARL0_hat-370.4)>1) {  
    L2=uniroot(Find_L2_Mood_EWMA_BSN,interval=c(2.2, 2.6),N=N,n=n,i=i,j=j,k_ic=k_ic,
               k_oc=k_oc,rho=rho,L1=L1,lambda=0.05,ARL=370.4,set=0,extendInt = "yes")$root
    RL=Find_L2_Mood_EWMA_BSN(N,n,i,j,k_ic,k_oc,rho,L1,L2,lambda=0.05,ARL0=370.4,set=1) 
    ARL0_hat=RL[1]
  }
  Ans=c(Ans,round(c(L1,L2),4),round(c(RL),2))
}
out=matrix(Ans,ncol=5,byrow=T)
end=Sys.time()
end-start


#### EVRM chart ARL1 ####
#Used the L1, L2 we found above to calculate ARL1

n=c(rep(4,15), rep(8,15), rep(16,15))
i=rep(10,45)
j=rep(5,45)
k_ic=c(rep(0.1,5), rep(1,5), rep(2,5), rep(0.1,5), rep(1,5), rep(2,5), rep(0.1,5), rep(1,5), rep(2,5))
k_oc=c(0.05,0.15,0.2,0.251,0.303, 0.485,1.543,2.110,2.692,3.288, 0.950,3.12,4.265,5.385,6.3, 
       0.05,0.15,0.2,0.251,0.303, 0.485,1.543,2.110,2.692,3.288, 0.950,3.12,4.265,5.385,6.3, 
       0.05,0.15,0.2,0.251,0.303, 0.485,1.543,2.110,2.692,3.288, 0.950,3.12,4.265,5.385,6.3)
rho=rep(0.8,45)

### ARL0=370.4
# N=50 
# L1=c(rep(2.547,15), rep(2.529,15), rep(2.508,15))
# L2=c(rep(2.484,15), rep(2.506,15), rep(2.542,15))

# N=100 
# L1=c(rep(2.548,15), rep(2.531,15), rep(2.516,15))
# L2=c(rep(2.490,15), rep(2.500,15), rep(2.518,15))

# N=500 
# L1=c(rep(2.551,15), rep(2.540,15), rep(2.529,15))
# L2=c(rep(2.470,15), rep(2.498,15), rep(2.515,15))


### ARL0=200
# N=50 
# L1=c(rep(2.255,15), rep(2.238,15), rep(2.225,15))
# L2=c(rep(2.296,15), rep(2.343,15), rep(2.346,15))

# N=100 
# L1=c(rep(2.259,15), rep(2.243,15), rep(2.236,15))
# L2=c(rep(2.273,15), rep(2.298,15), rep(2.310,15))

# N=500 
# L1=c(rep(2.270,15), rep(2.249,15), rep(2.238,15))
# L2=c(rep(2.268,15), rep(2.309,15), rep(2.319,15))

#df=data.frame(n,i,j,k_ic,k_oc,rho,L1,L2)

N=500
L1=c(rep(2.270,15), rep(2.249,15), rep(2.238,15))
L2=c(rep(2.268,15), rep(2.309,15), rep(2.319,15))

#ARL1
start= Sys.time() 
arl1=c()
for (w in 1:length(n)) {
  x=Find_L2_Mood_EWMA_BSN(N,n=n[w],i=i[w],j=j[w],k_ic=k_ic[w],k_oc=k_oc[w]
                          ,rho=rho[w],L1=L1[w],L2=L2[w],lambda=0.05,ARL0=370.4,set=1)
  arl1=c(arl1,x[1])
}
ARL1=matrix(round(arl1,2),ncol=3,byrow = F)
end= Sys.time()
end-start
ARL1
#write.csv(ARL1,file="C:/Users/WEIYU/Desktop/政治大學/品管研究計畫/Meeting/EWMA Mood chart/ARL1/ARL1_200_BSN_N_500_w1.csv")



