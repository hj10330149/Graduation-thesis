# Sign based EWMA variance of ratio control chart (EVRS control chart)

#### EVRS chart function ####
#The EVRS function.
#When set=0 it output ARL0_hat-ARL*2. This is used for "uniroot" to find L1 and L2
#When set!=0 it output ARL0_hat

EWMA_Sign_L1=function(n,L1,lambda,p0,p1,ARL,set){
  EWMAS=c()
  Z=c()
  RL=c()
  mean=0.5*n*p0
  var=0.5*n*p0*(1-p0)
  for (j in 1:10000) {
    S=rbinom(5000,0.5*n,p1) 
    EWMAS[1]=lambda*S[1]+(1-lambda)*mean
    Z[1]=(EWMAS[1]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^2))
    i=1
    while (Z[i]< L1) {
      i=i+1
      EWMAS[i]=lambda*S[i]+(1-lambda)*EWMAS[i-1]
      Z[i]=(EWMAS[i]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^(2*i)))
      if (i>=5000) {break}
    }
    RL=c(RL,i)
  }
  count=sum(RL==5000)
  if (count>0) {
    ARL0_hat=mean(RL[-which(RL==5000)])
  }else{
    ARL0_hat=mean(RL)
  }
  if (set==0) {
    return(ARL0_hat-ARL*2)
  }else{
    return(ARL0_hat)
  }
}
EWMA_Sign_L2=function(n,L1,L2,lambda,p0,p1,ARL,set){
  EWMAS=c()
  Z=c()
  RL=c()
  mean=0.5*n*p0
  var=0.5*n*p0*(1-p0)
  for (j in 1:10000) {
    S=rbinom(5000,0.5*n,p1) 
    EWMAS[1]=lambda*S[1]+(1-lambda)*mean
    Z[1]=(EWMAS[1]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^2))
    i=1
    while (Z[i]<L1 && Z[i]>-L2) {
      i=i+1
      EWMAS[i]=lambda*S[i]+(1-lambda)*EWMAS[i-1]
      Z[i]=(EWMAS[i]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^(2*i)))
      if (i>=5000) {break}
    }
    RL=c(RL,i)
  }
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
    return(ARL0_hat-ARL)
  }else{
    return(c(ARL0_hat,MRL,SDRL))
  }
}


#### EVRS chart L1 L2 and estimate ARL0 under ARL0=370.4 ####
#Change the difference p0 to find the estimate ARL0 under ARL0=370.4

# Bivariate skew normal p0=c(0.317, 0.316, 0.310)
# Bivariate gamma p0=c(0.318, 0.312, 0.291)
start=Sys.time()
for (p in c(0.317, 0.316, 0.310) ) {
  Ans=c()
  for (n in c(4,8,16,32)) {
    ARL0_hat=0 
    while (abs(ARL0_hat-740.8)>1) {  
      L1=uniroot(EWMA_Sign_L1, c(2.5, 2.8),n=n,lambda=0.05,p0=p,p1=p,ARL=370.4,set=0,extendInt = "yes")$root
      ARL0_hat=EWMA_Sign_L1(n=n,L1=L1,lambda=0.05,p0=p,p1=p,set=1)
    }
    while (abs(ARL0_hat-370.4)>1) {  
      L2=uniroot(EWMA_Sign_L2, c(2.2, 2.5),n=n,L1=L1,lambda=0.05,p0=p,p1=p,ARL=370.4,set=0,extendInt = "yes")$root
      RL=EWMA_Sign_L2(n=n,L1=L1,L2=L2,lambda=0.05,p0=p,p1=p,set=1) 
      ARL0_hat=RL[1]
    }
    Ans=c(Ans,round(c(L1,L2),4),round(c(RL),2))
  }
  out=matrix(Ans,ncol=5,byrow=T)
  path=paste0('C:/Users/WEIYU/Desktop/政治大學/品管研究計畫/Meeting/EWMA Sign chart/ARL0/ARL370_BG_p', p,'.csv')
  write.csv(out,file=path)
}
end=Sys.time()
end-start


#### EVRS chart L1 L2 and estimate ARL0 under ARL0=200 ####
#Change the difference p0 to find the estimate ARL0 under ARL0=200

# Bivariate skew normal p0=c(0.317, 0.310, 0.293)
# Bivariate gamma p0=c(0.318, 0.312, 0.291)
# Change the excel name when the distribution change
start=Sys.time()
for (p in c(0.316)) {
  Ans=c()
  for (n in c(4,8,16,32)) {
    ARL0_hat=0 
    while (abs(ARL0_hat-400)>1) {  
      L1=uniroot(EWMA_Sign_L1, c(2.1, 2.4),n=n,lambda=0.05,p0=p,p1=p,ARL=200,set=0,extendInt = "yes")$root
      ARL0_hat=EWMA_Sign_L1(n=n,L1=L1,lambda=0.05,p0=p,p1=p,set=1)
    }
    while (abs(ARL0_hat-200)>1) {  
      L2=uniroot(EWMA_Sign_L2, c(2.1, 2.4),n=n,L1=L1,lambda=0.05,p0=p,p1=p,ARL=200,set=0,extendInt = "yes")$root
      RL=EWMA_Sign_L2(n=n,L1=L1,L2=L2,lambda=0.05,p0=p,p1=p,set=1) 
      ARL0_hat=RL[1]
    }
    Ans=c(Ans,round(c(L1,L2),4),round(c(RL),2))
  }
  out=matrix(Ans,ncol=5,byrow=T)
  path=paste0('C:/Users/WEIYU/Desktop/政治大學/品管研究計畫/Meeting/EWMA Sign chart/ARL0/ARL200_BG_p', p,'.csv')
  write.csv(out,file=path)
}
end=Sys.time()
end-start



#### EVRS chart ARL1 under ARL0=370.4 ####
#Use different L1, L2, p0 and p1 to calculate the ARL1 under ARL0=370.4

### Bivariate Skew Normal ARL0=370.4
# p0=0.317
# L1=c(2.5861, 2.5775, 2.5431, 2.5340)
# L2=c(2.4175, 2.4458, 2.4835, 2.4998)
# p1=c(0.0465, 0.5032, 0.6145, 0.6873, 0.7376)

# p0=0.316
# L1=c(2.5851, 2.5712, 2.5396, 2.5307)
# L2=c(2.4201, 2.4553, 2.4881, 2.5140)
# p1=c(0.0455, 0.5002, 0.6098, 0.6787, 0.7258)

# p0=0.310
# L1=c(2.5934, 2.5819, 2.5449, 2.5336)
# L2=c(2.4124, 2.4462, 2.4786, 2.5011)
# p1=c(0.0458, 0.4879, 0.5903, 0.6518, 0.6874)


### Bivariate gamma ARL0=370.4
# p0=0.318
# L1=c(2.5877, 2.5769, 2.5374, 2.5307)
# L2=c(2.4186, 2.4563, 2.4929, 2.5053)
# p1=c(0.0476, 0.5053, 0.6159, 0.6886, 0.7389)

# p0=0.316
# L1=c(2.5908, 2.5726, 2.5405, 2.5306)
# L2=c(2.4075, 2.4615, 2.4784, 2.5194)
# p1=c(0.0457, 0.4998, 0.6069, 0.6756, 0.7215)

# p0=0.312
# L1=c(2.5908, 2.5790, 2.5432, 2.5351)
# L2=c(2.4127, 2.4519, 2.4682, 2.5058)
# p1=c(0.0483, 0.4859, 0.5849, 0.6438, 0.6779)


p0=0.316
L1=c(2.5908, 2.5726, 2.5405, 2.5306)
L2=c(2.4075, 2.4615, 2.4784, 2.5194)
p1=c(0.0457, 0.4998, 0.6069, 0.6756, 0.7215)


lambda=0.05
n=c(4,8,16,32)
ARL1=c()
for (j in 1:4){
  for (p1 in c(0.0457, 0.4998, 0.6069, 0.6756, 0.7215) ) {
    RL=EWMA_Sign_L2(n=n[j],L1=L1[j],L2=L2[j],lambda=0.05,p0=p0,p1=p1,set=1)
    ARL1=c(ARL1, round(RL[1],4))
  }
}
matrix=matrix(ARL1, nrow=5, ncol=4, byrow=F)
path=paste0('C:/Users/WEIYU/Desktop/政治大學/品管研究計畫/Meeting/EWMA Sign chart/ARL1/ARL1_BG_p0_',p0,'.csv')
write.csv(matrix,file=path)






#### EVRS chart ARL1 under ARL0=200 ####
#Use different L1, L2, p0 and p1 to calculate the ARL1 under ARL0=200

### Bivariate Skew Normal ARL0=200
# p0=0.317
# L1=c(2.2895, 2.2706, 2.2466, 2.2463)
# L2=c(2.2168, 2.2422, 2.2922, 2.3014)
# p1=c(0.0465, 0.5032, 0.6145, 0.6873, 0.7376)

# p0=0.316
# L1=c(2.2866, 2.2731, 2.2488, 2.2408)
# L2=c(2.2213, 2.2507, 2.2887, 2.3006)
# p1=c(0.0455, 0.5002, 0.6098, 0.6787, 0.7258)

# p0=0.310
# L1=c(2.2908, 2.2778, 2.2533, 2.2446)
# L2=c(2.2131, 2.2553, 2.2697, 2.3004)
# p1=c(0.0458, 0.4879, 0.5903, 0.6518, 0.6874)


### Bivariate gamma ARL0=200
# p0=0.318
# L1=c(2.2863, 2.2704, 2.2480, 2.2454)
# L2=c(2.2324, 2.2517, 2.2807, 2.2984)
# p1=c(0.0476, 0.5053, 0.6159, 0.6886, 0.7389)

# p0=0.316
# L1=c(2.2864, 2.2782, 2.2515, 2.2414)
# L2=c(2.2206, 2.2454, 2.2752, 2.3049)
# p1=c(0.0457, 0.4998, 0.6069, 0.6756, 0.7215)

# p0=0.312
# L1=c(2.2890, 2.2809, 2.2527, 2.2472)
# L2=c(2.2154, 2.2491, 2.2778, 2.2934)
# p1=c(0.0483, 0.4859, 0.5849, 0.6438, 0.6779)


p0=0.312
L1=c(2.2890, 2.2809, 2.2527, 2.2472)
L2=c(2.2154, 2.2491, 2.2778, 2.2934)
p1=c(0.0483, 0.4859, 0.5849, 0.6438, 0.6779)


lambda=0.05
n=c(4,8,16,32)
ARL1=c()
for (j in 1:4){
  for (p1 in c(0.0483, 0.4859, 0.5849, 0.6438, 0.6779) ) {
    RL=EWMA_Sign_L2(n=n[j],L1=L1[j],L2=L2[j],lambda=0.05,p0=p0,p1=p1,set=1)
    ARL1=c(ARL1, round(RL[1],4))
  }
}
matrix=matrix(ARL1, nrow=5, ncol=4, byrow=F)
path=paste0('C:/Users/WEIYU/Desktop/政治大學/品管研究計畫/Meeting/EWMA Sign chart/ARL1/ARL200/ARL1_200_BG_p0_',p0,'.csv')
write.csv(matrix,file=path)




