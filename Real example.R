# Real example

#### Semicinductor Data ####

realdata_IC=read.csv("C:/Users/WEIYU/Desktop/政治大學/品管研究計畫/Meeting/Real example/secom_icdata.csv")
realdata_OC=read.csv("C:/Users/WEIYU/Desktop/政治大學/品管研究計畫/Meeting/Real example/secom_ocdata.csv")

#### IC_R and OC_R ####
#Calculate IC_R and OC_R 

# IC_R
U=realdata_IC$V477[1:240]
V=realdata_IC$V478[1:240]
R=matrix(U/V, ncol=4, byrow=T)
R=data.frame(num=c(1:60),R1=R[,1],R2=R[,2],R3=R[,3],R4=R[,4])
#write.csv(R,'C:/Users/WEIYU/Desktop/政治大學/品管研究計畫/Meeting/實際資料/Excel/IC_R.csv')

# OC_R after adjusting the order of data
new=c(
  0.07860579, 0.10202354, 0.05659306, 0.08386750, 42107,
  0.55468822, 0.29772874, 0.42623377, 0.36413102, 33729,
  0.46537446, 0.12403892, 0.69106590, 0.29906535, 34097,
  0.33115259, 0.77457120, 0.16182476, 1.22224427, 35673,
  0.48358243, 0.23174015, 0.13928095, 0.17295723, 16087,
  0.42817402, 0.19466018, 0.85372395, 0.33389512, 30579,
  0.26141629, 0.40270675, 0.58530584, 0.36518171, 30853,
  0.34674190, 0.64029983, 0.44350342, 0.27911949, 33229,
  0.26673691, 0.08808623, 0.06852969, 0.14376873, 26905,
  0.35415184, 0.37656826, 0.22331038, 0.05079054, 28485,
  0.78588221, 0.11424528, 0.13439937, 0.30219238, 28303,
  0.32906246, 0.24324360, 0.39713349, 0.26657092, 17853,
  0.32370121, 0.32761341, 0.13632644, 0.38482665, 22863,
  0.18365746, 0.69065261, 0.19579510, 0.53735346, 26555,
  0.16661412, 0.23005886, 0.45199169, 0.85846473, 26345,
  0.11358504, 0.20322167, 0.21226386, 0.48363453, 17837,
  0.09968841, 0.10609497, 0.14762340, 0.21354348, 15387,
  0.10684180, 0.16184980, 0.58814305, 0.22056345, 20293,
  0.31502327, 0.07066136, 0.12495762, 0.21158368, 21923,
  0.11323926, 0.43881921, 0.28452358, 0.23309248, 20391,
  0.19612648, 0.23911530, 0.13508241, 0.41880955, 14109,
  0.14045102, 0.25982709, 0.31778623, 0.18347700, 10485,
  0.16713619, 0.23653972, 0.38807586, 0.21445525,  9981,
  0.12178084, 0.20755143, 0.24436652, 0.27420596,  9213,
  0.23776349, 0.14760573, 0.10883630, 0.19082319,  8945,
  0.28403129, 0.18916601, 0.19370839, 0.17643720,  4459
)

new_R_oc=matrix(new,ncol=5,byrow=T)[,1:4]
#### Zero state EVRS chart ####
#Plot the IC and OC zero state EVRS chart 

##   IC part
U=realdata_IC$V477[1:240]
V=realdata_IC$V478[1:240]
R=U/V  
Y=c()
for (i in c(1:(length(R)/2))) {
  Y[i]=(R[2*i-1]-R[2*i])^2/2
}
Y=matrix(Y,ncol=2,byrow=T)
p0=sum(Y>var(R))/length(Y) #p0=0.2333, var(R)=0.01994
I=(Y>var(R))*1
B=rowSums(I) #B~Bin(0.5n,p0)
EWMA=c()
Z=c()
n=4;lambda=0.05
mean=0.5*n*p0
var=0.5*n*p0*(1-p0)
EWMA[1]=lambda*B[1]+(1-lambda)*mean
Z[1]=(EWMA[1]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^2))
for (i in 2:length(B)) {
  EWMA[i]=lambda*B[i]+(1-lambda)*EWMA[i-1]
  Z[i]=(EWMA[i]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^(2*i)))
}

#作圖
interpolation=function(x,up,low,up_L,low_L){
  L=low_L+(up_L-low_L)*(x-low)/(up-low) #內差法interpolation反推L
  return(round(L,4))
}
L1=interpolation(x=p0,up=0.2,low=0.3,up_L=2.7167,low_L=2.5995)
L2=interpolation(x=p0,up=0.2,low=0.3,up_L=2.3035,low_L=2.4004)
#L1=2.6300 
#L2=2.3752 
Z_col=Z
Z_col[Z<L1 | Z>(-L2)]=1
Z_col[Z>L1 | Z<(-L2)]=2
num=c(1:length(Z))
plot(num,      # X軸
     Z,       # Y軸
     main="",    # 圖表名稱
     xlab="Samples",    # X軸名稱
     ylab="EVRS Z statistic",
     ylim=c(-4,5),
     type="b",
     col=Z_col
)
abline(h =L1, untf = FALSE)
text(50,L1+0.5,paste("UCL=",L1))
abline(h =-L2, untf = FALSE)
text(50,-L2-0.5,paste("LCL=",-L2))


##  OC part
U_oc=realdata_OC$V477
V_oc=realdata_OC$V478
R_oc=U_oc/V_oc
#用new_R_oc
y1=(new_R_oc[,1]-new_R_oc[,2])^2/2
y2=(new_R_oc[,3]-new_R_oc[,4])^2/2
Y_oc=matrix(c(y1, y2),ncol=2,byrow=F)
p1=sum(Y_oc>var(R))/length(Y_oc) #0.3846
I_oc=(Y_oc>var(R))*1
B_oc=rowSums(I_oc)
EWMA_oc=c()
Z_oc=c()
n=4;lambda=0.05
mean=0.5*n*p0
var=0.5*n*p0*(1-p0)
EWMA_oc[1]=lambda*B_oc[1]+(1-lambda)*mean
Z_oc[1]=(EWMA_oc[1]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^2))
for (i in 2:length(B_oc)) {
  EWMA_oc[i]=lambda*B_oc[i]+(1-lambda)*EWMA_oc[i-1]
  Z_oc[i]=(EWMA_oc[i]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^(2*i)))
}

# 作圖
interpolation=function(x,up,low,up_L,low_L){
  L=low_L+(up_L-low_L)*(x-low)/(up-low) #內差法interpolation反推L
  return(round(L,4))
}
L1=interpolation(x=p0,up=0.2,low=0.3,up_L=2.7167,low_L=2.5995)
L2=interpolation(x=p0,up=0.2,low=0.3,up_L=2.3035,low_L=2.4004)


Z_col=Z_oc
Z_col[Z_oc<L1 | Z_oc>(-L2)]=1
Z_col[Z_oc>L1 | Z_oc<(-L2)]=2
num=c(1:length(Z_oc))
plot(num,      # X軸
     Z_oc,       # Y軸
     main="",    # 圖表名稱
     xlab="Samples",    # X軸名稱
     ylab="EVRS Z statistic",
     ylim=c(-4,5),
     type="b",
     col=Z_col
)
abline(h =L1, untf = FALSE)
text(20,L1+1.5,paste("UCL=",L1))
abline(h =-L2, untf = FALSE)
text(20,-L2+0.5,paste("LCL=",-L2))
#### Steady state EVRS chart ####
#Plot the steady state EVRS chart 

n=4;lambda=0.05;
U=realdata_IC$V477[1:240]
V=realdata_IC$V478[1:240]
R=matrix(U/V,ncol=4,byrow=T)
R_pool_R_oc=rbind(R,new_R_oc) #Combine IC 1~240 OC 241~344 data
y1=(R_pool_R_oc[,1]-R_pool_R_oc[,2])^2/2
y2=(R_pool_R_oc[,3]-R_pool_R_oc[,4])^2/2
Y=matrix(c(y1, y2),ncol=2,byrow=F)
p0=sum(Y[1:60,]>var(U/V))/length(Y[1:60,]) #0.2333
I=(Y>var(U/V))*1
B=rowSums(I) 

EWMA=c()
Z=c()
mean=0.5*n*p0
var=0.5*n*p0*(1-p0)
EWMA[1]=lambda*B[1]+(1-lambda)*mean
Z[1]=(EWMA[1]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^2))
for (i in 2:length(B)) {
  EWMA[i]=lambda*B[i]+(1-lambda)*EWMA[i-1]
  Z[i]=(EWMA[i]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^(2*i)))
}

#作圖
interpolation=function(x,up,low,up_L,low_L){
  L=low_L+(up_L-low_L)*(x-low)/(up-low) #內差法interpolation反推L
  return(round(L,4))
}
L1=interpolation(x=p0,up=0.2,low=0.3,up_L=2.7167,low_L=2.5995)
L2=interpolation(x=p0,up=0.2,low=0.3,up_L=2.3035,low_L=2.4004)

Z_col=Z
Z_col[Z<L1 | Z>(-L2)]=1
Z_col[Z>L1 | Z<(-L2)]=2
num=c(1:length(Z))
plot(num,      # X軸
     Z,       # Y軸
     main="",    # 圖表名稱
     xlab="Samples",    # X軸名稱
     ylab="EVRS Z statistic",
     ylim=c(-4,5),
     type="b",
     col=Z_col
)
abline(h =L1, untf = FALSE)
text(20,L1+0.5,paste("UCL=",L1))
abline(h =-L2, untf = FALSE)
text(20,-L2-0.5,paste("LCL=",-L2))
abline(v=60,col="red")

#### Zero state EVRM chart ####
#Plot the IC and OC zero state EVRM chart 

##   IC part
# New IC R data (random sample)
U=realdata_IC$V477[1:240]
V=realdata_IC$V478[1:240]
R=U/V
set.seed(123)
new_R_ic=matrix(NA, ncol=4,nrow=60)
for (i in 1:60) {
  new_R_ic[i,]=sample(R,4,replace = F)
}

# 求出 IC EWMA、ZEWMA
N=240
n=4
N_pool=N+n
mean=n*(N_pool^2-1)/12
var=N*n*(N_pool+1)*(N_pool^2-4)/180
lambda=0.05
M=c()
EWMA=c()
Z1=c()
R_pool=sort(c(R, new_R_ic[1,]))
r=c()
for (j in 1:4) {
  a=which(R_pool==new_R_ic[1,j])
  r=c(r, mean(a))
}
M[1]=sum((r-(N_pool+1)/2)^2)
EWMA[1]=lambda*M[1]+(1-lambda)*mean
Z1[1]=(EWMA[1]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^2))
for (i in 2:60) {
  R_pool=sort(c(R, new_R_ic[i,]))
  r=c()
  for (j in 1:4) {
    a=which(R_pool==new_R_ic[i,j])
    r=c(r, mean(a))
  }
  M[i]=sum((r-(N_pool+1)/2)^2)
  EWMA[i]=lambda*M[i]+(1-lambda)*EWMA[i-1]
  Z1[i]=(EWMA[i]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^(2*i)))
}

#找L1、L2 畫圖
interpolation=function(x,up,low,up_L,low_L){
  L=up_L-(up_L-low_L)*(up-x)/(up-low) #內差法interpolation反推L
  return(round(L,4))
}
L1=interpolation(x=240,up=100,low=500,up_L=2.548,low_L=2.551)
L2=interpolation(x=240,up=100,low=500,up_L=2.490,low_L=2.470)

Z_col=Z1
Z_col[Z1<L1]=1
Z_col[Z1>L1]=2
num=c(1:60)
plot(num,      # X軸
     Z1,       # Y軸
     xlab="Samples",    # X軸名稱
     ylab="EVRM Z statistic",
     ylim=c(-4,5),
     type="b",
     col=Z_col
)
abline(h =L1, untf = FALSE)
text(22,L1+0.6,paste("UCL=",L1))
abline(h =-L2, untf = FALSE)
text(22,-L2-0.5,paste("LCL=",-L2))




##   OC part
# 求出 OC EWMA、ZEWMA
N=240
n=4
lambda=0.05
M=c()
EWMA=c()
Z1=c()
R_pool=sort(c(R, new_R_oc[1,]))
N_pool=N+n
mean=n*(N_pool^2-1)/12
var=N*n*(N_pool+1)*(N_pool^2-4)/180
b1=c()
for (count in 1:N_pool) {  
  if (R_pool[count] %in% new_R_oc[1,]) {
    b1=c(b1,(count-(N_pool+1)/2)^2) 
  }
}
M[1]=sum(b1)
EWMA[1]=lambda*M[1]+(1-lambda)*mean
Z1[1]=(EWMA[1]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^2))
for (i in 2:nrow(new_R_oc)) {
  R_pool=sort(c(R, new_R_oc[i,]))
  b1=c()
  for (count in 1:N_pool) {  
    if (R_pool[count] %in% new_R_oc[i,]) {
      b1=c(b1,(count-(N_pool+1)/2)^2) 
    }
  }
  M[i]=sum(b1)
  EWMA[i]=lambda*M[i]+(1-lambda)*EWMA[i-1]
  Z1[i]=(EWMA[i]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^(2*i)))
}

#找L1、L2
interpolation=function(x,up,low,up_L,low_L){
  L=up_L-(up_L-low_L)*(up-x)/(up-low) #內差法interpolation反推L
  return(round(L,4))
}
L1=interpolation(x=240,up=100,low=500,up_L=2.548,low_L=2.551)
L2=interpolation(x=240,up=100,low=500,up_L=2.490,low_L=2.470)

Z_col=Z1
Z_col[Z1<L1]=1
Z_col[Z1>L1]=2
num=c(1:26)
plot(num,      # X軸
     Z1,       # Y軸
     xlab="Samples",    # X軸名稱
     ylab="EVRM Z statistic",
     ylim=c(-4,5),
     type="b",
     col=Z_col
)
abline(h =L1, untf = FALSE)
text(22,L1+0.6,paste("UCL=",L1))
abline(h =-L2, untf = FALSE)
text(22,-L2-0.5,paste("LCL=",-L2))

#### Zero state EVRT chart ####
#Plot the IC and OC zero state EVRT chart 

# R、R*
U=realdata_IC$V477[1:240] #n=240 because we need sample group m=60
V=realdata_IC$V478[1:240]
R=U/V

##   IC part
# New IC R data (random sample)　
set.seed(123)
new_R_ic=matrix(NA, ncol=4,nrow=60)
for (i in 1:60) {
  new_R_ic[i,]=sample(R,4,replace = F)
}

# 求出EWMA、ZEWMA
generate_rank=function(N_pool){
  a=c()  #此rank只適用偶數、且N_pool不能太小 QQ
  a[1]=1
  for (i in seq(2,(N_pool-1)/2,by=2)) {
    a[i]=a[i-1]+3
    a[i+1]=a[i]+1
  }
  
  b=c()
  b[1]=2
  b[2]=3
  for (j in seq(3,(N_pool-1)/2,by=2)) {
    b[j]=b[j-1]+3
    b[j+1]=b[j]+1
  }
  ans=c(a,N_pool,sort(b,decreasing = TRUE))
  return(ans)
}
N=240;n=4;lambda=0.05
N_pool=N+n
rank=generate_rank(N_pool)
mean=N*(N+n+1)/2
var=N*(N_pool^2-1)*(N_pool-N) / (12*(N_pool-1))
m=matrix(NA,ncol=N_pool,nrow=2)
Tukey=c()
EWMA=c()
Z2=c()
for (i in 1:60) {
  R_pool=sort(c(R,new_R_ic[i,]))
  m[1,]=R_pool
  m[2,]=rank
  r=c()
  for (j in 1:4) {
    a=m[2,(m[1,] %in% new_R_ic[i,j])]
    r=c(r, mean(a))
  }
  Tukey[i]=sum(m[2,])-sum(r)
}
EWMA[1]=lambda*Tukey[1]+(1-lambda)*mean
Z2[1]=(EWMA[1]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^2))
for (j in 2:60) {
  EWMA[j]=lambda*Tukey[j]+(1-lambda)*EWMA[j-1]
  Z2[j]=(EWMA[j]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^(2*j)))
}

#求出L1、L2 畫圖
interpolation=function(x,up,low,up_L,low_L){
  L=up_L-(up_L-low_L)*(up-x)/(up-low) #內差法interpolation反推L
  return(round(L,4))
}
L1=interpolation(x=240,up=100,low=500,up_L=2.485,low_L=2.493)
L2=interpolation(x=240,up=100,low=500,up_L=2.554,low_L=2.559)

Z_col=Z2
Z_col[Z2<L1]=1
Z_col[Z2>L1]=2
num=c(1:60)
plot(num,      # X軸
     Z2,       # Y軸
     xlab="Samples",    # X軸名稱
     ylab="EVRT Z statistic",
     ylim=c(-4,5),
     type="b",
     col=Z_col
)
abline(h =L1, untf = FALSE)
text(22,L1+0.5,paste("UCL=",L1))
abline(h =-L2, untf = FALSE)
text(22,-L2-0.5,paste("LCL=",-L2))





##   OC part
# 求出EWMA、ZEWMA
generate_rank=function(N_pool){
  a=c()  #此rank只適用偶數、且N_pool不能太小 QQ
  a[1]=1
  for (i in seq(2,(N_pool-1)/2,by=2)) {
    a[i]=a[i-1]+3
    a[i+1]=a[i]+1
  }
  
  b=c()
  b[1]=2
  b[2]=3
  for (j in seq(3,(N_pool-1)/2,by=2)) {
    b[j]=b[j-1]+3
    b[j+1]=b[j]+1
  }
  ans=c(a,N_pool,sort(b,decreasing = TRUE))
  return(ans)
}
N=240;n=4;lambda=0.05
N_pool=N+n
rank=generate_rank(N_pool)
mean=N*(N+n+1)/2
var=N*(N_pool^2-1)*(N_pool-N) / (12*(N_pool-1))
m=matrix(NA,ncol=N_pool,nrow=2)
Tukey=c()
EWMA=c()
Z2=c()
for (i in 1:26) {
  R_pool=sort(c(R,new_R_oc[i,]))
  m[1,]=R_pool
  m[2,]=rank
  Tukey[i]=sum(m[2,(m[1,] %in% R)])
}
EWMA[1]=lambda*Tukey[1]+(1-lambda)*mean
Z2[1]=(EWMA[1]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^2))
for (j in 2:26) {
  EWMA[j]=lambda*Tukey[j]+(1-lambda)*EWMA[j-1]
  Z2[j]=(EWMA[j]-mean)/sqrt(var*lambda/(2-lambda)*(1-(1-lambda)^(2*j)))
}

#求出L1、L2 畫圖
interpolation=function(x,up,low,up_L,low_L){
  L=up_L-(up_L-low_L)*(up-x)/(up-low) #內差法interpolation反推L
  return(round(L,4))
}
L1=interpolation(x=240,up=100,low=500,up_L=2.485,low_L=2.493)
L2=interpolation(x=240,up=100,low=500,up_L=2.554,low_L=2.559)

Z_col=Z2
Z_col[Z2<L1]=1
Z_col[Z2>L1]=2
num=c(1:26)
plot(num,      # X軸
     Z2,       # Y軸
     xlab="Samples",    # X軸名稱
     ylab="EVRT Z statistic",
     ylim=c(-4,5),
     type="b",
     col=Z_col
)
abline(h =L1, untf = FALSE)
text(22,L1+0.5,paste("UCL=",L1))
abline(h =-L2, untf = FALSE)
text(22,-L2-0.5,paste("LCL=",-L2))
