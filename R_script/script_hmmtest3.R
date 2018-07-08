library("rstan")
#setwd("~/project/data")

#data cleaning

#read data, log of tree ring width, 247 obs of 500 variables
logy=read.table("logy.txt",head=FALSE)

#read year, 1 obs of 500 variables 1496-1995
year=read.table("year.txt",head=FALSE)

#read tree age, 247 obs of 500 variables
age=read.table("age.txt",head=FALSE)

#read temperature, 1 obs of 500 variables, first 417 are missing
temp=read.table("temp.txt",head=FALSE)

#rescale temperature, centering at 0 with sd 0.5
stemp=scale(as.numeric(temp))/2 

#substitute missing value as 0
mlogy=logy
mlogy[is.na(logy)]=0

#rescale age, substitute missing value as 0
sage=apply(age,2,function(x) (x-100)/100)
msage=sage
msage[is.na(msage)]=0

#start and end index of tree ring observation, vector of 247 elements
s=apply(sage,1,function(x) min(which(is.na(x)=="FALSE")))
e=apply(sage,1,function(x) max(which(is.na(x)=="FALSE")))



# getrowIndex=function(x){
#   index=which(x!=0)
#   index
# }
#nt_1=which(e>=418)
#get tree index per year
indexpool=apply(age[,1:500],2,function(x) {which(x!=0)})
#get how many trees' record available per year
len=as.numeric(do.call(rbind, lapply(indexpool, length)))
# get index matrix
b=matrix(0,nrow=247,500)
for(i in 1:500){
  b[1:len[i],i]=unlist(indexpool[i]) }



#preparation for rstan
standata=list(K=3,span=500,nt=247,start=s,end=e,
              age=as.matrix(msage),N_mis=417,y=as.matrix(mlogy),
              x_obs=as.vector(stemp)[418:500],alpha=c(1,1,1),b=b, len=len)

inits_stan=function(){
  mu_x=rnorm(3,0,0.5)
  sigmay=runif(1,0.3,0.4)
  sigmax=runif(1,0.4,0.6)
  sigbeta0=runif(1,0.4,0.6)
  #sigbeta1=runif(1,0.5,0.7)
  sigbeta2=runif(1,0.6,0.8)
  sigmaeta=runif(1,0.3,0.4)
  betamu0=runif(1,6,6.5)
  beta0=rnorm(247,6.5,0.5)
  betamu1=runif(1,-0.5,-0.3)
  #beta1=rnorm(247,-0.4,0.1)
  betamu2=runif(1,0.4,0.5)
  PI=c(1/3,1/3,1/3)
  list(mu_x=mu_x,ky=sigmay,kx=sigmax,sigbeta0=sigbeta0,
  #sigbeta1=sigbeta1,
       sigbeta2=sigbeta2,keta=sigmaeta,betamu0=betamu0,betamu1=betamu1,
	   #beta1=beta1,
	   beta0=beta0,betamu2=betamu2,PI=PI
  )
}
inits=list()
for(i in 1:3)inits[[i]] = inits_stan()

hmm.test3<-stan(file="hmm_test3.stan",data=standata,chains=3,init=inits,iter=10000)

save(hmm3.test,file="hmmtest3.Rdata")
# 
# save(stanfit5,file="stanfit5.Rdata")
# 
# stan2coda <- function(fit) {
#   mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
# }
# 
# mcmcfit5=stan2coda(stanfit5)


# plot(stanfit3,par=c("lp__","mux[1]"))
# 
# params1 <- as.data.frame(extract(stanfit3, permuted=FALSE)[,1,])
# params2<-as.data.frame(extract(stanfit3, permuted=FALSE)[,2,])
# params3<-as.data.frame(extract(stanfit3, permuted=FALSE)[,3,])
# par(mar = c(4, 4, 0.5, 0.5))
# plot(params1$"mux[2]", params1$"mux[3]", col=2, pch=16, cex=0.8,
#      xlab="mu1", xlim=c(-3, 3), ylab="mu2", ylim=c(-3, 3))
# lines(0.08*(1:100) - 4, 0.08*(1:100) - 4, col="grey", lw=2)
# params1
#   