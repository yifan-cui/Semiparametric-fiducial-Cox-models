
library(survival)
library(Rmosek)
source("coxf_functions_final.R")
nsample <- 40 #30 #20 #10 

nparts = 10 #set to 1 for 20,30, and 
ipart = 10 #set 1:nparts and run each part in parallel


nmcmc <- 1000 #in n=40 #400 #otherwise
nburn <- 100 #in n=40 40 #otherwise

p <- 2
#K <- 20
parT <- c(-0.5,0)
nu <- 1 #h_0(t)=1/2t^{-1/2}
nsimu <- 200
# competing methods
mse1 <- matrix(NA,2,nsimu)
cov1 <- matrix(NA,2,nsimu)
length1 <- matrix(NA,2,nsimu)
mse2 <- matrix(NA,2,nsimu)
cov2 <- matrix(NA,2,nsimu)
length2 <- matrix(NA,2,nsimu)
bd <- 4
flag <- 0  
length.record <- rep(NA,nsimu)
mle1.point <- rep(NA,nsimu)
mle1.upper <- rep(NA,nsimu)
mle1.lower <- rep(NA,nsimu) 
mle2.point <- rep(NA,nsimu) 
mle2.upper <- rep(NA,nsimu)
mle2.lower <- rep(NA,nsimu)
beta1.point <- rep(NA,nsimu)
beta1.upper <- rep(NA,nsimu)
beta1.lower <- rep(NA,nsimu)
beta2.point <- rep(NA,nsimu)
beta2.upper <- rep(NA,nsimu)
beta2.lower <- rep(NA,nsimu)


cpart = nsimu %/% nparts
rpart = nsimu %% nparts 
iindex = 1:(cpart + as.integer(ipart<=rpart))+(ipart-1)*cpart + min(ipart-1, rpart)

for (ii in iindex){
  set.seed(2023*ii)
  print(ii)  
  #X <- matrix(runif(nsample*p),nsample,p)
  flag = TRUE
  while(flag){
    X <- matrix(rbinom(nsample*p,1,1/2),nsample,p)
    T <- cox_gen(X, nu, parT)
    C <- runif(nsample)*2
    Y <- pmin(T,C)
    ind <- (T<=C)
    flag = sum(ind)<=2
  }
  
  model <- coxph(Surv(Y, ind) ~ X)
  beta.mle <- model$coefficients
  mle1.point[ii] <-  model$coefficients[1]  
  mle1.upper[ii] <-  confint(model)[1,2]  
  mle1.lower[ii] <-  confint(model)[1,1]
  mle2.point[ii] <-  model$coefficients[2]  
  mle2.upper[ii] <-  confint(model)[2,2]  
  mle2.lower[ii] <-  confint(model)[2,1]
  
  # summary(model)
  # confint(model) 
  mse1[1,ii] <- (mle1.point[ii]-parT[1])^2
  mse2[1,ii] <- (mle2.point[ii]-parT[2])^2
  cov1[1,ii] <- (mle1.upper[ii]>=parT[1])&(mle1.lower[ii]<=parT[1])
  cov2[1,ii] <- (mle2.upper[ii]>=parT[2])&(mle2.lower[ii]<=parT[2])
  length1[1,ii] <- mle1.upper[ii] - mle1.lower[ii]
  length2[1,ii] <- mle2.upper[ii] - mle2.lower[ii]
  
  
  # Run Gibbs sampler
  q.initial = rep(NA,sum(ind))
  failure.grid0=sort(Y[ind==1])
  for (kk in 1:sum(ind)){
    kkk=which(failure.grid0[kk]==Y) 
    Y.risk0=which(Y>=failure.grid0[kk]) 
    q.initial[kk] = exp(X[kkk,]%*%beta.mle)/sum(exp(X[Y.risk0,]%*%beta.mle))
  }
  u.initial = runif(sum(ind),min=rep(0,sum(ind)),max=q.initial)
  particle <- list(u=u.initial,beta.sup=matrix(NA,ncol=p,nrow=sum(ind)),beta.opt=matrix(NA,ncol=p,nrow=1))
  mcmc.particle <- QGibbs.run(particle,Y,ind,X,nmcmc,nburn,bd)
  beta1.na <- NULL
  beta2.na <- NULL 
  for(nn in 1:nmcmc){
    beta1.na <- c(beta1.na,mcmc.particle[[nn]]$beta.opt[1])
    beta2.na <- c(beta2.na,mcmc.particle[[nn]]$beta.opt[2])
    #print(mcmc.particle[[nn]]$beta.sup[,1])
    #print(mcmc.particle[[nn]]$beta.sup[,2])
  }
  beta1 <- na.omit(beta1.na)
  beta2 <- na.omit(beta2.na)
  length.record[ii] <- length(beta1)
  if(length(beta1)!=0){
    beta1.point[ii] <-  sort(beta1)[length(beta1)*0.5]  
    beta1.upper[ii] <-  sort(beta1)[length(beta1)*0.975]  
    beta1.lower[ii] <-  sort(beta1)[length(beta1)*0.025]  
    beta2.point[ii] <-  sort(beta2)[length(beta2)*0.5]  
    beta2.upper[ii] <-  sort(beta2)[length(beta2)*0.975]  
    beta2.lower[ii] <-  sort(beta2)[length(beta2)*0.025]   
    mse1[2,ii] <- (beta1.point[ii]-parT[1])^2
    mse2[2,ii] <- (beta2.point[ii]-parT[2])^2
    cov1[2,ii] <- (beta1.upper[ii]>=parT[1])&(beta1.lower[ii]<=parT[1])
    cov2[2,ii] <- (beta2.upper[ii]>=parT[2])&(beta2.lower[ii]<=parT[2])
    length1[2,ii] <- beta1.upper[ii]-beta1.lower[ii]
    length2[2,ii] <- beta2.upper[ii]-beta2.lower[ii]
  } 
  else{
    flag <- flag+1
    print(ii*1000)
    mse1[2,ii] <- (mle1.point[ii]-parT[1])^2
    mse2[2,ii] <- (mle2.point[ii]-parT[2])^2
    cov1[2,ii] <- (mle1.upper[ii]>=parT[1])&(mle1.lower[ii]<=parT[1])
    cov2[2,ii] <- (mle2.upper[ii]>=parT[2])&(mle2.lower[ii]<=parT[2])
    length1[2,ii] <- mle1.upper[ii] - mle1.lower[ii]
    length2[2,ii] <- mle2.upper[ii] - mle2.lower[ii]
  }
  
}

rowMeans(mse1)
rowMeans(mse2)
rowMeans(cov1)
rowMeans(cov2)
rowMeans(length1)
rowMeans(length2)

#save(mse1,mse2,cov1,cov2,length1,length2,flag,length.record,file = "cox1n20.RData")
save.image(file=paste("cox1n",nsample,"p",ipart,".RData",sep=""))