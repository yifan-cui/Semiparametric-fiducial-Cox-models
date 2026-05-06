#setwd('~/Dropbox/Documents/Students/Old/UNC2018-2020/CuiYifan/Cox/revision')
library(foreach)
library(doFuture)
#library(progressr)

plan(multisession, workers = 100) # Define your parallel strategy


library(survival)
library(Rmosek)
source("coxf_functions_final.R")
nsample <- 20 #40 #30 #

nmcmc <- 1000 
nburn <- 100 
nsimu <- 200

isetup = 5 #simulation setup

p <- 2
#K <- 20
parT <- c(-1.,-0.5)+isetup*0.5
nu <- 1 

bd <- 4


#with_progress({
  results.list = foreach(ii=1:nsimu,
                         .options.future = list(packages = c("survival", "Rmosek"))) %dofuture%
    {
      set.seed(2023*ii)
      #print(ii)  
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
      mle1.point <-  model$coefficients[1]  
      mle1.upper <-  confint(model)[1,2]  
      mle1.lower <-  confint(model)[1,1]
      mle2.point <-  model$coefficients[2]  
      mle2.upper <-  confint(model)[2,2]  
      mle2.lower <-  confint(model)[2,1]
      
      
      # Run Gibbs sampler
      q.initial = rep(NA,sum(ind))
      failure.grid0=sort(Y[ind==1])
      for (kk in 1:sum(ind)){
        kkk=which(failure.grid0[kk]==Y) 
        Y.risk0=which(Y>=failure.grid0[kk]) 
        q.initial[kk] = exp(X[kkk,]%*%((abs(beta.mle)<bd)*beta.mle))/sum(exp(X[Y.risk0,]%*%
                                                                               ((abs(beta.mle)<bd)*beta.mle)))
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
      length.record <- length(beta1)
      if(length.record != 0){
        beta1.summary <-  quantile(beta1,c(.5,.975,.025))
        beta2.summary <-  quantile(beta2,c(.5,.975,.025))
      }else{
        beta1.summary <- c(mle1.point, mle1.upper, mle1.lower)
        beta2.summary <- c(mle2.point, mle2.upper, mle2.lower)
      }
      c( mle1.point,
         mle1.upper,
         mle1.lower,
         mle2.point,
         mle2.upper,
         mle2.lower,
         beta1.summary,
         beta2.summary, 
         length.record)
    }
#})
  
plan(sequential)#shut down the parallel workers

save.image(file=paste("CoxSetting",isetup,"n",nsample,".RData",sep=""))

results=do.call(cbind, results.list)
mle1.point=results[1,]
mle1.upper=results[2,]
mle1.lower=results[3,]
mle2.point=results[4,]
mle2.upper=results[5,]
mle2.lower=results[6,]
beta1.point=results[7,]  
beta1.upper=results[8,]
beta1.lower=results[9,]
beta2.point=results[10,]
beta2.upper=results[11,]
beta2.lower=results[12,] 
length.record = results[13,]
########

mse1 <- c(mean((mle1.point-parT[1])^2), mean((beta1.point-parT[1])^2))
mse2 <- c(mean((mle2.point-parT[2])^2), mean((beta2.point-parT[2])^2))
cov1 <- c(mean((mle1.upper>=parT[1])&(mle1.lower<=parT[1])),
          mean((beta1.upper>=parT[1])&(beta1.lower<=parT[1])))
cov2 <- c(mean((mle2.upper>=parT[2])&(mle2.lower<=parT[2])),
          mean((beta2.upper>=parT[2])&(beta2.lower<=parT[2])))
length1 <- c(mean(mle1.upper - mle1.lower),mean(beta1.upper-beta1.lower))
length2 <- c(mean(mle2.upper - mle2.lower), mean(beta2.upper-beta2.lower))

save.image(file=paste("CoxSetting",isetup,"n",nsample,".RData",sep=""))

print("beta1")
print(mse1)
print(cov1)
print(length1)

print("beta2")
print(mse2)
print(cov2)
print(length2)

