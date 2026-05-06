setwd('~/Dropbox/Documents/Students/Old/UNC2018-2020/CuiYifan/Cox/revision')
library(survival)
library(Rmosek)
source("coxf_functions_final.R")


select_protocol="HVTN 703/HPTN 081" 
#select_country="South Africa"
select_country="Malawi"
#select_country=c("Malawi","Zimbabwe")


#select_protocol="HVTN 704/HPTN 085"
#select_country="United States"
#select_country="Peru"

nmcmc <- 10000
nburn <- 100
bd <- 10
set.seed(1)
survival.data <- read.csv("survival.csv")
subject.master <- read.csv("subject_master.csv")

merged.subject.data  <- merge( subject.master, survival.data, "pub_id" )
dim( merged.subject.data )
table( merged.subject.data[ , "country" ], merged.subject.data[ , "hiv1event" ] )


my_protocol <- subset(merged.subject.data, protocol.x==select_protocol)
table(my_protocol[ , "country" ], my_protocol[ , "hiv1event80mf" ] )
my_data <- subset(my_protocol, country %in% select_country)
#my_data <- subset(my_protocol, country!=select_country)
#my_data <- my_protocol


X <- as.matrix(ifelse(my_data$rx_pool=="T1+T2",1,0))
Y <- my_data$hiv1survday

#ind <- my_data$hiv1event
ind <- ifelse(my_data$hiv1event80mf==1 | my_data$hiv1event80mf==2 , 1, 0)
ind[is.na(ind)]=0

model <- coxph(Surv(Y, ind) ~ X)
beta.mle <- as.matrix(model$coefficients)
mle.point <-  model$coefficients[1]  
mle.upper <-  confint(model)[1,2]  
mle.lower <-  confint(model)[1,1]

# Run Gibbs sampler
q.initial = rep(NA,sum(ind))
failure.grid0=sort(Y[ind==1])
failure.index0=order(Y)
failure.index0=failure.index0[ind[failure.index0]==1]
for (kk in 1:sum(ind)){
  kkk=failure.index0[kk]
  #kkk=which((failure.grid0[kk]==Y)&(ind==1))
  Y.risk0=which(Y>=failure.grid0[kk]) 
  q.initial[kk] = exp(X[kkk,]%*%((abs(beta.mle)<bd)*beta.mle))/sum(exp(X[Y.risk0,]%*%
                                                                         ((abs(beta.mle)<bd)*beta.mle)))
  #exp(X[kkk,]%*%beta.mle)/sum(exp(X[Y.risk0,]%*%beta.mle))
}
u.initial = runif(sum(ind),min=rep(0,sum(ind)),max=q.initial)
particle <- list(u=u.initial,beta.sup=matrix(NA,ncol=1,nrow=sum(ind)),beta.opt=matrix(NA,ncol=1,nrow=1))
mcmc.particle <- QGibbs.run(particle,Y,ind,X,nmcmc,nburn,bd)

##########
#load(file = paste("Africa-703-",paste0(select_country,collapse = ""),".RData",sep=""))

maxPlot=700
Lambda0.particle <- QLambda0(mcmc.particle,Y,ind,X,replperbeta=1,maxPlot=maxPlot)




beta.na <- NULL
for(nn in 1:nmcmc){
  beta.na <- c(beta.na,mcmc.particle[[nn]]$beta.opt[1])
}
beta <- na.omit(beta.na)

  beta.fid=quantile(beta,c(0.5,0.025,.975))
  beta.point <-  beta.fid[1]
  beta.upper <-  beta.fid[3]
  beta.lower <-  beta.fid[2]
  
  print(c(mle.lower,mle.point,mle.upper))
  print(c(beta.lower,beta.point,beta.upper))
  
  efficacyCI= 1-exp(quantile(beta,c(0.9,.95)))
  print(efficacyCI)
  
  efficacyPoint=1-exp(mean(beta[beta>-bd+.5]))
  print(efficacyPoint)
  
  pvalue=mean(beta>0)
  print(pvalue)
  
  L0timegrid=c(100,200,300,400,500,600)
  
  lambda0_100=numeric(0)
  for(nn in 1:length(Lambda0.particle$samples)){
    temp_func=Lambda0.particle$samples[[nn]]$Lambda0
    if(nn==1){
      curve(temp_func,from=0,to=maxPlot)
    }else{
      curve(temp_func,from=0,to=maxPlot,add=TRUE)
    }
    lambda0_100=c(lambda0_100,temp_func(100))
  }
  
  Lambda0_grid=numeric(0)
  for(nn in 1:length(Lambda0.particle$samples)){
    temp_func=Lambda0.particle$samples[[nn]]$Lambda0
    Lambda0_grid=rbind(Lambda0_grid,temp_func(L0timegrid))
  }
  #save.image("temp.RData")
  Lambda0_grid.fid=apply(Lambda0_grid,2,quantile,c(0.5,0.025,.975))
  Lambda0.point <-  Lambda0_grid.fid[1,]
  Lambda0.upper <-  Lambda0_grid.fid[3,]
  Lambda0.lower <-  Lambda0_grid.fid[2,]
  
  print(rbind(L0timegrid,Lambda0.upper,Lambda0.point,Lambda0.lower))
  
  save.image(file = paste("Africa-703-",paste0(select_country,collapse = ""),".RData",sep=""))
  #save.image(file = paste("Africa-703-Not",select_country,".RData",sep=""))
  
  
  
