cox_gen=function(predX, nu, parT){
  m=dim(predX)[2]
  NN=dim(predX)[1]
  #Z=cbind(matrix(1,NN,1), predX)
  Z=predX
  UT=runif(NN)
  TT=(-log(UT)/exp(Z%*%parT))^nu
  return(TT)
}

QGibbs.single<-function(particle,i,Y,ind,X,bd){
  #n=length(Y)
  #Y.sort=sort(Y)
  failure.grid=sort(Y[ind==1])
  m=sum(ind)
  preparticle=particle
  preparticle$u=particle$u[-i]  # i from 1 to m: ordered deaths
  
  pp=dim(X)[2]
  
  failure.index0=order(Y)
  failure.index0=failure.index0[ind[failure.index0]==1]
  iii=failure.index0[i]
  X.die = as.matrix(X[ind==1,])########
  X.die.reorder = as.matrix(X.die[order(Y[ind==1]),])########
  X.die.reorder.exi = as.matrix(X.die.reorder[-i,])########
  
  num.atrisk=rep(NA,m)
  for(k in 1:m){
    num.atrisk[k]=sum(Y>=failure.grid[k]) 
  }
  
  total.num.atrisk = sum(num.atrisk)
  
  ## Cone optimization 
  
  prob <- list(sense="min")
  # Objective coefficients   #  beta; s; t
  prob$c <- c(-X[iii,], diag(m)[i,], rep(0,total.num.atrisk))
  
  # Specify matrix 'A' 
  r1 <- cbind(X.die.reorder.exi,matrix(-diag(m)[-i,],nrow=m-1),matrix(0,m-1,total.num.atrisk))
  temp.t <- matrix(0,m,total.num.atrisk)
  temp.t[1,1:num.atrisk[1]] <- 1
  for(k in 2:m){
    temp.t[k,(sum(num.atrisk[1:(k-1)])+1):sum(num.atrisk[1:k])] <- 1
  }
  r2 <- cbind(matrix(0,m,pp+m),temp.t)
  prob$A <- Matrix(rbind(r1,r2))
  
  # Bound values for constraints
  prob$bc <- rbind(blc=c(log(preparticle$u),  rep(0,m)), 
                   buc=c(rep(bd-log(length(Y)*exp(-bd)),m-1), rep(1,m)))
  
  #print(prob$bc)
  #print(preparticle$u)
  
  # Bound values for variables
  prob$bx <- rbind(blx=c(rep(-bd,pp), rep(log(length(Y)*exp(-bd)),m),rep(0,total.num.atrisk)), 
                   bux=c(rep(bd,pp), rep(log(length(Y)*exp(bd)),m),rep(1,total.num.atrisk)))
  
  #print(prob$bx)
  # Specify the affine conic constraints.
  FE <- Matrix(nrow=0, ncol = pp+m+total.num.atrisk)
  
  Y.risk=which(Y>=failure.grid[1]) 
  for(jj in 1:num.atrisk[1]) {
    temp <- integer(total.num.atrisk)
    k <- jj
    temp[k] <- 1
    FE <- rbind(FE,
                rbind(c(rep(0,pp+m),temp),rep(0,pp+m+total.num.atrisk),c(X[Y.risk[jj],],-diag(m)[1,],rep(0,total.num.atrisk))) )
  }
  for(ll in 2:m) {
    Y.risk=which(Y>=failure.grid[ll]) 
    for(jj in 1:num.atrisk[ll]) {
      temp <- integer(total.num.atrisk)
      k <- sum(num.atrisk[1:(ll-1)])+jj
      temp[k] <- 1
      FE <- rbind(FE,
                  rbind(c(rep(0,pp+m),temp),rep(0,pp+m+total.num.atrisk),c(X[Y.risk[jj],],-diag(m)[ll,],rep(0,total.num.atrisk))) )
    }
  }  
  gE <- rep(c(0, 1, 0), total.num.atrisk)
  
  prob$F <- Matrix(FE)
  prob$g <- gE
  prob$cones <- matrix(list("PEXP", 3, NULL), nrow=3, ncol=total.num.atrisk)
  rownames(prob$cones) <- c("type","dim", "conepar")
  
  
  # Solve the problem
  invisible(capture.output(r <- mosek(prob)))
  #r <- mosek(prob, list(verbose=0))
  
  #print(i)
  # if (r$sol$itr$prosta == 'PRIMAL_INFEASIBLE') {
  # 
  #   print("primal infeasibility")
  # }
  
  if (r$sol$itr$prosta == 'DUAL_INFEASIBLE') {
    
    print("dual infeasibility")
  }
  
  # Return the solution
  #stopifnot(identical(r$response$code, 0))
  
  if (r$sol$itr$prosta == 'PRIMAL_INFEASIBLE') {
    
    print("primal infeasibility")
    beta.sup=rep(NA,pp)
  }
  else{
    solu=r$sol$itr$xx
    beta.sup=as.matrix(solu[1:pp])#######
    Y.risk=which(Y>=failure.grid[i]) 
    q.sup=exp(X[iii,]%*%beta.sup)/sum(exp(X[Y.risk,]%*%beta.sup))
    particle$u[i]=runif(1)*q.sup}
  #print(beta.sup)
  #print(q.sup) 
  
  
  particle$beta.sup[i,]=beta.sup
  return(particle)
}

QGibbs.rscan<-function(particle,Y,ind,X,bd){
  
  for (i in 1:(sum(ind))){
    particle=QGibbs.single(particle,i,Y,ind,X,bd)
  }
  
  
  
  failure.grid=sort(Y[ind==1])
  m=sum(ind)
  pp=dim(X)[2]
  w=rnorm(pp,0,1)
  X.die = as.matrix(X[ind==1,]) #########
  X.die.reorder = as.matrix(X.die[order(Y[ind==1]),]) #######
  
  
  num.atrisk=rep(NA,m)
  for(k in 1:m){
    num.atrisk[k]=sum(Y>=failure.grid[k]) 
  }
  
  total.num.atrisk = sum(num.atrisk)
  
  ## Cone optimization 
  
  prob <- list(sense="min")
  # Objective coefficients   #  beta; s; t
  prob$c <- c(-w, rep(0,m+total.num.atrisk))
  
  # Specify matrix 'A' 
  r1 <- cbind(X.die.reorder,-diag(m),matrix(0,m,total.num.atrisk))
  temp.t <- matrix(0,m,total.num.atrisk)
  temp.t[1,1:num.atrisk[1]] <- 1
  for(k in 2:m){
    temp.t[k,(sum(num.atrisk[1:(k-1)])+1):sum(num.atrisk[1:k])] <- 1
  }
  r2 <- cbind(matrix(0,m,pp+m),temp.t)
  prob$A <- Matrix(rbind(r1,r2))
  
  
  prob$bc <- rbind(blc=c(log(particle$u),  rep(0,m)), 
                   buc=c(rep(bd-log(length(Y)*exp(-bd)),m), rep(1,m)))
  
  prob$bx <- rbind(blx=c(rep(-bd,pp), rep(log(length(Y)*exp(-bd)),m),rep(0,total.num.atrisk)), 
                   bux=c(rep(bd,pp), rep(log(length(Y)*exp(bd)),m),rep(1,total.num.atrisk)))
  
  FE <- Matrix(nrow=0, ncol = pp+m+total.num.atrisk)
  
  Y.risk=which(Y>=failure.grid[1]) 
  for(jj in 1:num.atrisk[1]) {
    temp <- integer(total.num.atrisk)
    k <- jj
    temp[k] <- 1
    FE <- rbind(FE,
                rbind(c(rep(0,pp+m),temp),rep(0,pp+m+total.num.atrisk),c(X[Y.risk[jj],],-diag(m)[1,],rep(0,total.num.atrisk))) )
  }
  for(ll in 2:m) {
    Y.risk=which(Y>=failure.grid[ll]) 
    for(jj in 1:num.atrisk[ll]) {
      temp <- integer(total.num.atrisk)
      k <- sum(num.atrisk[1:(ll-1)])+jj
      temp[k] <- 1
      FE <- rbind(FE,
                  rbind(c(rep(0,pp+m),temp),rep(0,pp+m+total.num.atrisk),c(X[Y.risk[jj],],-diag(m)[ll,],rep(0,total.num.atrisk))) )
    }
  }  
  gE <- rep(c(0, 1, 0), total.num.atrisk)
  
  prob$F <- Matrix(FE)
  prob$g <- gE
  prob$cones <- matrix(list("PEXP", 3, NULL), nrow=3, ncol=total.num.atrisk)
  rownames(prob$cones) <- c("type","dim", "conepar")
  
  
  # Solve the problem
  invisible(capture.output(r <- mosek(prob)))
  #r <- mosek(prob, list(verbose=0))
  
  #print(i)
  # if (r$sol$itr$prosta == 'PRIMAL_INFEASIBLE') {
  # 
  #   print("primal infeasibility")
  # }
  
  if (r$sol$itr$prosta == 'DUAL_INFEASIBLE') {
    
    print("dual infeasibility")
  }
  
  # Return the solution
  #stopifnot(identical(r$response$code, 0))
  
  if (r$sol$itr$prosta == 'PRIMAL_INFEASIBLE') {
    
    print("primal infeasibility")
    beta.opt=rep(NA,pp)
  }
  else{
    solu=r$sol$itr$xx
    beta.opt=solu[1:pp]}
  
  particle$beta.opt=beta.opt
  
  
  return(particle)
}

QGibbs.run<-function(particle,Y,ind,X,nmcmc,nburn=0,bd){
  v.particle = list()
  for(i in 1:nburn){
    #print(i-nburn)
    particle = QGibbs.rscan(particle,Y,ind,X,bd)
    # sample beta subject U\leq q
  }
  for(i in 1:nmcmc){
    #print(i)
    particle = QGibbs.rscan(particle,Y,ind,X,bd)
    v.particle[[i]] = particle
  }
  return(v.particle)
}

QLambda0<-function(particle,Y,ind,X,replperbeta=1,maxPlot=NA){
  #get ordered failure times
  t=unique(sort(c(0,Inf,Y[which(ind==1)])))
  m=length(t)-1
  ltime=t(diff(outer(t,as.vector(Y),pmin)))
  
  approxt=pmin(t,max(maxPlot,2*max(Y),na.rm=TRUE))
  
  lambda0list <- list(t=t,X=X,Y=Y,ind=ind, samples=list())
  
  for (i in 1:length(particle)){
    #print(i)
    beta=particle[[i]]$beta.opt
    if (sum(is.na(beta))==0){
      imultiple=exp(X%*%beta)
      lmat=sweep(ltime, 1, imultiple, FUN = "*")
      li=colSums(lmat)
      
      liuse=li
      liuse[m]=max(li[m-1],2*li[m])
      
      
      for (j in 1:replperbeta){
        
        lambda0=rexp(liuse)
        
        approxL0=cumsum(c(0,diff(approxt)*lambda0))
        Lambda0=approxfun(approxt,approxL0,method = "linear")
        
        lambda0list$samples <- c(lambda0list$samples,list(list(beta=beta,lambda0=lambda0,Lambda0=Lambda0)))
      }
    }
  }
  return(lambda0list)
}