
library(Rcpp) ;library(rrpack)
sourceCpp('test.cpp')
tau = .1  # in the prior
Iters = 15000
burnin = 1000
# random data generation
n = 100  # samples
ntest = n*.2
l = 8   # response
p = 12   # predictors
r = 2    # true rank
simdata <- rrr.sim3(n = n + ntest , p = p, q.mix = c(0, 0, l), intercept = rep(0,l), nrank = r, mis.prop = 0)
family <- simdata$family 
control = list(epsilon = 1e-4, sv.tol = 1e-2, maxit = 1000, trace =FALSE,gammaC0 = 1.1, plot.cv = F,conv.obj=TRUE)
alpha = 0.99

mala = lmc = rrr =  list(); acceptrate = c()


for (ss in 72:100) {
  simdata <- rrr.sim3(n = n + ntest , p = p, q.mix = c(0, 0, l), intercept = rep(0,l), nrank = r, mis.prop = 0)
  X = simdata$X[1:n,]; tX = t(X)
  C = simdata$C
  Y = simdata$Y[1:n,] ; tY = t(Y); tXY = tX%*%Y
  xtest = simdata$X[-(1:n),] ; ytest = simdata$Y[-(1:n) ,]
  
  ### mRRR
  fit.cv.mrrr <- cv.mrrr(Y , X, family = family,control = control, penstr = list(penaltySVD = "rankCon"))
  hatc = coef(fit.cv.mrrr$fit)[-1,]
  
  # ell_32 1.69   p_120 1.6993
  h = 1/(p*l)^1.73   # for n =500, 1.83  n 100 1.832
  ystar = diag(p)*tau^2 ; B_mala = matrix(data=0,nr=p,nc=l)
  a = 0  ; M =  hatc
  for(s in 1:Iters){
    MtM = tcrossprod(M);  tam1 = solve(ystar + MtM, M );  XM = eigenMapMatMult(X , M)
    hxtyxm = h* (tXY - tX %*% exp(XM) ) * alpha
    tam = M + hxtyxm - h*(p+l+2)*tam1 + sqrt(2*h)*matrix(rnorm(p*l),nr = p)
    tattam = tcrossprod(tam);  tam2 = solve(ystar +tattam , tam );  Xtam = eigenMapMatMult(X, tam)
    pro.tam =  sum( Y*Xtam - exp(Xtam) )* alpha - 0.5*(p+l+2)*determinant(ystar + tattam)[[1]][1]
    pro.M =  sum( Y*XM - exp(XM) )* alpha - 0.5*(p+l+2)*determinant(ystar + MtM)[[1]][1]
    tran.m = -sum((M-tam - h*(tXY - tX%*%exp(Xtam) )* alpha + h*(p+l+2)*tam2 )^2)/(4*h)
    tran.tam = -sum((tam - M - hxtyxm + h*(p+l+2)*tam1 )^2)/(4*h)
    pro.trans = pro.tam+tran.m-pro.M-tran.tam
    if(log(runif(1)) <= pro.trans){  M = tam; a = a+1  } 
    if (s>burnin)B_mala = B_mala+ M/(Iters-burnin)
  } 
  print(ss);print(acceptrate[ss] <- a/Iters)  
  while (a/Iters <0.2 | a/Iters >0.6 ) {
    if(a/Iters <0.2){h = h*runif(1,max = .99,min = .8) }
    if(a/Iters >0.6){h = h*runif(1,max = 1.2,min = 1.01) }
    ystar = diag(p)*tau^2 ; B_mala = matrix(data=0,nr=p,nc=l)
    a = 0  ; M =  hatc
    for(s in 1:Iters){
      MtM = tcrossprod(M);  tam1 = solve(ystar + MtM, M );  XM = eigenMapMatMult(X , M)
      hxtyxm = h* (tXY - tX %*% exp(XM) )* alpha
      tam = M + hxtyxm - h*(p+l+2)*tam1 + sqrt(2*h)*matrix(rnorm(p*l),nr = p)
      tattam = tcrossprod(tam);  tam2 = solve(ystar +tattam , tam );  Xtam = eigenMapMatMult(X, tam)
      pro.tam =  sum( Y*Xtam - exp(Xtam) )* alpha - 0.5*(p+l+2)*determinant(ystar + tattam)[[1]][1]
      pro.M =  sum( Y*XM - exp(XM) )* alpha - 0.5*(p+l+2)*determinant(ystar + MtM)[[1]][1]
      tran.m = -sum((M-tam - h*(tXY - tX%*%exp(Xtam) )* alpha + h*(p+l+2)*tam2 )^2)/(4*h)
      tran.tam = -sum((tam - M - hxtyxm + h*(p+l+2)*tam1 )^2)/(4*h)
      pro.trans = pro.tam+tran.m-pro.M-tran.tam
      if(log(runif(1)) <= pro.trans){  M = tam; a = a+1  } 
      if (s>burnin)B_mala = B_mala+ M/(Iters-burnin)
    } 
    print(ss);print(acceptrate[ss] <- a/Iters) 
  }
  ### LMC
  M1 = hatc ; h_lmc =  h/n
  B_lmc = matrix(data= 0, nr=p, nc=l)
  for(s in 1:Iters){
    tam = solve(ystar + tcrossprod(M1),M1);  XM = eigenMapMatMult(X , M1)
    M1 = M1 + h_lmc * (tXY - tX %*% exp(XM) )* alpha - h_lmc *(p+l+2)*tam + sqrt(2*h_lmc)*rnorm(p*l)
    if(s>burnin) B_lmc =  B_lmc + M1/(Iters-burnin)
  }
  rrr[[ss]] = c( mean((hatc - C)^2 ),mean((X%*%hatc - X%*%C)^2 ), mean((ytest-exp(xtest%*%hatc))^2/mean(ytest^2)) )
  mala[[ss]] = c( mean((B_mala - C)^2 ),mean((X%*%B_mala - X%*%C)^2 ), mean((ytest-exp(xtest%*%B_mala))^2/mean(ytest^2)) )
  lmc[[ss]] = c( mean((B_lmc - C)^2 ),mean((X%*%B_lmc - X%*%C)^2 ) , mean((ytest-exp(xtest%*%B_lmc))^2/mean(ytest^2)) )
}

save.image(file = 'pois_n100q8p12_r7_alpha_0_1.rda')



