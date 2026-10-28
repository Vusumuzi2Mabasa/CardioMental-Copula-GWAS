# Independent expected-Hessian integration using base R Gaussian quadrature.
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))args[1] else 'v2'
library(jsonlite)
expected<-fromJSON(file.path(root,'results/information_controls.json'),simplifyVector=FALSE)
n<-80L;J<-matrix(0,n,n)
for(i in 1:(n-1)) {J[i,i+1]<-sqrt(i);J[i+1,i]<-sqrt(i)}
gh<-eigen(J,symmetric=TRUE);nodes<-gh$values;weights<-gh$vectors[1,]^2
grid<-expand.grid(age=seq_len(n),exposure=c(-.5,.5),z=seq_len(n),d=c(0,1))
age<-pmax(-2.5,pmin(2.5,nodes[grid$age]));z<-nodes[grid$z];d<-grid$d;e<-grid$exposure
mass0<-weights[grid$age]*.5*weights[grid$z]
errors<-c();marginal.errors<-c()
for(env in c('independent','constant','varying')) {
  b0<-expected$information[[env]]$binary_intercept
  eta<-b0+.25*age+.4*e
  rho<-if(env=='constant') rep(.4,length(age)) else if(env=='varying') tanh(atanh(.3)+.6*age) else rep(0,length(age))
  prob<-pnorm((2*d-1)*(eta+rho*z)/sqrt(1-rho^2));mass<-mass0*prob
  ll<-function(mu,shift) dnorm(z-mu,log=TRUE)+pnorm((2*d-1)*(eta+shift+rho*(z-mu))/sqrt(1-rho^2),log.p=TRUE)
  h<-2e-4;base<-ll(0,0)
  H11<--sum(mass*(ll(h,0)-2*base+ll(-h,0)))/h^2
  H22<--sum(mass*(ll(0,h)-2*base+ll(0,-h)))/h^2
  H12<--sum(mass*(ll(h,h)-ll(h,-h)-ll(-h,h)+ll(-h,-h)))/(4*h^2)
  got<-matrix(c(H11,H12,H12,H22),2,2)*.375
  wanted<-do.call(rbind,lapply(expected$information[[env]]$joint,unlist))
  errors<-c(errors,max(abs(got-wanted)))
  # Independently integrate the covariance of the two marginal scores and
  # finite-difference their expected sensitivity, rather than reuse Python's
  # closed-form E(rho*w) expression.
  pm<-pnorm(eta);phi<-dnorm(eta)
  sm<-cbind(z,ifelse(d==1,phi/pm,-phi/(1-pm)))
  Vm<-crossprod(sm,sm*mass)*.375
  marginal.ll<-function(shift) pnorm((2*d-1)*(eta+shift),log.p=TRUE)
  Hb<--sum(mass*(marginal.ll(h)-2*marginal.ll(0)+marginal.ll(-h)))/h^2
  sensitivity<-.375*diag(c(1,Hb))
  Im<-t(sensitivity)%*%solve(Vm,sensitivity)
  wanted.m<-do.call(rbind,lapply(expected$information[[env]]$marginal,unlist))
  marginal.errors<-c(marginal.errors,max(abs(Im-wanted.m)))
}
stopifnot(max(errors)<2e-6,max(marginal.errors)<2e-6)
result<-list(passed=TRUE,R=as.character(getRversion()),quadrature_nodes=80,expected_Hessian_max_error=max(errors),
  marginal_Godambe_max_error=max(marginal.errors),
  approach='Independent Jacobi eigenvalue quadrature, finite-difference likelihood Hessian, marginal score covariance and sensitivity')
write_json(result,file.path(root,'results/information_R_validation.json'),auto_unbox=TRUE,pretty=TRUE)
print(result)
