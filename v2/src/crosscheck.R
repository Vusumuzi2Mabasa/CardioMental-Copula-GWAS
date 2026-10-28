args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else 'v2'
library(jsonlite)
dat <- read.csv(file.path(root,'results/oracle_fixture.csv'))
expected <- fromJSON(file.path(root,'results/oracle_fixture_expected.json'),simplifyVector=FALSE)
x <- cbind(1,dat$age,dat$exposure)
w <- cbind(1,dat$age)
alpha <- c(atanh(.3),.6)
ll <- function(v, allknown, fixed) {
  theta <- if(allknown) fixed else c(v[3:9],alpha)
  gg <- dat$g-if(allknown) .5 else mean(dat$g)
  mu <- drop(x%*%theta[1:3])+v[1]*gg
  eta <- drop(x%*%theta[4:6])+v[2]*gg
  z <- (dat$y-mu)/exp(theta[7]);rho <- tanh(drop(w%*%alpha))
  dnorm(z,log=TRUE)-theta[7]+pnorm((2*dat$d-1)*(eta+rho*z)/sqrt(1-rho^2),log.p=TRUE)
}
errors <- c();fit.error <- 0
for(k in seq_along(expected)) {
  ref <- expected[[k]];allknown <- ref$all_known;theta <- unlist(ref$theta)
  if(!allknown) {
    start <- c(lm(dat$y~dat$age+dat$exposure)$coefficients,qnorm(mean(dat$d)),0,0,0)
    fit <- optim(start,function(nu)-sum(ll(c(0,0,nu),FALSE,theta)),method='BFGS',control=list(reltol=1e-12,maxit=2000))
    stopifnot(fit$convergence==0)
    fit.error <- max(abs(fit$par-theta[1:7]));stopifnot(fit.error<2e-4)
  }
  v <- if(allknown) c(0,0) else c(0,0,theta[1:7])
  h <- 1e-4*(1+abs(v));ss <- matrix(0,nrow(dat),length(v))
  for(j in seq_along(v)) {a<-b<-v;a[j]<-a[j]+h[j];b[j]<-b[j]-h[j];ss[,j]<-(ll(a,allknown,theta)-ll(b,allknown,theta))/(2*h[j])}
  H <- optimHess(v,function(a)-sum(ll(a,allknown,theta)),control=list(ndeps=rep(1e-4,length(v))))
  u <- colSums(ss[,1:2]);info <- H[1:2,1:2];psi<-ss[,1:2]
  if(!allknown) {a<-H[1:2,-(1:2)]%*%solve(H[-(1:2),-(1:2)]);u<-u-drop(a%*%colSums(ss[,-(1:2)]));info<-info-a%*%H[-(1:2),1:2];psi<-psi-ss[,-(1:2)]%*%t(a)}
  tm<-drop(t(u)%*%solve(info,u));tr<-drop(t(u)%*%solve(crossprod(psi),u))
  errors<-c(errors,abs(pchisq(tm,2,lower.tail=FALSE)-ref$model_p),abs(pchisq(tr,2,lower.tail=FALSE)-ref$robust_p))
}
stopifnot(max(errors)<2e-5)
result<-list(passed=TRUE,R=as.character(getRversion()),max_score_p_error=max(errors),restricted_fit_max_parameter_error=fit.error,approach='Independent base R likelihood finite differences BFGS and expanded Hessian')
write_json(result,file.path(root,'results/R_validation.json'),auto_unbox=TRUE,pretty=TRUE)
print(result)
