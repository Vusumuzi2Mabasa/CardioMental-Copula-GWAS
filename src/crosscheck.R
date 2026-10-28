# Independent base-R likelihood, numerical derivatives and score calculation.
# Run with Rscript --vanilla. No personal startup files or external data.
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) args[1] else "."
d <- read.csv(file.path(root,"results/synthetic_fixture.csv"))
ref <- jsonlite::fromJSON(file.path(root,"results/fixture_python.json"), simplifyVector=FALSE)
x <- cbind(1,d$age,d$exposure); g <- d$g0-mean(d$g0)
y <- d$continuous; outcome <- d$binary

ll <- function(theta,x,w) {
  p <- ncol(x)
  mu <- drop(x %*% theta[seq_len(p)])
  eta <- drop(x %*% theta[p+seq_len(p)])
  sigma <- exp(theta[2*p+1])
  rho <- if (ncol(w)) tanh(drop(w %*% theta[2*p+1+seq_len(ncol(w))])) else rep(0,length(y))
  z <- (y-mu)/sigma
  a <- (eta+rho*z)/sqrt(1-rho^2)
  dnorm(y,mu,sigma,log=TRUE)+pnorm((2*outcome-1)*a,log.p=TRUE)
}
fitR <- function(xx,w) {
  ols <- lm.fit(xx,y)
  prob <- glm.fit(xx,outcome,family=binomial(link="probit"))
  initial <- c(ols$coefficients,prob$coefficients,
               log(sqrt(mean(ols$residuals^2))),rep(0,ncol(w)))
  fn <- function(th) -sum(ll(th,xx,w))
  starts <- list(initial)
  if (ncol(w)) for (r in c(-.4,.4)) {
    z <- initial;z[2*ncol(xx)+2] <- atanh(r);starts[[length(starts)+1]] <- z
  }
  fits <- lapply(starts,function(z) optim(z,fn,method="BFGS",control=list(maxit=1500,reltol=1e-12)))
  opt <- fits[[which.min(vapply(fits,function(z) z$value,numeric(1)))]]
  stopifnot(opt$convergence==0)
  list(theta=opt$par,nll=opt$value,fn=fn)
}
scoresR <- function(theta,xx,w) {
  s <- matrix(NA_real_,length(y),length(theta))
  for (j in seq_along(theta)) {
    h <- 1e-5*(1+abs(theta[j]));a <- b <- theta
    a[j] <- a[j]+h;b[j] <- b[j]-h
    s[,j] <- (ll(a,xx,w)-ll(b,xx,w))/(2*h)
  }
  s
}
rows <- list()
for (k in 0:2) {
  w <- if (k==0) matrix(numeric(),length(y),0) else if (k==1) matrix(1,length(y),1) else cbind(1,d$age)
  f <- fitR(x,w)
  # Construct the expanded SNP-effect parameterization at beta=0.
  xx <- cbind(x,g);p <- ncol(x)
  th <- c(f$theta[seq_len(p)],0,f$theta[p+seq_len(p)],0,f$theta[-seq_len(2*p)])
  s <- scoresR(th,xx,w)
  h <- optimHess(th,function(z) -sum(ll(z,xx,w)),control=list(ndeps=rep(2e-5,length(th))))
  bi <- c(p+1,2*p+2);ei <- setdiff(seq_along(th),bi)
  a <- h[bi,ei,drop=FALSE] %*% solve(h[ei,ei,drop=FALSE])
  u <- colSums(s)[bi]-drop(a %*% colSums(s)[ei])
  info <- h[bi,bi,drop=FALSE]-a %*% h[ei,bi,drop=FALSE]
  e <- s[,bi,drop=FALSE]-s[,ei,drop=FALSE] %*% t(a)
  robust <- crossprod(e)
  pm <- pchisq(drop(t(u)%*%solve(info,u)),2,lower.tail=FALSE)
  pr <- pchisq(drop(t(u)%*%solve(robust,u)),2,lower.tail=FALSE)
  target <- ref$null_checks[[k+1]]
  err <- max(abs(c(pm-target$model_p,pr-target$robust_p)))
  stopifnot(abs(f$nll-target$nll)<1e-4,err<3e-5)
  rows[[k+1]] <- data.frame(dependence_columns=k,nll=f$nll,
      nll_error=abs(f$nll-target$nll),model_p=pm,robust_p=pr,p_error=err)
}
ff <- fitR(cbind(x,g),cbind(1,d$age))
ss <- scoresR(ff$theta,cbind(x,g),cbind(1,d$age))
hh <- optimHess(ff$theta,ff$fn,control=list(ndeps=rep(2e-5,length(ff$theta))))
hi <- solve(hh);vv <- hi %*% crossprod(ss) %*% hi
ix <- c(4,8)
effect_error <- max(abs(ff$theta[ix]-unlist(ref$effects)))
se_error <- max(abs(sqrt(diag(vv)[ix])-unlist(ref$sandwich_se)))
stopifnot(effect_error<2e-4,se_error<2e-4)
write.csv(do.call(rbind,rows),file.path(root,"results/R_score_crosscheck.csv"),row.names=FALSE)
jsonlite::write_json(list(passed=TRUE,engine=R.version.string,
    engine_md5=unname(tools::md5sum(file.path(root,"src/model.py"))),
    max_score_p_error=max(vapply(rows,function(z) z$p_error,numeric(1))),
    full_effect_max_error=effect_error,sandwich_se_max_error=se_error,
    calculation="Independent base-R likelihood, BFGS fits and numerical derivatives"),
    file.path(root,"results/R_validation.json"),pretty=TRUE,auto_unbox=TRUE)
cat("Independent R validation passed\n")
