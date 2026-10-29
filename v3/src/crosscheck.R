# Independent R lm/glm and the established coloc package; never source personal startup files.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else 'v3'
local <- file.path(root,'work','Rlib')
if(dir.exists(local)) .libPaths(c(normalizePath(local),.libPaths()))
stopifnot(requireNamespace('jsonlite',quietly=TRUE),requireNamespace('coloc',quietly=TRUE))
files <- list.files(file.path(root,'work','r_fixtures'),pattern='json$',full.names=TRUE)
stopifnot(length(files)>0)
errs <- c(ols=0,logistic=0,ols_robust=0,logistic_robust=0,posterior=0)
fit_count <- 0
for(file in files){
  z <- jsonlite::fromJSON(file);X<-z$X;G<-z$G;y<-z$y;d<-z$d
  for(j in seq_len(ncol(G))){
    A<-cbind(X,G[,j]);p<-ncol(A)
    ol<-lm.fit(A,y);iv<-solve(crossprod(A));res<-ol$residuals
    ol_v<-sum(res^2)/(length(y)-p)*iv[p,p]
    ol_i<-as.vector(A%*%iv[,p])*res
    gl<-glm.fit(A,d,family=binomial(),control=glm.control(epsilon=1e-10,maxit=100))
    stopifnot(gl$converged)
    w<-gl$fitted.values*(1-gl$fitted.values);vi<-solve(crossprod(A,A*w))
    gl_i<-as.vector(A%*%vi[,p])*(d-gl$fitted.values)
    errs['ols']<-max(errs['ols'],abs(ol$coefficients[p]-z$c_beta[j]),abs(ol_v-z$c_var[j]))
    errs['logistic']<-max(errs['logistic'],abs(gl$coefficients[p]-z$b_beta[j]),abs(vi[p,p]-z$b_var[j]))
    errs['ols_robust']<-max(errs['ols_robust'],abs(sum(ol_i^2)-z$c_robust[j]))
    errs['logistic_robust']<-max(errs['logistic_robust'],abs(sum(gl_i^2)-z$b_robust[j]))
    fit_count<-fit_count+2
  }
  ids<-paste0('s',seq_len(ncol(G)))
  cc<-suppressWarnings(coloc::coloc.abf(
    list(beta=z$split_c_beta,varbeta=z$split_c_var,snp=ids,type='quant',sdY=1),
    list(beta=z$split_b_beta,varbeta=z$split_b_var,snp=ids,type='cc'),
    p1=1e-4,p2=1e-4,p12=1e-5))
  pp<-as.numeric(cc$summary[paste0('PP.H',0:4,'.abf')])
  errs['posterior']<-max(errs['posterior'],max(abs(pp-z$pp)))
}
stopifnot(max(errs)<2e-6)
out<-list(passed=TRUE,fixtures=length(files),independent_fits=fit_count,max_absolute_errors=as.list(errs),R=R.version.string,coloc=as.character(packageVersion('coloc')),jsonlite=as.character(packageVersion('jsonlite')))
jsonlite::write_json(out,file.path(root,'results','tests','R_validation.json'),pretty=TRUE,auto_unbox=TRUE)
print(out)
