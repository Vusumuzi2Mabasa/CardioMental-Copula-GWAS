# Independent R optimisation of bivariate restricted likelihood, including binary working fits.
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args)) args[1] else 'v4'
library(jsonlite)
files<-file.path(root,'results','tests','reml_fixture.json')
more<-list.files(file.path(root,'results','polygenic_fixtures'),pattern='g0.0.*json$',full.names=TRUE)
more<-c(more,list.files(file.path(root,'results','polygenic_validation','polygenic_fixtures'),pattern='g0.0.*json$',full.names=TRUE))
files<-c(files,more);stopifnot(length(files)>0)
maxll<-0;maxrg<-0
for(file in files){
  z<-fromJSON(file);l<-z$eigenvalues;y<-z$y;nn<-length(l)
  fn<-function(a,null=FALSE){
    if(null) a<-c(a[1:2],0,a[3:5])
    gc<-exp(2*a[1]);gm<-exp(2*a[2]);gcm<-exp(a[1]+a[2])*tanh(a[3])
    ec<-exp(2*a[4]);em<-exp(2*a[5]);ecm<-exp(a[4]+a[5])*tanh(a[6])
    v1<-l*gc+ec;v2<-l*gm+em;cross<-l*gcm+ecm;det<-v1*v2-cross^2
    if(any(!is.finite(det))||any(det<=0)) return(1e50)
    sum(log(det)+(v2*y[,1]^2-2*cross*y[,1]*y[,2]+v1*y[,2]^2)/det)/2
  }
  total<-colMeans(y*y);initial<-c(log(sqrt(.4*total)),0,log(sqrt(.6*total)),0)
  fit<-optim(initial,fn,method='BFGS',control=list(maxit=2000,reltol=1e-12))
  null<-optim(initial[c(1,2,4,5,6)],function(a) fn(a,TRUE),method='BFGS',control=list(maxit=2000,reltol=1e-12))
  stopifnot(fit$convergence==0,null$convergence==0)
  llerr<-max(abs(fit$value+z$expected$loglik_alt),abs(null$value+z$expected$loglik_null))
  rgerr<-abs(tanh(fit$par[3])-z$expected$genetic_correlation)
  maxll<-max(maxll,llerr);maxrg<-max(maxrg,rgerr)
  stopifnot(llerr<1e-4,rgerr<5e-4)
}
out<-list(passed=TRUE,fixtures=length(files),max_loglikelihood_error=maxll,max_genetic_correlation_error=maxrg,R=R.version.string)
write_json(out,file.path(root,'results','tests','R_reml_validation.json'),pretty=TRUE,auto_unbox=TRUE)
print(out)
