# Established-method controls for sharing, opposing signs and signed LD harmonisation.
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args)) args[1] else 'v4'
.libPaths(c('v3/work/Rlib',.libPaths()))
library(coloc);library(susieR);library(jsonlite)
m<-40;N<-10000;ids<-paste0('s',seq_len(m));R<-.6^abs(outer(seq_len(m),seq_len(m),'-'));dimnames(R)<-list(ids,ids)
zc<-8*R[,7];zs<-9*R[,7];zd<-9*R[,32]
run<-function(z1,z2,rr){
  a<-list(beta=z1/sqrt(N),varbeta=rep(1/N,m),type='quant',sdY=1,N=N,snp=ids,LD=rr)
  b<-list(beta=z2/sqrt(N),varbeta=rep(1/N,m),type='cc',N=N,snp=ids,LD=rr)
  s1<-suppressWarnings(runsusie(a,L=5,coverage=.95,maxit=1000,repeat_until_convergence=FALSE))
  s2<-suppressWarnings(runsusie(b,L=5,coverage=.95,maxit=1000,repeat_until_convergence=FALSE))
  stopifnot(s1$converged,s2$converged,length(s1$sets$cs)>0,length(s2$sets$cs)>0)
  cc<-coloc.susie(s1,s2,p1=1e-4,p2=1e-4,p12=1e-5)
  list(pp=max(cc$summary$PP.H4.abf),pip1=s1$pip,pip2=s2$pip)
}
shared<-run(zc,zs,R);opposing<-run(zc,-zs,R);distinct<-run(zc,zd,R)
signs<-rep(1,m);signs[7]<--1;Rf<-R*outer(signs,signs)
flipped<-run(zc*signs,zs*signs,Rf)
err<-max(abs(shared$pip1-flipped$pip1),abs(shared$pip2-flipped$pip2),abs(shared$pp-flipped$pp),abs(shared$pp-opposing$pp))
stopifnot(shared$pp>.95,opposing$pp>.95,distinct$pp<.05,err<1e-7)
out<-list(passed=TRUE,shared_PP4=shared$pp,opposing_PP4=opposing$pp,distinct_PP4=distinct$pp,max_harmonisation_error=err,coloc=as.character(packageVersion('coloc')),susieR=as.character(packageVersion('susieR')))
dir.create(file.path(root,'results','tests'),recursive=TRUE,showWarnings=FALSE)
write_json(out,file.path(root,'results','tests','R_finemap_validation.json'),pretty=TRUE,auto_unbox=TRUE)
print(out)
