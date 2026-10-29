# Established multi-causal SuSiE fine-mapping and coloc, using signed aligned LD.
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args)) args[1] else 'v4'
phase<-if(length(args)>1) args[2] else 'cal';smoke<-length(args)>2&&args[3]=='smoke'
.libPaths(c('v3/work/Rlib',.libPaths()))
library(jsonlite);library(coloc);library(susieR)
lock<-fromJSON(file.path(root,'protocol','source_lock.json'))
for(p in names(lock$md5)) stopifnot(unname(tools::md5sum(p))==lock$md5[[p]])
input<-file.path(root,'work',if(smoke) 'fm_smoke' else 'finemap_inputs')
out<-file.path(root,if(smoke) 'work/fm_smoke_results' else 'results/finemap_chunks')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
files<-list.files(input,pattern=paste0('^',phase,'_.*json$'),full.names=TRUE)
stopifnot(length(files)>0)
for(file in files){
  stem<-sub('\\.json$','',basename(file));dest<-file.path(out,paste0(stem,'.csv'))
  metapath<-sub('csv$','json',dest)
  if(file.exists(dest)&&file.exists(metapath)){
    old<-fromJSON(metapath);stopifnot(old$identity==lock$identity,old$input_md5==unname(tools::md5sum(file)),old$output_md5==unname(tools::md5sum(dest)));next
  }
  z<-fromJSON(file,simplifyVector=FALSE);rr<-list();piprows<-list();csrows<-list();ridx<-0
  for(rec in z$records){
    m<-length(rec$c_beta);snps<-paste0('s',0:(m-1));asmat<-function(a) do.call(rbind,lapply(a,unlist))
    for(condition in c('matched_reference','in_sample')){
      Rc<-if(condition=='matched_reference') asmat(z$Rref) else asmat(rec$Rc)
      Rb<-if(condition=='matched_reference') asmat(z$Rref) else asmat(rec$Rb)
      dimnames(Rc)<-dimnames(Rb)<-list(snps,snps)
      row<-list(population=z$population,region=z$region,regime=z$regime,phase=z$phase,scenario=z$scenario,replicate=rec$replicate,LD=condition,status='ok',max_pp4=0,signal_pairs=0,c_sets=0,b_sets=0,c_covered=0,b_covered=0,c_false_sets=0,b_false_sets=0,c_max_cs=0,b_max_cs=0)
      ctrue<-as.integer(unlist(rec$true_c))+1L;btrue<-as.integer(unlist(rec$true_b))+1L
      ans<-tryCatch({
        d1<-list(beta=unlist(rec$c_beta),varbeta=unlist(rec$c_var),snp=snps,type='quant',sdY=1,N=rec$Nc,LD=Rc)
        d2<-list(beta=unlist(rec$b_beta),varbeta=unlist(rec$b_var),snp=snps,type='cc',N=rec$Nb,LD=Rb)
        cfit<-suppressWarnings(coloc::runsusie(d1,L=5,coverage=.95,maxit=1000,tol=1e-3,repeat_until_convergence=FALSE))
        bfit<-suppressWarnings(coloc::runsusie(d2,L=5,coverage=.95,maxit=1000,tol=1e-3,repeat_until_convergence=FALSE))
        stopifnot(cfit$converged,bfit$converged)
        for(trait in c('c','b')){
          fit<-if(trait=='c') cfit else bfit;truth<-if(trait=='c') ctrue else btrue
          cs<-fit$sets$cs
          row[[paste0(trait,'_sets')]]<-length(cs)
          covered<-if(length(cs)) unique(unlist(cs)) else integer()
          row[[paste0(trait,'_covered')]]<-sum(truth%in%covered)
          row[[paste0(trait,'_false_sets')]]<-sum(vapply(cs,function(s) !any(s%in%truth),logical(1)))
          row[[paste0(trait,'_max_cs')]]<-if(length(cs)) max(lengths(cs)) else 0
          common<-data.frame(population=z$population,region=z$region,regime=z$regime,phase=z$phase,scenario=z$scenario,replicate=rec$replicate,LD=condition,trait=trait)
          piprows[[length(piprows)+1]]<-cbind(common[rep(1,m),],marker=0:(m-1),pip=fit$pip,true_causal=as.integer(seq_len(m)%in%truth))
          if(length(cs)) for(k in seq_along(cs)) csrows[[length(csrows)+1]]<-cbind(common,set=names(cs)[k],markers=paste(cs[[k]]-1L,collapse=';'),size=length(cs[[k]]),covers_truth=as.integer(any(cs[[k]]%in%truth)))
        }
        if(length(cfit$sets$cs)&&length(bfit$sets$cs)){
          cc<-suppressWarnings(coloc::coloc.susie(cfit,bfit,p1=1e-4,p2=1e-4,p12=1e-5))
          if(is.data.frame(cc$summary)&&nrow(cc$summary)>0){row$signal_pairs<-nrow(cc$summary);row$max_pp4<-max(cc$summary$PP.H4.abf)}
        }
        row
      },error=function(e){row$status<-paste0('failure: ',conditionMessage(e));row})
      ridx<-ridx+1;rr[[ridx]]<-as.data.frame(ans,stringsAsFactors=FALSE)
    }
  }
  if(length(piprows)){
    con<-gzfile(file.path(out,paste0(stem,'_pip.csv.gz')),'wt');write.csv(do.call(rbind,piprows),con,row.names=FALSE);close(con)
  }
  if(length(csrows)) write.csv(do.call(rbind,csrows),file.path(out,paste0(stem,'_sets.csv')),row.names=FALSE)
  tmp<-paste0(dest,'.tmp');write.csv(do.call(rbind,rr),tmp,row.names=FALSE);stopifnot(file.rename(tmp,dest))
  write_json(list(identity=lock$identity,input_md5=unname(tools::md5sum(file)),output_md5=unname(tools::md5sum(dest)),rows=length(rr)),metapath,pretty=TRUE,auto_unbox=TRUE)
  cat(stem,'records',length(rr),'\n')
}
write_json(list(R=R.version.string,coloc=as.character(packageVersion('coloc')),susieR=as.character(packageVersion('susieR')),L=5,coverage=.95,LD_shrinkage=.05),file.path(out,'software.json'),pretty=TRUE,auto_unbox=TRUE)
