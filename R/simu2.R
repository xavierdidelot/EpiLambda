simul_pcoal2=function(rt,vt=rt,nt,nt1=round(rt*nt),n=nt1,type='inclusive',nrep=1e3,method='rejection')
{
  p=rep(0,n)
  for (i in 1:nrep) {
    if (method=='rejection') {
      ks=0
      while (sum(ks)!=n){
        if (vt==rt) ks=rpois(nt,rt)
        else ks=rnbinom(nt,rt*rt/(vt-rt),rt/vt)
      }
    } else if (method=='multinomial') {
      if (vt==rt) ks=rmultinom(1, n, rep(1/nt,nt))
      else ks=rmultinom(1, n, rgamma(nt,shape=rt*rt/(vt-rt),scale=1))#scale doesn't actually matter
    }
    if (type=='inclusive') p[1]=p[1]+1
    else p[1]=p[1]+sum(ks==1)/n
    for (j in 2:n) {
      v=0:(j-1)
      p0=0
      for (k in ks) p0=p0+prod(k-v)*ifelse(type=='inclusive',k>=j,k==j)
      p[j]=p[j]+p0/prod(n-v)
    }
  }
  p=p/nrep
  return(p)
}

compare_simu2=function(rt=1.5,vt=10,nt=10,nt1=20,n=15,type='inclusive',nrep=1e3,method='multinomial')
{
  s=simul_pcoal(rt=rt,vt=vt,nt=nt,nt1=nt1,n=n,type=type,nrep=nrep,method=method)
  s2=simul_pcoal2(rt=rt,vt=vt,nt=nt,nt1=nt1,n=n,type=type,nrep=nrep,method=method)
  if (type=='inclusive') {
    if (rt==vt) p=pois_inclusive(1:nt1,nt)
    else {
      r=rt^2/(vt-rt)
      p=sapply(1:nt1,function(k) negbin_inclusive(k,nt,r))
    }
  } else {
    if (rt==vt) p=pois_exclusive(1:n,n,nt)
    else {
      r=rt^2/(vt-rt)
      p=sapply(1:nt1,function(k) negbin_exclusive(k,n,nt,r))
    }
  }

  ret=cbind(s,s2,p,s-p, s2-p)
  colnames(ret)=c('Simulated','Simulated_2', 'Calculated','Difference', 'Difference_2')
  return(ret)
}
