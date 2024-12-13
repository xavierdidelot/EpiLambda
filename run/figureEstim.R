rm(list=ls())
library(EpiLambda)
library(ggplot2)
library(patchwork)
set.seed(1)


data=data.frame()
for (i in 1:100) {
  print(i)
  nt=runif(1,100,500)
  r=runif(1,0.01,2)
  t=new_simtree(100,nt=nt,r=r)
  p1=new_mle(t,r=r,ntbounds=c(100,500),rbounds=c(0,2))
  p2=new_mle(t,nt=nt,ntbounds=c(100,500),rbounds=c(0,2))
  p=new_mle(t,ntbounds=c(100,500),rbounds=c(0,2))
  data=rbind(data,data.frame(nt=nt,r=r,ont=p1,or=p2,ent=p[1],er=p[2]))
}

p1<-ggplot(data,aes(x=nt,y=ont))+geom_point()
p2<-ggplot(data,aes(x=r,y=or))+geom_point()
p3<-ggplot(data,aes(x=nt,y=ent))+geom_point()
p4<-ggplot(data,aes(x=r,y=er))+geom_point()
pdf('figureEstim.pdf')
p1+p2+p3+p4+plot_layout(ncol = 2, nrow = 2, guides = "collect")+plot_annotation(tag_levels = 'A')
dev.off()
system('open figureEstim.pdf')

set.seed(3)
nt=200
r=0.5
t=new_simtree(100,nt=nt,r=r)
ntbounds=c(100,500)
rbounds=c(0.05,2)
ss=seq(0,1,0.05)
dat=data.frame()
for (nts in ss*(ntbounds[2]-ntbounds[1])+ntbounds[1]) for (rs in ss*(rbounds[2]-rbounds[1])+rbounds[1]) {
  v=new_loglik(t,nts,rs)
  dat=rbind(dat,data.frame(x=rs,y=nts,z=v))
}
ggplot(dat,aes(x = x, y = y, z = z)) +xlab('Dispersion parameter')+ylab('Population size')+
   scale_x_continuous(limits=c(0,2), expand = c(0, 0))+ scale_y_continuous(limits=c(100,500), expand = c(0, 0))+geom_contour_filled()
