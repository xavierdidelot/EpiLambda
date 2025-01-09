rm(list=ls())
library(EpiLambda)
library(ggplot2)
library(patchwork)
set.seed(1)

ntbounds=c(100,500)
rbounds=c(0.01,2)

data=data.frame()
for (i in 1:100) {
  print(i)
  nt=runif(1,ntbounds[1],ntbounds[2])
  r=runif(1,rbounds[1],rbounds[2])
  t=omega_simtree(100,nt=nt,r=r)
  p1=omega_mle(t,r=r,ntbounds=ntbounds,rbounds=c(0,2))
  p2=omega_mle(t,nt=nt,ntbounds=ntbounds,rbounds=c(0,2))
  p=omega_mle(t,ntbounds=ntbounds,rbounds=c(0,2))
  data=rbind(data,data.frame(nt=nt,r=r,ont=p1,or=p2,ent=p[1],er=p[2]))
}

p1<-ggplot(data,aes(x=nt,y=ont))+geom_point()+xlab('True population size')+ylab('Estimated population size')+geom_smooth( )+scale_x_continuous(limits=c(100,500))+ scale_y_continuous(limits=c(100,500))
p2<-ggplot(data,aes(x=r,y=or))+geom_point()+xlab('True dispersion')+ylab('Estimated dispersion')+geom_smooth( )+scale_x_continuous(limits=c(0,2))+ scale_y_continuous(limits=c(0,2))
p3<-ggplot(data,aes(x=nt,y=ent))+geom_point()+xlab('True population size')+ylab('Estimated population size')+geom_smooth( )+scale_x_continuous(limits=c(100,500))+ scale_y_continuous(limits=c(100,500))
p4<-ggplot(data,aes(x=r,y=er))+geom_point()+xlab('True dispersion')+ylab('Estimated dispersion')+geom_smooth( )+scale_x_continuous(limits=c(0,2))+ scale_y_continuous(limits=c(0,2))

set.seed(3)
nt=200
r=0.5
t=omega_simtree(100,nt=nt,r=r)
rbounds=c(0.05,2)

ss=seq(0,1,0.05)
dat=data.frame()
for (nts in ss*(ntbounds[2]-ntbounds[1])+ntbounds[1]) for (rs in ss*(rbounds[2]-rbounds[1])+rbounds[1]) {
  v=omega_loglik(t,nts,rs)
  dat=rbind(dat,data.frame(x=rs,y=nts,z=v))
}
p5<-ggplot(dat,aes(x = x, y = y, z = z)) +labs(x='Dispersion parameter',y='Population size',fill='Likelihood')+
  scale_x_continuous(limits=c(0,rbounds[2]), expand = c(0, 0))+ scale_y_continuous(limits=ntbounds, expand = c(0, 0))+geom_contour_filled()


pdf('figureEstim.pdf',8,12)
layout <- "
AB
CD
EF"
p1+p2+p3+p4+p5+ guide_area()+plot_layout(design=layout, guides = "collect")+plot_annotation(tag_levels = 'A')
dev.off()
system('open figureEstim.pdf')

