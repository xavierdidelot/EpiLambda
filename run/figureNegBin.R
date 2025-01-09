rm(list=ls())
library(EpiLambda)
library(ggplot2)
library(patchwork)

data=data.frame()
rt=2
data=rbind(data,data.frame(x=0:10,y=dpois(0:10,rt),Distribution='Poisson'))
rs=c(0.1,1,10,100)
for (i in 1:length(rs))
  data=rbind(data,data.frame(x=0:10,y=dnbinom(0:10,size=rs[i],mu=rt),Distribution=paste0('NegBin(r=',rs[i],')')))
data$Distribution=factor(data$Distribution,levels=unique(data$Distribution))#forces order to remain as input
plot1<-ggplot(data, aes(x=x)) + scale_x_continuous(breaks=seq(0,10,1)) +theme(panel.grid.minor.x = element_blank() )+
  geom_line(aes(y = y, colour = Distribution)) + geom_point(aes(y=y,colour=Distribution))+
  xlab('Offspring number')+ylab('Probability')


nt=20
n=10
data=data.frame()
v=sapply(1:10, function(k) exp(pois_inclusive(k=k,nt=nt,log=T)))
data=rbind(data,data.frame(x=1:10,y=v,Distribution='Poisson'))
for (i in 1:length(rs)) {
  v=sapply(1:10, function(k) exp(negbin_inclusive(k=k,nt=nt,r=rs[i],log=T)))
  data=rbind(data,data.frame(x=1:10,y=v,Distribution=paste0('NegBin(r=',rs[i],')')))}
data$Distribution=factor(data$Distribution,levels=unique(data$Distribution))#forces order to remain as input
plot2<-ggplot(data, aes(x=x)) + scale_x_continuous(breaks=seq(0,10,1)) +theme(panel.grid.minor.x = element_blank() )+
  geom_line(aes(y = y, colour = Distribution)) + geom_point(aes(y=y,colour=Distribution))+
  xlab('Size of multimerger event')+ylab('Inclusive Probability') + scale_y_log10(breaks=10^-(0:12),limits=c(1e-12,1))+theme(panel.grid.minor.y = element_blank() )

data=data.frame()
v=sapply(1:10, function(k) exp(pois_exclusive(k=k,n=n,nt=nt,log=T)))
data=rbind(data,data.frame(x=1:10,y=v,Distribution='Poisson'))
for (i in 1:length(rs)) {
  v=sapply(1:10, function(k) exp(negbin_exclusive(k=k,n=n,nt=nt,r=rs[i],log=T)))
  data=rbind(data,data.frame(x=1:10,y=v,Distribution=paste0('NegBin(r=',rs[i],')')))}
data$Distribution=factor(data$Distribution,levels=unique(data$Distribution))#forces order to remain as input
plot3<-ggplot(data, aes(x=x)) + scale_x_continuous(breaks=seq(0,10,1)) +theme(panel.grid.minor.x = element_blank() )+
  geom_line(aes(y = y, colour = Distribution)) + geom_point(aes(y=y,colour=Distribution))+
  xlab('Size of multimerger event')+ylab('Exclusive Probability') + scale_y_log10(breaks=10^-(0:12),limits=c(1e-12,1))+theme(panel.grid.minor.y = element_blank() )

pdf('figureNegBin.pdf',7,7)
plot1+plot2+plot3+guide_area()+plot_layout(ncol = 2, nrow = 2, guides = "collect")+plot_annotation(tag_levels = 'A')
dev.off()
system('open figureNegBin.pdf')
