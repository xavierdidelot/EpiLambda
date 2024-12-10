rm(list=ls())
library(EpiLambda)
library(ggplot2)
library(patchwork)

rt=1
data=data.frame()
v=sapply(1:10, function(k) pois_inclusive(k=k,nt=10))
data=rbind(data,data.frame(x=1:10,y=v,PopSize='10',Type='Inclusive'))
v=sapply(1:10, function(k) pois_exclusive(k=k,n=10,nt=10))
data=rbind(data,data.frame(x=1:10,y=v,PopSize='10',Type='Exclusive'))
v=sapply(1:10, function(k) pois_inclusive(k=k,nt=20))
data=rbind(data,data.frame(x=1:10,y=v,PopSize='20',Type='Inclusive'))
v=sapply(1:10, function(k) pois_exclusive(k=k,n=10,nt=20))
data=rbind(data,data.frame(x=1:10,y=v,PopSize='20',Type='Exclusive'))
v=sapply(1:10, function(k) pois_inclusive(k=k,nt=30))
data=rbind(data,data.frame(x=1:10,y=v,PopSize='30',Type='Inclusive'))
v=sapply(1:10, function(k) pois_exclusive(k=k,n=10,nt=30))
data=rbind(data,data.frame(x=1:10,y=v,PopSize='30',Type='Exclusive'))
data$Type=factor(data$Type,levels=unique(data$Type))#forces order to remain as input
data$PopSize=factor(data$PopSize,levels=unique(data$PopSize))#forces order to remain as input
plot1<-ggplot(data, aes(x=x)) + scale_x_continuous(breaks=seq(0,10,1)) +theme(panel.grid.minor.x = element_blank() )+
  geom_line(aes(y = y, colour = PopSize,group=interaction(PopSize,Type))) + geom_point(aes(y=y,colour=PopSize,shape=Type))+
  xlab('Size of multimerger event')+ylab('Probability') + scale_y_log10(breaks=10^-(0:15),limits=c(1e-15,1))+theme(panel.grid.minor.y = element_blank() )


pdf('figurePois.pdf',7,7)
plot1
dev.off()
system('open figurePois.pdf')
