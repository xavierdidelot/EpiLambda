rm(list=ls())
library(EpiLambda)
library(ggplot2)
library(patchwork)

data=data.frame()
rt=2
data=rbind(data,data.frame(x=0:10,y=dpois(0:10,rt),Distribution='Poisson'))
rs=c(1,2,10,100)
for (i in 1:length(rs))
  data=rbind(data,data.frame(x=0:10,y=dnbinom(0:10,size=rs[i],mu=rt),Distribution=paste0('NegBin(r=',rs[i],')')))
data$Distribution=factor(data$Distribution,levels=unique(data$Distribution))#forces order to remain as input
plot1<-ggplot(data, aes(x=x)) +
  geom_line(aes(y = y, colour = Distribution)) + geom_point(aes(y=y,colour=Distribution))+
  xlab('Offspring number')+ylab('Probability')


nt=10
data=data.frame()
v=sapply(2:10, function(n) pois_pnt(n,nt=nt))
data=rbind(data,data.frame(x=2:10,y=v,Distribution='Poisson'))
for (i in 1:length(rs)) {
  v=sapply(2:10, function(n) negbin_pnt(n,nt=nt,r=rs[i]))
  data=rbind(data,data.frame(x=2:10,y=v,Distribution=paste0('NegBin(r=',rs[i],')')))}
data$Distribution=factor(data$Distribution,levels=unique(data$Distribution))#forces order to remain as input
plot2<-ggplot(data, aes(x=x)) +
  geom_line(aes(y = y, colour = Distribution)) + geom_point(aes(y=y,colour=Distribution))+
  xlab('Size of multimerger event')+ylab('Probability') + scale_y_log10()

pdf('figure.pdf',10,5)
plot1+plot2+plot_layout(ncol = 2, nrow = 1, guides = "collect")+plot_annotation(tag_levels = 'A')
dev.off()
system('open figure.pdf')
