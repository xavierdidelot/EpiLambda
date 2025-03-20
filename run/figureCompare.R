library(EpiLambda)
library(ggplot2)
library(viridis)
rm(list=ls())
n=10

alphas=c(0.5,1,1.5)
data=data.frame()
for (i in 1:12) {
  if (i<=3) {
    alpha=alphas[i]
    d=beta_psize(n=n,alpha=alpha)
    nam=paste0('Beta(alpha=',as.character(alpha),')')
  }
  else
  {
    if (i<=6) nt=10
    if (i>6 && i<=9) nt=20
    if (i>9) nt=30
    r=c(10,0.1,1)[i%%3+1]
    d=omega_psize(n=n,nt=nt,r=r)
    nam=paste0('Omega(N=',as.character(nt),',r=',as.character(r),')')
  }
  data=rbind(data.frame(k=2:n,alpha=nam,i=i,p=d[1:(n-1)]),data)
}
limits=c(-4,0)
data$Probability=pmin(0,pmax(-4,log10(data$p)))
data$alpha=factor(data$alpha,levels=unique(data$alpha))#forces order to remain as input
pdf('figureCompare.pdf')
#ggplot(data, aes(k,p,fill=alpha))+geom_bar(stat="identity",position='dodge')+ scale_x_continuous(breaks=seq(2,n,1))+xlab('Size of next event')+ylab('Probability of event')
ggplot(data, aes(k, alpha, fill=Probability)) + geom_tile()+scale_fill_viridis(limits=limits,breaks=limits[1]:limits[2],labels=c('<0.0001','0.001','0.01','0.1','1'))+
  scale_x_continuous(breaks=seq(2,n,1),expand = c(0, 0))+scale_y_discrete(expand = c(0, 0))+xlab('Size of next event')+ylab('')+
  geom_segment(aes(x = -1, xend = 10.5, y = 3.5, yend = 3.5), linetype="dashed", color = "red")+
  geom_segment(aes(x = -1, xend = 10.5, y = 6.5, yend = 6.5), linetype="dashed", color = "red")+
  geom_segment(aes(x = -1, xend = 10.5, y = 9.5, yend = 9.5), linetype="dashed", color = "red")+
  coord_cartesian(xlim = c(1.5, 10.5), clip="off")
dev.off()
system('open figureCompare.pdf')
