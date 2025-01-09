library(EpiLambda)
library(ggplot2)
library(patchwork)
set.seed(0)
data=data.frame()

rs=c(0.1,1,10,100)
n=20
nt=30
for (i in 1:10000) {
  r=rs[i%%4+1]
  t=omega_simtree(n=n,nt=nt,r=r)
  tmrca=max(dist.nodes(t)[n+1,])
  bralen=sum(t$edge.length)
  stem=sum(t$edge.length[-(1:n)])/sum(t$edge.length)
  nmulti=Ntip(t)-Nnode(t)-1
  largest=max(table(t$edge[,1]))
  data=rbind(data,data.frame(tmrca=tmrca,bralen=bralen,stem=stem,nmulti=nmulti,largest=largest,r=r))
}
p1 <- ggplot(data, aes(factor(r), tmrca)) + geom_violin()+geom_boxplot(width=.1,outliers = F)+xlab('Dispersion parameter')+ylab('Time to the most recent common ancestor')
p2 <- ggplot(data, aes(factor(r), stem)) + geom_violin()+geom_boxplot(width=.1,outliers = F)+xlab('Dispersion parameter')+ylab('Stemminess')
p3 <- ggplot(data, aes(factor(r), nmulti)) + geom_violin()+geom_boxplot(width=.1,outliers = F)+xlab('Dispersion parameter')+ylab('Number of multimergers')
p4 <- ggplot(data, aes(factor(r), largest)) + geom_violin()+geom_boxplot(width=.1,outliers = F)+xlab('Dispersion parameter')+ylab('Size of the largest multimerger')
pdf('figureStats.pdf',10,10)
p3+p4+p1+p2+plot_layout(ncol = 2, nrow = 2, guides = "collect")+plot_annotation(tag_levels = 'A')
dev.off()
system('open figureStats.pdf')
