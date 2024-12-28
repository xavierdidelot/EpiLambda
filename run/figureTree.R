library(EpiLambda)
library(ggplot2)
library(ggtree)
library(patchwork)
set.seed(0)
pdf('figureTree.pdf',10,10)
ps=list()
rs=c(0.1,1,5,10)
for (i in 1:4) {
  t=omega_simtree(n=20,nt=100,r=rs[i])
  p <- ggtree(t) + theme_tree2() + scale_x_continuous(labels = abs)+ geom_tiplab()+xlab("Generations since present")
  p = revts(p)
  ps[[i]]=p
}
ps[[1]]+ps[[2]]+ps[[3]]+ps[[4]]+plot_layout(ncol = 2, nrow = 2, guides = "collect")+plot_annotation(tag_levels = 'A')
dev.off()
system('open figureTree.pdf')
