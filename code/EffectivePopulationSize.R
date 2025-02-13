# Effective Population Size--Phylodyn

rm(list=ls())

library(phylodyn)
library(ggtree)
library(treeio)
library(beastio)
library(phytools)
library(ggplot2)
library(scales)

# Take the first sample obtained on day 3 (S1_3d) as an example

M1D3_an <- read.csv("~/.../data/S1_3d_annotation.csv",header = T) 
M1D3_an <- data.frame(M1D3_an)

M1D3_reduced <- read.tree('~/.../data/IQ-TREE/S1_3d.nwk')
M1D3_reducedfor <- as.data.frame(fortify(M1D3_reduced))

M1D3_reducedHei <- nodeHeights(M1D3_reduced)
M1D3_reducedHei_coor <- M1D3_reduced$edge
nM1D3_reduced_1 <- c(M1D3_reducedHei[,1],M1D3_reducedHei[,2])
nM1D3_reduced_2 <- c(M1D3_reducedHei_coor[,1],M1D3_reducedHei_coor[,2])
nM1D3_reduced <- data.frame(nM1D3_reduced_1,nM1D3_reduced_2)
colnames(nM1D3_reduced) <- c('data','name')
M1D3_reducednew <- nM1D3_reduced[order(nM1D3_reduced$name),]
M1D3_reducedfinal <- unique(M1D3_reducednew)

colnames(M1D3_reducedfinal) <- c('height','nodename')
M1D3_reducedforn <- cbind(M1D3_reducedfor,M1D3_reducedfinal)
M1D3_reducedbranlab <- M1D3_reduced$tip.label
M1D3_reducedforntip <- M1D3_reducedforn[1:length(M1D3_reducedbranlab),]
M1D3_mut <- vector()
M1D3_hei <- vector()

for (i in 1:length(M1D3_reducedbranlab)){
  
  id <- as.character(M1D3_reducedbranlab[i])
  idata <- M1D3_an[M1D3_an$Name.of.tips==id,]
  mut <- as.numeric(idata[2])
  heiuse <- M1D3_reducedforntip[M1D3_reducedforntip$label==id,]$height
  M1D3_mut <- append(M1D3_mut,mut)
  M1D3_hei <- append(M1D3_hei,heiuse)
  
}

M1D3_mutrate <- data.frame(M1D3_mut,M1D3_hei)

relation <- lm(M1D3_mut~M1D3_hei+0)
coef(relation)

# Coefficient: 124.9242

M1D3_reducedforntip$depth <- 124.9242*M1D3_reducedforntip$height
M1D3_reducedfornnode <- M1D3_reducedforn[2953:nrow(M1D3_reducedforn),]
M1D3_reducedfornnode$depth <- 124.9242*M1D3_reducedfornnode$height

M1D3_reducedforn$depth <- 124.9242*M1D3_reducedforn$height
tmpd <- abs(M1D3_reducedforn$depth-max(M1D3_reducedforn$depth))
M1D3samp_time <- sort(tmpd[1:length(M1D3_reducedbranlab)])
M1D3coal_time <- sort(tmpd[2953:length(tmpd)])

M1D3_reducedforn$coldepth <- 124.9242*(abs(M1D3_reducedforn$height-max(M1D3_reducedforn$height)))

# Input data for Phylodyn

M1D3list <- summarize_phylo(M1D3_reduced)
M1D3list[[1]] <- M1D3samp_time
M1D3list[[3]] <- M1D3coal_time

# Identify the largest generation of cell division based on the height of phylogenetic tree

max(M1D3_reducedforn$depth)

M1D3_cond <- BNPR(data = M1D3list, lengthout = 10)
M1D3_eps <- data.frame(M1D3_cond['x'],M1D3_cond['effpop'],M1D3_cond['effpop975'],M1D3_cond['effpop025'])
M1D3_eps$div <- M1D3_eps$effpop975-M1D3_eps$effpop025

M1D3_epslab <- seq(1, trunc(max(M1D3_eps$x)),length.out=5)
M1D3_epslabs <- as.character(rev(M1D3_epslab))

# plot

z1 <- ggplot(M1D3_eps,aes(x = x)) + geom_ribbon(aes(ymin=effpop025,ymax=effpop975),fill='#f5cac3')+geom_line(aes(y=effpop),linetype="solid",color="#f28482",size=1)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position = 'none')+scale_x_reverse(breaks=M1D3_epslab,labels=M1D3_epslabs)+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(1,2000))+ylab(NULL)+xlab(NULL)+ggtitle('S1_3d')
z1+theme_bw()+theme(panel.grid=element_blank(),axis.text.y = element_text(size = 12,face = "bold",color='black'),axis.text.x = element_text(size = 12,color='black'),axis.title.y = element_text(size = 16,vjust=2),axis.title.x = element_text(size = 16,vjust=-0.5),legend.position="none",plot.title = element_text(size=13, face="bold",hjust=0.5))