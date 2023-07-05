# Lineage Diversity--Time-Scaled Phylogenetic Tree

rm(list=ls())

library(ggtree)
library(treeio)
library(phytools)
library(phangorn)
library(ggplot2)
library(ggpubr)

load('~/.../data/Phylotime/Phylotime_WholeMidgut.RData')

# Take the time-scaled phylogenetic tree of first sample obtained on day 3 (S1_3d) as an example

M1D33_reducedfor <- as.data.frame(fortify(S1_3d))

M1D33_reducedHei <- nodeHeights(S1_3d)
M1D33_reducedHei_coor <- S1_3d$edge
nM1D33_reduced_1 <- c(M1D33_reducedHei[,1],M1D33_reducedHei[,2])
nM1D33_reduced_2 <- c(M1D33_reducedHei_coor[,1],M1D33_reducedHei_coor[,2])
nM1D33_reduced <- data.frame(nM1D33_reduced_1,nM1D33_reduced_2)
colnames(nM1D33_reduced) <- c('data','name')
M1D33_reducednew <- nM1D33_reduced[order(nM1D33_reduced$name),]
M1D33_reducedfinal <- unique(M1D33_reducednew)

colnames(M1D33_reducedfinal) <- c('height','nodename')
M1D33_reducedforn <- cbind(M1D33_reducedfor,M1D33_reducedfinal)
M1D33_reducedforn[M1D33_reducedforn$height==2,]
s1_3d_use <- M1D33_reducedforn[M1D33_reducedforn$height>2,]
s1_3d_use_tip <- s1_3d_use[s1_3d_use$isTip=='TRUE',]$label
s1_3d_use_node <- s1_3d_use[s1_3d_use$isTip=='FALSE',]$label

# From every branch to root
branlab_1 <- S1_3d$tip.label
nodelist_1 <- vector("list", length(branlab_1))
names(nodelist_1) <- branlab_1
rootno_1 <- getRoot(S1_3d)

for (i in 1:length(branlab_1)) {
  
  parlab <- M1D33_reducedfor[M1D33_reducedfor$label==branlab_1[i],]$parent
  templi <- vector()
  templi <- append(templi,branlab_1[i])
  
  while (parlab!=rootno_1) {
    
    parlab1 <- M1D33_reducedfor[M1D33_reducedfor$node==parlab,]$parent
    nodel <- M1D33_reducedfor[M1D33_reducedfor$node==parlab,]$label
    templi <- append(templi,nodel)
    parlab <- parlab1
    
  }
  
  templi <- append(templi,"root")
  nodelist_1[branlab_1[i]] <- list(templi)
  
}

table(2==lengths(nodelist_1)) 
nodelist_1_tip <- nodelist_1[s1_3d_use_tip]

# From every node to branch
totallist_1 <- vector("list", length(S1_3d$node.label))
names(totallist_1) <- S1_3d$node.label

for (i in 1:length(branlab_1)) {
  
  idata <- nodelist_1[[i]]
  
  if (length(idata)!=2) {
    
    for (j in 2:(length(idata)-1)) {
      
      jdata <- idata[j]
      tmp <- c(as.vector(unlist(totallist_1[jdata])),branlab_1[i])
      totallist_1[jdata] <- list(tmp)
      
    }
    
  }
  
}

totallist_1_node <- totallist_1[s1_3d_use_node]

# Find the highest node
progenitornode_1_new <- vector()
idata <- s1_3d_use_node
tmpdata <- unique(as.vector(unlist(totallist_1_node[s1_3d_use_node])))
jdata <- nodelist_1[tmpdata]

for (j in 1:length(jdata)) {
  
  nodedata <- jdata[[j]]
  tmpdata_1 <- intersect(nodedata,idata)
  progenitornode_1_new <- append(progenitornode_1_new,tmpdata_1[length(tmpdata_1)])
  
}

progenitornode_1_new <- unique(progenitornode_1_new)
progenitornode_1_newtip <- as.vector(unlist(totallist_1_node[progenitornode_1_new]))
progenitornode_1_newtip_other <- setdiff(s1_3d_use_tip,progenitornode_1_newtip)
s1_3d_lineage <- c(rep(1,length(progenitornode_1_newtip_other)),lengths(totallist_1_node[progenitornode_1_new]))

# Lineage richness 

s1_3d_q0 <- length(s1_3d_lineage)

# Shannon diversity

s1_3d_q1fre <- s1_3d_lineage/sum(s1_3d_lineage)
s1_3d_q1 <- exp(-sum(s1_3d_q1fre * log(s1_3d_q1fre)))

# Reciprocal of the maximum lineage frequency

s1_3d_qInf <- 1/(max(s1_3d_q1fre))