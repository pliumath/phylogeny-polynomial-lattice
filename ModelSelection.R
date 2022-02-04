library(treenomial)
library(ape)
library(apTreeshape)
library(phyloTop)
library(RPANDA)
library(dplyr)
library(phangorn)
library(phytools)
library(treetop)
library(cluster)
library(TreeSim)
library(castor)
library(e1071)
library(abc)
library(FNN)
library(KernelKnn)
library(tictoc)


summstats<- function(tree){
  
  #topology
  st1 <- colless.phylo(tree)
  st2 <- sackin.phylo(tree)
  st3 <- max(widths(tree)/max(node.depth(tree,method=2)))
  W<-widths(tree)
  st4 <- max(abs(W- c(0,W[1:length(W)-1])))
  st5 <- max(ladderSizes(tree)[[1]])
  st6 <- ILnumber(tree)
  st7 <- stairs(tree)
  
  #branch length
  rh <- cbind(c(tree$edge[,1],tree$edge[,2]),c(nodeHeights(tree)[,1],nodeHeights(tree)[,2]))
  rh <- as.data.frame(rh[!duplicated(rh),])
  rh <- rh[order(rh$V1),]
  
  rbl <- as.data.frame(rbind(cbind(tree$edge[,2],tree$edge.length),c(tree$Nnode+2,0)))
  rbl <- rbl[order(rbl$V1),]
  
  r <- cbind(rbl,rh[,2])
  names(r) <- c('V1','V2','V3')
  rex <- r[1:(tree$Nnode+1),]
  rin <- r[(tree$Nnode+3):(2*tree$Nnode+1),]
  
  sb1 <- max(rex[,3])
  sb2 <- min(rex[,3])
  sb3 <- mean(tree$edge.length)
  sb4 <- median(tree$edge.length)
  sb5 <- var(tree$edge.length)
  sb6 <- mean(rex[,2])
  sb7 <- median(rex[,2])
  sb8 <- var(rex[,2])
  sb9 <- c(mean(rin[rin[,3]<(sb1/3),2]),mean(rin[rin[,3]<(2*sb1/3) & rin[,3]>=(sb1/3),2]),mean(rin[rin[,3]>=(2*sb1/3),2]))
  sb10 <- c(median(rin[rin[,3]<(sb1/3),2]),median(rin[rin[,3]<(2*sb1/3) & rin[,3]>=(sb1/3),2]),median(rin[rin[,3]>=(2*sb1/3),2]))
  sb11 <- c(var(rin[rin[,3]<(sb1/3),2]),var(rin[rin[,3]<(2*sb1/3) & rin[,3]>=(sb1/3),2]),var(rin[rin[,3]>=(2*sb1/3),2]))
  sb12 <-sb9/sb6
  sb13 <-sb10/sb7
  sb14 <- sb11/sb8
  
  #ltt
  rltt <- ltt(tree,plot = FALSE)
  rl <- cbind(rltt$ltt,rltt$times)
  sl1 <- max(rl[,1])
  sl2 <- min(rl[rl[,1]==sl1,2])
  sl3 <- sl1/sl2
  sl4 <- (rl[nrow(rl),1]-sl1)/(max(rl[,2])-sl2)
  sl5 <- sl3/sl4
  
  sumdowntime <- 0
  sumdownnum <- 0
  for (i in 2:(nrow(rl)-2)) {
    
    if (rl[i,1]>rl[i+1,1] & rl[i+1,1]>rl[i+2,1]) {
      sumdowntime <- sumdowntime + rl[i+2,2]-rl[i,2]
      sumdownnum <- sumdownnum + 1
    }
    
  }
  sl6 <- sumdowntime/sumdownnum
  
  sut1 <- 0
  sun1 <- 0
  sut2 <- 0
  sun2 <- 0
  sut3 <- 0
  sun3 <- 0
  for (i in  2:(nrow(rl)-2)) {
    
    if (rl[i,1]<rl[i+1,1] & rl[i+1,1]<rl[i+2,1] & rl[i,2]<(sb1/3)) {
      sut1 <- sut1 + rl[i+2,2]-rl[i,2]
      sun1 <- sun1+ 1
    }
    
    if (rl[i,1]<rl[i+1,1] & rl[i+1,1]<rl[i+2,1] & rl[i,2]<(2*sb1/3) & rl[i,2]>=(sb1/3)) {
      sut2 <- sut2 + rl[i+2,2]-rl[i,2]
      sun2 <- sun2+ 1
    }
    
    if (rl[i,1]<rl[i+1,1] & rl[i+1,1]<rl[i+2,1] & rl[i,2]>=(2*sb1/3)) {
      sut3 <- sut3 + rl[i+2,2]-rl[i,2]
      sun3 <- sun3 + 1
    }
    
  }
  sl7 <- c(sut1/sun1,sut2/sun2,sut3/sun3)
  
  
  ss <- c(st1,st2,st3,st4,st5,st6,st7,sb1,sb2,sb3,sb4,sb5,sb6,sb7,sb8,sb9,sb10,sb11,sb12,sb13,sb14,sl1,sl2,sl3,sl4,sl5,sl6,sl7)
  names(ss) <- NULL
  ss[is.nan(ss)] <- 0
  ss[is.na(ss)] <- 0
  ss[is.infinite(ss)] <- 0
  return(ss)
  
}

###################################################
#generate trees

#loading data
exptrees <- read.tree("explosiveRadiation-tips500.nwk")
tratrees <- read.tree("traitEvolution-tips500.nwk")

numTrees <- 200
numSamps <- 100

MSTrees <- list()

for (i in 1:100) {
  
  print(i)
  
  #select trees
  samps <- sample(1:5000, numTrees, replace=FALSE)
  setr <- exptrees[samps]
  samps <- sample(1:5000, numTrees, replace=FALSE)
  sttr <- tratrees[samps]
  
  #select tests
  tsamps <- sample(5001:10000, numSamps, replace=FALSE)
  setest <- exptrees[tsamps]
  tsamps <- sample(5001:10000, numSamps, replace=FALSE)
  trtest <- tratrees[tsamps]
  
  Trees <- list("ExpRad" = setr,  "TraEvo" = sttr, "TestEx" = setest, "TestTr" = trtest)
  Trees <- lapply(unlist(Trees, recursive = FALSE, use.names = FALSE), as.phylo)
  names(Trees) <- c(rep("ExpRad", numTrees), rep("TraEvo", numTrees),rep("TestEx", numSamps),rep("TestTr", numSamps))
  
  MSTrees[[i]] <- Trees
  
}

save(MSTrees,file = "MSTrees.Rdata")


###################################################
#SS-Knn
load("MSTrees.Rdata")

MSSS <- matrix(0,nrow = 200,ncol = 200)

tic("SS")
for (i in 1:100) {
  
  print(i)
  
  A <- MSTrees[[i]]
  cls <- c(rep(1,200),rep(2,200))
  
  SS <- t(sapply(A, summstats))
  DSS <- as.matrix(dist(SS))
  
  ResSS <- distMat.KernelKnn(DSS,TEST_indices = 401:600, cls, regression = FALSE, Levels = c(1,2))
  
  MSSS[,(2*i-1):(2*i)] <- ResSS
  
}
Time <- toc()

save(MSSS,Time,file = "MSSS.Rdata") 


#################################
#LM-Knn
load("MSTrees.Rdata")

MSLM <- matrix(0,nrow = 200,ncol = 200)

tic("LM")
for (i in 1:100) {
  
  print(i)
  
  A <- MSTrees[[i]]
  cls <- c(rep(1,200),rep(2,200))
  
  DLM <- JSDtree(A)
  
  ResLM <- distMat.KernelKnn(DLM,TEST_indices = 401:600, cls, regression = FALSE, Levels = c(1,2))
  
  MSLM[,(2*i-1):(2*i)] <- ResLM
  
}
Time <- toc()

save(MSLM,Time,file = "MSLM.Rdata") 

#################################
#Latt-Knn

load("MSTrees.Rdata")

MSL <- matrix(0,nrow = 200,ncol = 200)

tic(paste("Lat",toString(k),sep = ""))
for (i in 1:100) {
  
  print(i)
  
  A <- MSTrees[[i]]
  cls <- c(rep(1,200),rep(2,200))
  
  L <- treeToLattice(A)
  L <- lapply(L, function(x){x[x[,3]==0 & x[,4]!=0, 3] <- .Machine$double.eps; return(x)})
  DL <- lattToDistMat(L)
  
  ResL <- distMat.KernelKnn(DL,TEST_indices = 401:600, cls, regression = FALSE, Levels = c(1,2))
  
  MSL[,(2*i-1):(2*i)] <- ResL
  
}
Time <- toc()

save(MSL,Time,file = "MSL.Rdata") 


#################################
#data and plots

library(treenomial)
library(ape)
library(apTreeshape)
library(ggplot2)
library(scales)
library(egg)
library(grid)
library(Rtsne)
library(phyloTop)
library(castor)
library(RPANDA)
library(reshape2)

load("MSSS.Rdata")
tss <- Time$toc - Time$tic
load("MSLM.Rdata")
tlm <- Time$toc - Time$tic
load("MSL.Rdata")
tl <- Time$toc - Time$tic

RSS <- colSums(MSSS[101:200,seq(1,200,2)])/200 + colSums(MSSS[1:100,seq(2,200,2)])/200 
RLM <- colSums(MSLM[101:200,seq(1,200,2)])/200 + colSums(MSLM[1:100,seq(2,200,2)])/200
RL <- colSums(MSL[101:200,seq(1,200,2)])/200 + colSums(MSL[1:100,seq(2,200,2)])/200

R <- as.data.frame(c(RL,RLM,RSS))
R[,2] <- c(rep("Lattice",100),rep("Lewitus-Morlon",100),rep("Summary Statistics",100))

names(R) <- c("M","D")

p <- ggplot(R, aes(x=D, y=M)) + 
  geom_jitter(shape=16, position=position_jitter(0.1),size=0.75) + geom_violin(alpha =0.5) +
  labs(x = "Distance", y = "Misclassification rate") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Misclassification rate in model selection for ER and TE trees by KNN classification") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)




df <- data.frame(c(tl,tlm,tss))
df[,2] <-  c("Lattice","Lewitus-Morlon","Summary Statistics")
names(df) <- c("time","dis")
p <- ggplot(data=df, aes(x=dis, y=time)) +
  geom_bar(stat="identity",width=0.5*3/5) +
  geom_text(aes(label=time,family="serif"), vjust=-0.3, size=3.5)+
  labs(x = "Distance", y = "Duration (s)") +
  ylim(0,50000)+
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Duration of the model selection experiment for ER and TE trees by KNN classification") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)



AA <- colSums(MSL[1:100,seq(1,200,2)])
AB <- colSums(MSL[1:100,seq(2,200,2)]) 
BA <- colSums(MSL[101:200,seq(1,200,2)]) 
BB <- colSums(MSL[101:200,seq(2,200,2)]) 

M <- matrix(c(mean(AA),mean(AB) ,mean(BA) ,mean(BB)),nrow = 2,ncol = 2) ####
rownames(M) <-  c("Explosive Radiation","Trait Evolution")
colnames(M) <-  c("Explosive Radiation","Trait Evolution")

MT <- melt(M)


p <- ggplot(data = MT, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high='#0072BD') +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) + 
  labs(x = "Generating model", y = "Estimated model", color = "Model:")+ theme(legend.position = "none") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Average percentages of estimated models for ER and TE trees with lattice distance") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12)) + theme(axis.text.y = element_text(angle = 90, hjust = .5))


tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)




#################################
#################################
#SUPPs


AA <- colSums(MSLM[1:100,seq(1,200,2)])
AB <- colSums(MSLM[1:100,seq(2,200,2)]) 
BA <- colSums(MSLM[101:200,seq(1,200,2)]) 
BB <- colSums(MSLM[101:200,seq(2,200,2)]) 

M <- matrix(c(mean(AA),mean(AB) ,mean(BA),mean(BB)),nrow = 2,ncol = 2)
rownames(M) <-  c("Explosive Radiation","Trait Evolution")
colnames(M) <-  c("Explosive Radiation","Trait Evolution")

MT <- melt(M)


p <- ggplot(data = MT, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high='#D95319') +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) + 
  labs(x = "Generating model", y = "Estimated model", color = "Model:")+ theme(legend.position = "none") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Average percentages of estimated models for ER and TE trees with Lewitus-Morlon distance") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12)) + theme(axis.text.y = element_text(angle = 90, hjust = .5))


tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)




AA <- colSums(MSSS[1:100,seq(1,200,2)])
AB <- colSums(MSSS[1:100,seq(2,200,2)]) 
BA <- colSums(MSSS[101:200,seq(1,200,2)]) 
BB <- colSums(MSSS[101:200,seq(2,200,2)]) 

M <- matrix(c(mean(AA),mean(AB) ,mean(BA),mean(BB)),nrow = 2,ncol = 2)
rownames(M) <-  c("Explosive Radiation","Trait Evolution")
colnames(M) <-  c("Explosive Radiation","Trait Evolution")

MT <- melt(M)


p <- ggplot(data = MT, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high='#EDB120') +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) + 
  labs(x = "Generating model", y = "Estimated model", color = "Model:")+ theme(legend.position = "none") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Average percentages of estimated models for ER and TE trees with summary statistics distance") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12)) + theme(axis.text.y = element_text(angle = 90, hjust = .5))


tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)

