library(ape)
library(apTreeshape)
library(phyloTop)
library(RPANDA)
library(dplyr)
library(phangorn)
library(phytools)
library(treetop)
library(cluster)
library(treenomial)
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
  # sb2 <- min(rex[,3])
  # sb3 <- mean(tree$edge.length)
  # sb4 <- median(tree$edge.length)
  # sb5 <- var(tree$edge.length)
  # sb6 <- mean(rex[,2])
  # sb7 <- median(rex[,2])
  # sb8 <- var(rex[,2])
  # sb9 <- c(mean(rin[rin[,3]<(sb1/3),2]),mean(rin[rin[,3]<(2*sb1/3) & rin[,3]>=(sb1/3),2]),mean(rin[rin[,3]>=(2*sb1/3),2]))
  # sb10 <- c(median(rin[rin[,3]<(sb1/3),2]),median(rin[rin[,3]<(2*sb1/3) & rin[,3]>=(sb1/3),2]),median(rin[rin[,3]>=(2*sb1/3),2]))
  # sb11 <- c(var(rin[rin[,3]<(sb1/3),2]),var(rin[rin[,3]<(2*sb1/3) & rin[,3]>=(sb1/3),2]),var(rin[rin[,3]>=(2*sb1/3),2]))
  # sb12 <-sb9/sb6
  # sb13 <-sb10/sb7
  # sb14 <- sb11/sb8
  
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
  
  
  # ss <- c(st1,st2,st3,st4,st5,st6,st7,sb1,sb2,sb3,sb4,sb5,sb6,sb7,sb8,sb9,sb10,sb11,sb12,sb13,sb14,sl1,sl2,sl3,sl4,sl5,sl6,sl7)
  ss <- c(st1,st2,st3,st4,st5,st6,st7,sl1,sl2,sl3,sl4,sl5,sl6,sl7)
  names(ss) <- NULL
  ss[is.nan(ss)] <- 0
  ss[is.na(ss)] <- 0
  ss[is.infinite(ss)] <- 0
  return(ss)
  
}

misrate <- function(K,numTrees){
  
  R1 <- sum(K != c(rep(1,numTrees),rep(2,numTrees),rep(3,numTrees)))/(3*numTrees)
  R2 <- sum(K != c(rep(1,numTrees),rep(3,numTrees),rep(2,numTrees)))/(3*numTrees)
  R3 <- sum(K != c(rep(2,numTrees),rep(1,numTrees),rep(3,numTrees)))/(3*numTrees)
  R4 <- sum(K != c(rep(2,numTrees),rep(3,numTrees),rep(1,numTrees)))/(3*numTrees)
  R5 <- sum(K != c(rep(3,numTrees),rep(1,numTrees),rep(2,numTrees)))/(3*numTrees)
  R6 <- sum(K != c(rep(3,numTrees),rep(2,numTrees),rep(1,numTrees)))/(3*numTrees)
  
  return(min(c(R1,R2,R3,R4,R5,R6)))
  
}

#####################################################

#generate a list of data sets of tree shapes
TreeShapes <- list()

for (i in 1:100) {
  
  print(i)
  
  numTrees = 100
  numTips = 100
  
  sim1 <- rtreeshape(numTrees, tip.number = numTips, model = "yule")
  sim2 <- rtreeshape(numTrees, tip.number = numTips, model = "aldous")
  sim3 <- rtreeshape(numTrees, tip.number = numTips, model = "pda")
  
  modelTrees <- list("Yule" = sim1,  "Aldous" = sim2, "PDA" = sim3)
  modelTrees <- lapply(unlist(modelTrees, recursive = FALSE, use.names = FALSE), as.phylo)
  names(modelTrees) <- c(rep("Yule", numTrees), rep("Aldous", numTrees), rep("PDA", numTrees))
  
  TreeShapes[[i]] <- lapply(modelTrees,function(x){x$edge.length <- rep(1,1,numTips*2-2);return(x)})
  
}

save(TreeShapes,file = "TreeShapes.Rdata")


#clustering the tree shapes with label D2 distance
load("TreeShapes.Rdata")

RCK <- c()

tic("CK")
for (i in 1:100) {
  
  print(i)
  
  A <- TreeShapes[[i]]
  
  DCK <-vecMultiDistUnlab(A)
  KCK <- pam(DCK,3)
  KCK <- KCK$clustering
  RCK <- cbind(RCK,misrate(KCK,100))
  
}
time <- toc()

save(RCK,time,file = "RCK_TS.Rdata")


#clustering the tree shapes with the lattice distance
load("TreeShapes.Rdata")

RL <- c()

tic("Lat")
for (i in 1:100) {
  
  print(i)
  
  A <- TreeShapes[[i]]
  
  L <- treeToLattice(A)
  DL <- lattToDistMat(L)
  
  KL <- pam(DL,3)
  KL <- KL$clustering
  RL <- cbind(RL,misrate(KL,100))
  
}
time <- toc()

save(RL,time,file = "RL_TS.Rdata")


#clustering the tree shapes with LM distance
load("TreeShapes.Rdata")

RLM <- c()

tic("LM")
for (i in 1:100) {
  
  print(i)
  
  A <- TreeShapes[[i]]
  
  DLM <- JSDtree(A)
  KLM <- pam(DLM,3)
  KLM <- KLM$clustering
  RLM <- cbind(RLM,misrate(KLM,100))
  
}
time <- toc()

save(RLM,time,file = "RLM_TS.Rdata")


#clustering the tree shapes with the polynomial distance
load("TreeShapes.Rdata")

RP <- c()

tic("Poly")
for (i in 1:100) {
  
  print(i)
  
  A <- TreeShapes[[i]]
  
  P <- treeToPoly(A)
  DP <- polyToDistMat(P)
  
  KP <- pam(DP,3)
  KP <- KP$clustering
  RP <- cbind(RP,misrate(KP,100))
  
}
time <- toc()

save(RP,time,file = "RP_TS.Rdata")


#clustering the tree shapes with SS-Euc distance
load("TreeShapes.Rdata")

RSS <- c()

tic("SS")
for (i in 1:100) {
  
  print(i)
  A <- TreeShapes[[i]]
  
  SS <- t(sapply(A, summstats))
  
  DSS <- as.matrix(dist(SS))
  
  KSS <- pam(DSS,3)
  KSS <- KSS$clustering
  RSS <- cbind(RSS,misrate(KSS,100))
  
}
time <- toc()

save(RSS,time,file = "RSS_TS.Rdata")

#####################################################

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


#MCR

load("RCK_TS.Rdata")
R1<-RCK
load("RL_TS.Rdata")
R2<-RL
load("RLM_TS.Rdata")
R3<-RLM
load("RP_TS.Rdata")
R4<-RP
load("RSS_TS.Rdata")
R5<-RSS

R <- as.data.frame(t(rbind(R1,R2,R3,R4,R5)))
names(R) <- c("Label D2","Lattice","Lewitus-Morlon","Polynomial","Summary Statistics")

R <- melt(R)
names(R) <- c("D","M")

p <- ggplot(R, aes(x=D, y=M)) +
  geom_jitter(shape=16, position=position_jitter(0.1),size=0.75) + geom_violin(alpha =0.5) +
  labs(x = "Distance", y = "Misclassification rate") +
  ylim(0,0.5)+
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Misclassification rate of clustering beta-splitting tree shapes by k-medoids") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)



#MDS
fit <- cmdscale(DP, k = 3)

methodResults <- data.frame(fit, "color" = names(A))
p <- ggplot(methodResults, aes(x = X1, y = X2, color = color)) +
  geom_point() +
  scale_color_manual(values=c('#0072BD', '#D95319', '#EDB120'))+
  labs(x = "Dim. 1", y = "Dim. 2", color = "Model:") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Visualization of polynomial distances between beta-splitting tree shapes") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))


tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)




#time
load("RCK_TS.Rdata")
R1<-time
t1 <- R1$toc-R1$tic
load("RL_TS.Rdata")
R2<-time
t2 <- R2$toc-R2$tic
load("RLM_TS.Rdata")
R3<-time
t3 <- R3$toc-R3$tic
load("RP_TS.Rdata")
R4<-time
t4 <- R4$toc-R4$tic
load("RSS_TS.Rdata")
R5<-time
t5 <- R5$toc-R5$tic

df <- data.frame(c(t1,t2,t3,t4,t5))
df[,2] <- c("Label D2","Lattice","Lewitus-Morlon","Polynomial","Summary Statistics")
names(df) <- c("time","dis")

p <- ggplot(data=df, aes(x=dis, y=time)) +
  geom_bar(stat="identity",width=0.5) +
  geom_text(aes(label=time,family="serif"), vjust=-0.3, size=3.5)+
labs(x = "Distance", y = "Duration (s)") +
  ylim(0,1000)+
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Duration of the experiment of clustering beta-splitting tree shapes") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)

#####################################################
#####################################################
#SUPPs


fit <- cmdscale(DCK, k = 3)

methodResults <- data.frame(fit, "color" = names(A))
p <- ggplot(methodResults, aes(x = X1, y = X2, color = color)) +
  geom_point() +
  scale_color_manual(values=c('#0072BD', '#D95319', '#EDB120'))+
  labs(x = "Dim. 1", y = "Dim. 2", color = "Model:") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Visualization of label d2 distances between beta-splitting tree shapes") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))


tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)






fit <- cmdscale(DL, k = 3)

methodResults <- data.frame(fit, "color" = names(A))
p <- ggplot(methodResults, aes(x = X1, y = X2, color = color)) +
  geom_point() +
  scale_color_manual(values=c('#0072BD', '#D95319', '#EDB120'))+
  labs(x = "Dim. 1", y = "Dim. 2", color = "Model:") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Visualization of lattice distances between beta-splitting tree shapes") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))


tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)







fit <- cmdscale(DLM, k = 3)

methodResults <- data.frame(fit, "color" = names(A))
p <- ggplot(methodResults, aes(x = X1, y = X2, color = color)) +
  geom_point() +
  scale_color_manual(values=c('#0072BD', '#D95319', '#EDB120'))+
  labs(x = "Dim. 1", y = "Dim. 2", color = "Model:") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Visualization of Lewitus-Morlon distances between beta-splitting tree shapes") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))


tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)




fit <- cmdscale(DSS, k = 3)

methodResults <- data.frame(fit, "color" = names(A))
p <- ggplot(methodResults, aes(x = X1, y = X2, color = color)) +
  geom_point() +
  scale_color_manual(values=c('#0072BD', '#D95319', '#EDB120'))+
  labs(x = "Dim. 1", y = "Dim. 2", color = "Model:") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Visualization of summary statistics distances between beta-splitting tree shapes") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))


tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)
