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



lattostats <- function(L){
  
  S <- sapply(L, dim)
  m <- max(S[2,])
  
  if (m == 5) {
    
    M <- sapply(L, function(x){
      
      y <- cbind(x[,5],x[,3])
      y <- y[order(y[,1]),]
      y <- c(y[,1],y[,2])
      
      return(y)
      
    })
    
    M <- t(M)
    
  }
  
  if (m == 6) {
    
    M <- sapply(L, function(x){
      
      n1 <- dim(x)[1]
      n2 <- dim(x)[2]
      
      if (n2 == 6) {
        
        y <- cbind(x[,5],x[,6],x[,3])
        y <- y[order(y[,1]),]
        y <- y[order(y[,2]),]
        y <- c(y[,1],y[,2],y[,3])
        
      }
      
      if (n2 ==5) {
        
        y <- cbind(x[,5],rep(0,n1),x[,3])
        y <- y[order(y[,1]),]
        y <- c(y[,1],y[,2],y[,3])
        
      }
      
      return(y)
      
    })
    
    M <- t(M)
    
  }
  
  return(M)
  
}


meanrelerr <- function(R){
  
  RE <- abs(R-TruePars)/TruePars
  MRE <- matrix(colMeans(RE),nrow = 100,ncol = 10)
  
}

#################################
#generate datasets

numtrees <- 400
numestim <- 100
ne <- numtrees - numestim

for (k in 1:10) {
  
  Paras <- list()
  BDTrees <- list()
  
  for (i in 1:100) {
    
    print(i)
    
    para <- runif(numtrees,1,8)
    trees <- list()
    
    j <- 1
    while (j <= numtrees) {
      
      x <- para[j]
      tree <- sim.bd.taxa(n=50*k, numbsim=1 ,mu=1, lambda=x, complete = FALSE)
      tree <- tree[[1]]
      
      if (tree$Nnode == 50*k - 1) {
        
        trees[[j]] <- as.phylo(tree)
        j <- j + 1
        
      }
      
    }
    
    Paras[[i]] <- para
    BDTrees[[i]] <- trees
    
  }
  
  save(Paras,BDTrees,file = paste("Set",toString(k),".Rdata",sep = ""))
  
}


###################################################
#SS-Knn
for (k in 1:10) {
  load(paste("Set",toString(k),".Rdata",sep = ""))
  
  ResSS <- matrix(0,nrow = 100,ncol = 100) 
  
  tic(paste("SS",toString(k),sep = ""))
  
  for (i in 1:100) {
    
    print(i)
    
    A <- BDTrees[[i]] 
    para <- Paras[[i]] 
    
    SS <- t(sapply(A, summstats))
    DSS <- as.matrix(dist(SS))
    
    EstSS <- distMat.KernelKnn(DSS,TEST_indices = (ne+1):numtrees, para[1:ne], regression = TRUE)
    
    ResSS[,i] <- EstSS 
    
  }
  
  Time <- toc() 
  
  save(ResSS,Time,file = paste("ResSS",toString(k),".Rdata",sep = "")) 
}

#################################
#LM-Knn
for (k in 1:10) {
  
  load(paste("Set",toString(k),".Rdata",sep = ""))
  
  ResLM <- matrix(0,nrow = 100,ncol = 100) ##
  
  tic(paste("LM",toString(k),sep = ""))
  
  for (i in 1:100) {
    
    print(i)
    
    A <- BDTrees[[i]] 
    para <- Paras[[i]] 
    
    DLM <- JSDtree(A)
    
    EstLM <- distMat.KernelKnn(DLM,TEST_indices = (ne+1):numtrees, para[1:ne], regression = TRUE)
    
    ResLM[,i] <- EstLM 
    
  }
  
  Time <- toc()
  save(ResLM,Time,file = paste("ResLM",toString(k),".Rdata",sep = "")) 
  
}

#################################
#Latt-Knn
for (k in 1:10) {
  
  load(paste("Set",toString(k),".Rdata",sep = ""))
  
  ResL <- matrix(0,nrow = 100,ncol = 100) ##
  
  tic(paste("Lat",toString(k),sep = ""))
  
  for (i in 1:100) {
    
    print(i)
    
    A <- BDTrees[[i]] 
    para <- Paras[[i]] 
    
    L <- treeToLattice(A)
    DL <- lattToDistMat(L)
    
    EstL <- distMat.KernelKnn(DL,TEST_indices = (ne+1):numtrees, para[1:ne], regression = TRUE)
    
    ResL[,i] <- EstL
    
  }
  
  Time <- toc() 
  save(ResL,Time,file = paste("ResL",toString(k),".Rdata",sep = ""))
}



###################################################
#SS-abc
for (k in 1:10) {
  load(paste("Set",toString(k),".Rdata",sep = ""))
  
  ResSS <- list()
  
  tic(paste("SS",toString(k),sep = ""))
  
  for (i in 1:100) {
    
    print(i)
    
    A <- BDTrees[[i]] 
    para <- Paras[[i]] 
    
    SS <- t(sapply(A, summstats))
    DSS <- as.matrix(dist(SS))
    
    EstSS <- lapply((ne+1):numtrees, function(x){abc(SS[x,],para[1:ne],SS[1:ne,],tol = .1,method = "rejection")})
    
    ResSS[[i]] <- EstSS 
    
  }
  
  Time <- toc() 
  
  save(ResSS,Time,file = paste("ABCSS",toString(k),".Rdata",sep = "")) 
}



#################################
#Latt-abc
for (k in 1:10) {
  
  load(paste("Set",toString(k),".Rdata",sep = ""))
  
  ResL <- list()
  
  tic(paste("Lat",toString(k),sep = ""))
  
  for (i in 1:100) {
    
    print(i)
    
    A <- BDTrees[[i]] 
    para <- Paras[[i]] 
    
    L <- treeToLattice(A)
    
    M <- lattostats(L)
    
    EstL <- lapply((ne+1):numtrees, function(x){abc(M[x,],para[1:ne],M[1:ne,],tol = .1,method = "rejection")})
    
    ResL[[i]] <- EstL
    
  }
  
  Time <- toc() 
  save(ResL,Time,file = paste("ABCL",toString(k),".Rdata",sep = ""))
}

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

#true parameters
TruePars <- matrix(0,nrow = 100,ncol = 1000)

for (k in 1:10) {
  
  load(paste("Set",toString(k),".Rdata",sep = ""))
  
  for (i in 1:100) {
    
    para <- Paras[[i]] 
    TruePars[,(k-1)*100+i] <- para[(ne+1):numtrees]
    
  }
  
}

#Knn
RSS <- matrix(0,nrow = 100,ncol = 1000)
tss <- c()

for (k in 1:10) {
  
  load(paste("ResSS",toString(k),".Rdata",sep = ""))
  
  RSS[,((k-1)*100+1):(k*100)] <- ResSS
  tss[k] <- Time$toc - Time$tic
  
}

RLM <- matrix(0,nrow = 100,ncol = 1000)
tlm <- c()

for (k in 1:10) {
  
  load(paste("ResLM",toString(k),".Rdata",sep = ""))
  
  RLM[,((k-1)*100+1):(k*100)] <- ResLM
  tlm[k] <- Time$toc - Time$tic
  
}

RL <- matrix(0,nrow = 100,ncol = 1000)
tl <- c()

for (k in 1:10) {
  
  load(paste("ResL",toString(k),".Rdata",sep = ""))

  RL[,((k-1)*100+1):(k*100)] <- ResL
  tl[k] <- Time$toc - Time$tic
  
}



MRESS <- meanrelerr(RSS)
MRELM <- meanrelerr(RLM)
MREL <- meanrelerr(RL)


ntips <- rep(seq(50,500,50),3)
dists <- c(rep("Lattice",10),rep("Lewitus-Morlon",10),rep("Summary Statistics",10))

me <- c(apply(MREL,2,function(x){return(quantile(x,0.5))}),apply(MRELM,2,function(x){return(quantile(x,0.5))}),apply(MRESS,2,function(x){return(quantile(x,0.5))}))
mi <- me - c(qnorm(0.975)*apply(MREL,2,sd)/sqrt(100),qnorm(0.975)*apply(MRELM,2,sd)/sqrt(100),qnorm(0.975)*apply(MRESS,2,sd)/sqrt(100))
ma <- me - c(qnorm(0.975)*apply(MREL,2,sd)/sqrt(100),qnorm(0.975)*apply(MRELM,2,sd)/sqrt(100),qnorm(0.975)*apply(MRESS,2,sd)/sqrt(100))

Data <- as.data.frame(cbind(ntips,me,mi,ma))
Data[,5] <- dists

pal <- c('#0072BD','#D95319','#EDB120','#7E2F8E')

p <- ggplot(Data,aes(x=ntips,y=me)) +
  geom_smooth(aes(ymin=mi,ymax=ma,fill=dists,color=dists))+
  scale_color_manual(values=pal[c(1,2,3)])+
  scale_fill_manual(values=pal[c(1,2,3)])+
  labs(x = "Number of tips", y = "Average relative error",fill = "Distance:", color = "Distance:")+
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12)) + 
  ggtitle("Average relative error in parameter estimation by KNN regression") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)


#
df <- data.frame(cbind(ntips,c(tl,tlm,tss)))
df[1:10,3] <- df[1:10,1] - 13.5
df[11:20,3] <- df[11:20,1]
df[21:30,3] <- df[21:30,1] + 13.5
df[,4] <- dists

names(df) <- c("ntips","time","tx","dists")

p <- ggplot(data=df, aes(x=ntips, y=time, fill = dists)) +
  geom_bar(stat="identity",width=40, position=position_dodge()) +
  geom_text(aes(x = tx,label=time,family="serif"), angle = -90 , hjust= 1.05, size=2.5)+
  scale_fill_manual(values=pal[c(1,2,3)])+
  labs(x = "Number of tips", y = "Duration (s)",fill = "Distance:") +
  ylim(0,30000) +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Duration of the experiment of estimating parameter by KNN regression") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)


#
ind <- 998
TP <- rep(TruePars[,ind],3)
EP <- c(RL[,ind],RLM[,ind],RSS[,ind])

R <- as.data.frame(cbind(TP,EP))

R[,3] <- c(rep("Lattice",100),rep("Lewitus-Morlon",100),rep("Summary Statistics",100))

p <- ggplot(R, aes(x = TP, y = EP, color = V3)) +
  geom_abline() +
  geom_point() +
  xlim(1,8)+
  ylim(1,8)+
  scale_color_manual(values=c('#0072BD', '#D95319', '#EDB120'))+
  labs(x = "Real value", y = "Estimated value", color = "Distance:") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Real vs. estimated values of birth rate for 500-tip trees") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))


tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)




###################################################
###################################################
#SUPPs


RL <- matrix(0,nrow = 100,ncol = 1000)
tl <- c()

for (k in 1:10) {
  
  load(paste("ResL",toString(k),".Rdata",sep = ""))
  
  RL[,((k-1)*100+1):(k*100)] <- ResL
  tl[k] <- Time$toc - Time$tic
  
}


RSSabc <- matrix(0,nrow = 100,ncol = 1000)
tssabc <- c()

for (k in 1:10) {
  
  load(paste("ABCSS",toString(k),".Rdata",sep = ""))
  
  Res <-  matrix(0,nrow = 100,ncol = 100)
  for (j in 1:100) {
    print(j)
    y <- ResSS[[j]]
      for (i in 1:100) {
        Res[i,j] <- mean(y[[i]]$unadj.values)
      }
  }
  
  RSSabc[,((k-1)*100+1):(k*100)] <- Res
  tssabc[k] <- Time$toc - Time$tic
  
}



RLabc <- matrix(0,nrow = 100,ncol = 1000)
tlabc <- c()

for (k in 1:10) {
  
  load(paste("ABCL",toString(k),".Rdata",sep = ""))
  
  Res <-  matrix(0,nrow = 100,ncol = 100)
  for (j in 1:100) {
    print(j)
    y <- ResL[[j]]
    for (i in 1:100) {
      Res[i,j] <- mean(y[[i]]$unadj.values)
    }
  }
  
  RLabc[,((k-1)*100+1):(k*100)] <- Res
  tlabc[k] <- Time$toc - Time$tic
  
}




MRESS <- meanrelerr(RSSabc)
MRELM <- meanrelerr(RLabc)
MREL <- meanrelerr(RL)


ntips <- rep(seq(50,500,50),3)
dists <- c(rep("Lattice-ABC",10),rep("Lattice-KNN",10),rep("Summary Statistics-ABC",10))

me <- c(apply(MRELM,2,function(x){return(quantile(x,0.5))}),apply(MREL,2,function(x){return(quantile(x,0.5))}),apply(MRESS,2,function(x){return(quantile(x,0.5))}))
mi <- me - c(qnorm(0.975)*apply(MREL,2,sd)/sqrt(100),qnorm(0.975)*apply(MRELM,2,sd)/sqrt(100),qnorm(0.975)*apply(MRESS,2,sd)/sqrt(100))
ma <- me - c(qnorm(0.975)*apply(MREL,2,sd)/sqrt(100),qnorm(0.975)*apply(MRELM,2,sd)/sqrt(100),qnorm(0.975)*apply(MRESS,2,sd)/sqrt(100))

Data <- as.data.frame(cbind(ntips,me,mi,ma))
Data[,5] <- dists

pal <- c('#7E2F8E','#0072BD','#77AC30','#D95319')

p <- ggplot(Data,aes(x=ntips,y=me)) +
  geom_smooth(aes(ymin=mi,ymax=ma,fill=dists,color=dists))+
  scale_color_manual(values=pal[c(1,2,3)])+
  scale_fill_manual(values=pal[c(1,2,3)])+
  labs(x = "Number of tips", y = "Average relative error",fill = "Method:", color = "Method:")+
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12)) + 
  ggtitle("Average relative error in parameter estimation by diffenerent methods") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)




df <- data.frame(cbind(ntips,c(tlabc,tl,tssabc)))
df[1:10,3] <- df[1:10,1] - 13.5
df[11:20,3] <- df[11:20,1]
df[21:30,3] <- df[21:30,1] + 13.5
df[,4] <- dists

names(df) <- c("ntips","time","tx","dists")

p <- ggplot(data=df, aes(x=ntips, y=time, fill = dists)) +
  geom_bar(stat="identity",width=40, position=position_dodge()) +
  geom_text(aes(x = tx,label=time,family="serif"), angle = -90 , hjust= 1.05, size=2.5)+
  scale_fill_manual(values=pal[c(1,2,3)])+
  labs(x = "Number of tips", y = "Duration (s)",fill = "Method:") +
  ylim(0,8000) +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Duration of the experiment of estimating parameter by different methods") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)

