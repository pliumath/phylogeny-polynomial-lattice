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

treelength <- function(x){
  
  edgeinfo <- cbind(x$edge,x$edge.length)
  
  paths <- nodepath(x)
  
  pathlength <- function(p){
    
    l <- sapply(1:(length(p)-1), function(i){edgeinfo[edgeinfo[,1] == p[i] & edgeinfo[,2] == p[i+1],3]})
    
    return(sum(l))
    
  }
  
  m <- max(sapply(paths,pathlength))
  
  return(m)
  
}

normalizetree <- function(x){
  
  edgeinfo <- cbind(x$edge,x$edge.length)
  
  paths <- nodepath(x)
  
  pathlength <- function(p){
    
    l <- sapply(1:(length(p)-1), function(i){edgeinfo[edgeinfo[,1] == p[i] & edgeinfo[,2] == p[i+1],3]})
    
    return(sum(l))
    
  }
  
  m <- max(sapply(paths,pathlength))
  
  x$edge.length <- x$edge.length * (0.5/m)
  
  return(x)
  
}

genbd <- function(para,numtips){
  
  n <- length(para)
  sims <- list()
  
  for (i in 1:n) {
    
    x <- para[i]
    j <- 1
    
    while (j == 1) {
      
      tree <- sim.bd.taxa(n=ceiling(x*numtips), numbsim=1 ,mu=1, lambda=x, complete = TRUE)
      tree <- tree[[1]]
      tree <- drop.tip(tree,1:ceiling(x*numtips))
      
      if (tree$Nnode > numtips - 1) {
        
        tree <- drop.random(tree,tree$Nnode + 1 - numtips)
        sims[[i]] <- as.phylo(tree)
        j <- 0
        
      }
      
    }
  }
  
  return(sims)
  
}

#####################################
#####################################
load("who_clade_trees.Rdata")

N <- as.phylo(cladeA1b135N)

Res <- c()

numTips <- N$Nnode + 1
n <- 100

k <- 1

tic("A1b135N")
for (i in 1:100) {
  
  print(i)
  
  para <- runif(n,1,10)
  
  A <- genbd(para,numTips)
  A[[101]] <- N
  
  L <- treeToLattice(A)
  L <- lapply(L, function(x){x[x[,3]==0 & x[,4]!=0, 3] <- .Machine$double.eps; return(x)})
  DL <- lattToDistMat(L)
  
  Res[i] <- distMat.KernelKnn(DL,TEST_indices = 101, para, regression = TRUE)
  
}
Time <- toc()

save(Res,Time,file = paste("Flu",toString(k),".Rdata",sep = "")) 

############################################3
N <- as.phylo(cladeA1b135K)

Res <- c()

numTips <- N$Nnode + 1
n <- 100

k <- 2

tic("A1b135K")
for (i in 1:100) {
  
  print(i)
  
  para <- runif(n,1,10)
  
  A <- genbd(para,numTips)
  A[[101]] <- N
  
  L <- treeToLattice(A)
  L <- lapply(L, function(x){x[x[,3]==0 & x[,4]!=0, 3] <- .Machine$double.eps; return(x)})
  DL <- lattToDistMat(L)
  
  Res[i] <- distMat.KernelKnn(DL,TEST_indices = 101, para, regression = TRUE)
  
}
Time <- toc()

save(Res,Time,file = paste("Flu",toString(k),".Rdata",sep = "")) 



######################################
N <- as.phylo(clade3c3.B)

Res <- c()

numTips <- N$Nnode + 1
n <- 100

k <- 3

tic("3c3.B")
for (i in 1:100) {
  
  print(i)
  
  para <- runif(n,1,10)
  
  A <- genbd(para,numTips)
  A[[101]] <- N
  
  L <- treeToLattice(A)
  L <- lapply(L, function(x){x[x[,3]==0 & x[,4]!=0, 3] <- .Machine$double.eps; return(x)})
  DL <- lattToDistMat(L)
  
  Res[i] <- distMat.KernelKnn(DL,TEST_indices = 101, para, regression = TRUE)
  
}
Time <- toc()

save(Res,Time,file = paste("Flu",toString(k),".Rdata",sep = "")) 

###############################

N <- as.phylo(cladeA3)

Res <- c()

numTips <- N$Nnode + 1
n <- 100

k <- 4

tic("A3")
for (i in 1:100) {
  
  print(i)
  
  para <- runif(n,1,10)
  
  A <- genbd(para,numTips)
  A[[101]] <- N
  
  L <- treeToLattice(A)
  L <- lapply(L, function(x){x[x[,3]==0 & x[,4]!=0, 3] <- .Machine$double.eps; return(x)})
  DL <- lattToDistMat(L)
  
  Res[i] <- distMat.KernelKnn(DL,TEST_indices = 101, para, regression = TRUE)
  
}
Time <- toc()

save(Res,Time,file = paste("Flu",toString(k),".Rdata",sep = "")) 

###############################

Hunt <- read.tree("Hunt.nwk")
H <- as.phylo(midpoint(Hunt))

Res <- c()

numTips <- H$Nnode + 1
n <- 100

tic("Hunt")
for (i in 1:100) {
  
  print(i)
  
  para <- runif(n,1,10)
  
  A <- genbd(para,numTips)
  A[[101]] <- N
  
  A <- lapply(A, normalizetree)
  
  L <- treeToLattice(A)
  DL <- lattToDistMat(L)
  
  Res[i] <- distMat.KernelKnn(DL,TEST_indices = 101, para, regression = TRUE)
  
}
Time <- toc()

save(Res,Time,file = "Hunt.Rdata") 

###############################


Novitsky <- read.nexus("Novitsky.nex")
N <- as.phylo(midpoint(Novitsky))

Res <- c()

numTips <- N$Nnode + 1
n <- 100

tic("Nov")
for (i in 1:100) {
  
  print(i)
  
  para <- runif(n,1,10)
  
  A <- genbd(para,numTips)
  A[[101]] <- N
  
  A <- lapply(A, normalizetree)
  
  L <- treeToLattice(A)
  DL <- lattToDistMat(L)
  
  Res[i] <- distMat.KernelKnn(DL,TEST_indices = 101, para, regression = TRUE)
  
}
Time <- toc()

save(Res,Time,file = "Nov.Rdata") 

###############################

Wolf <- read.tree("Wolf.nwk")
W <- as.phylo(midpoint(Wolf))

Res <- c()

numTips <- W$Nnode + 1
n <- 100

tic("Wolf")
for (i in 1:100) {
  
  print(i)
  
  para <- runif(n,1,10)
  
  A <- genbd(para,numTips)
  A[[101]] <- N
  
  A <- lapply(A, normalizetree)
  
  L <- treeToLattice(A)
  DL <- lattToDistMat(L)
  
  Res[i] <- distMat.KernelKnn(DL,TEST_indices = 101, para, regression = TRUE)
  
}
Time <- toc()

save(Res,Time,file = "Wolf.Rdata") 


###############################
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

load("Hunt.Rdata")
RH <- Res
th <- Time$toc-Time$tic

load("Nov.Rdata")
RN <- Res
tn <- Time$toc-Time$tic

load("Wolf.Rdata")
RW <- Res
tw <- Time$toc-Time$tic

load("Flu1.Rdata")
R1 <- Res
t1 <- Time$toc-Time$tic
Time$msg

load("Flu2.Rdata")
R2 <- Res
t2 <- Time$toc-Time$tic
Time$msg

load("Flu3.Rdata")
R3 <- Res
t3 <- Time$toc-Time$tic
Time$msg

load("Flu4.Rdata")
R4 <- Res
t4 <- Time$toc-Time$tic
Time$msg

R <- as.data.frame(cbind(RH,RN,RW,R1,R2,R3,R4))
names(R) <- c("Hunt","Novitsky","Wolf","A1B/135N","A1B/135K","3c3.B","A3")

R <- melt(R)
names(R) <- c("D","M")


p <- ggplot(R, aes(x=D, y=M)) + 
  geom_jitter(shape=16, position=position_jitter(0.1),size=0.75) + geom_violin(alpha =0.5) +
  labs(x = "HIV tree / influenza clade", y = "Estimated value") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + 
  ggtitle("Estimated basic reproduction numbers of HIV trees and WHO influenza clades") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)





df <- data.frame(c(th,tn,tw,t1,t2,t3,t4))
df[,2] <- c("Hunt","Novitsky","Wolf","A1B/135N","A1B/135K","3c3.B","A3")
names(df) <- c("time","dis")
p <- ggplot(data=df, aes(x=dis, y=time)) +
  geom_bar(stat="identity",width=0.5*7/5) +
  geom_text(aes(label=time,family="serif"), vjust=-0.3, size=3.5)+
  labs(x = "HIV tree / influenza clade", y = "Duration (s)") +
  scale_x_discrete(limits=c("Hunt","Novitsky","Wolf","A1B/135N","A1B/135K","3c3.B","A3")) +
  ylim(0,70000) +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Duration of the experiments of estimating basic reproduction numbers") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)




#################################
#################################
load("hivtrees.Rdata")

Hunt <- lapply(Hunt, normalizetree)
Novitsky <- lapply(Novitsky, normalizetree)
Wolf <- lapply(Wolf, normalizetree)

numTrees = 100
numTips = 500

HIVTrees <- list("Hunt" = Hunt,  "Novitsky" = Novitsky, "Wolf" = Wolf)
HIVTrees <- lapply(unlist(HIVTrees, recursive = FALSE, use.names = FALSE), as.phylo)
names(HIVTrees) <- c(rep("Hunt", numTrees), rep("Novitsky", numTrees), rep("Wolf", numTrees))

A <- HIVTrees

#Lat
L <- treeToLattice(A)

DL <- lattToDistMat(L)

KL <- pam(DL,3)
KLL <- KL$clustering

fit <- cmdscale(DL, k = 3)

methodResults <- data.frame(fit, "color" = names(HIVTrees))
p <- ggplot(methodResults, aes(x = X1, y = X2, color = color)) +
  geom_point() +
  scale_color_manual(values=c('#0072BD', '#D95319', '#EDB120'))+
  labs(x = "Dim. 1", y = "Dim. 2", color = "Dataset:") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Visualization of lattice distances between normalized HIV trees") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))


tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)



#Poly
DP <- treeToDistMat(A)

KP <- pam(DP,3)
KPP <- KP$clustering

fit <- cmdscale(DP, k = 3)

methodResults <- data.frame(fit, "color" = names(HIVTrees))
p <- ggplot(methodResults, aes(x = X1, y = X2, color = color)) +
  geom_point() +
  scale_color_manual(values=c('#0072BD', '#D95319', '#EDB120'))+
  labs(x = "Dim. 1", y = "Dim. 2", color = "Dataset:") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Visualization of polynomial distances between HIV tree shapes") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))


tp <- set_panel_size(p,
                     width  = unit(16, "cm"),
                     height = unit(9, "cm"))

grid.draw(tp)




