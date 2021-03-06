#--------------Creating a network from scratch using Igraph-----------------------------
#code by Nick Fountain-Jones 9 Dec 2019

library(igraph)
library(spectralGraphTopology)
library(quadprog)
library(pals)
library(ggplot2)
#create graph object

#mating season (Fig 2a Hamemede et al 2009) - not weighted in this case
gMating <- graph(edges=c(1,201,5, 1,19, 1,15, 1,21, 1,16, 1,23, 1,22, 1,8, 2,19, 2,10, 2,26, 2,20, 2,12, 2,8, 3,19, 3,25, 3,16, 3,4,
                   3,18, 3,7, 3,14, 3,10, 3,11, 3,24, 3,9, 4,13, 4,9, 4,18, 4,16, 4,7, 4,14, 4,25,
                   5,19,  5,10, 5,8, 5,16, 6,12, 7,25, 7,19, 7,9, 7,18, 7,23, 7,12, 7,26, 7,22, 7,16, 
                   7,8, 7,21, 7,11, 7,17, 7,14, 7,10, 7,24,7,20, 7,15, 8,21 ,8,10, 8,23, 8,17, 8,19, 8,20, 8,16,
                   9,18, 9,13, 9,25, 9,21, 9,16, 9,10, 9,15, 10,17, 10,25, 10,20, 10,15, 10,23, 11,14, 11,25, 11,16,
                   11,20, 12,22, 12,26, 12,17, 12,16, 12,20, 13,24, 13,16, 13,14, 13,20, 13,25, 14,25, 14,24, 
                   14,20, 15,17, 15,18, 15,16, 15,20, 16,23, 17,23, 17,20, 17,26, 18,21, 18,19, 18,20, 19,20, 19,25,
                   19,23, 20,25, 20,23, 20,26, 20,22, 21,25, 22,23, 23,26, 26,27),  n=27, directed=F)

layout <- layout_nicely(gMating)

plot(gMating, layout=layout)
str(gMating)

#Degree to compare to Hamede et al
Mating.degrees <- as.data.frame(degree(gMating)); colnames(Mating.degrees) <- 'Degree'
Mating.degrees$Degree <- as.numeric(Mating.degrees$Degree)
summary(Mating.degrees$Degree)#10.5 in Hamede et al. But unclear what encounter length used

# Now, plot it!
ggplot(Mating.degrees, aes(x = Degree )) +
  geom_histogram(breaks=seq(0,30, by=2), col= 'black',fill='grey') + 
  labs(title="Degree Distribution",x="Degree", y="Frequency", y="count") +
  theme_bw()


gAfterMating <- graph(edges=c(1,2, 1,20, 1,10,1,8, 1,15, 1,7,  1,5, 1,23, 2,8, 2,12, 2,10, 2,15, 2,25, 3,9, 3,11, 
                              3,19, 3,4, 3,14, 3,10, 3,8, 3,15, 4,25, 4,21, 4,22, 4,18, 4,7, 4,24, 4,13, 4,19,
                              5,23, 5,6, 5,20, 5,17, 6,17, 6,27,7,15, 7,12, 7,10, 8, 10, 8,13, 8,20, 8,15, 8,12, 8,18,
                              9,13, 9,25, 9,15, 9,14, 9,16, 10, 15, 10,20, 10,19, 10,13, 10,12, 10,14, 11,13, 11,14, 
                              12,23, 12,20, 12,22, 13,24, 13,14, 13,25, 13,18, 14,24, 14,23, 14,17, 14,22, 16,20, 16,25, 
                              17,22, 17,27, 17,24, 17,23, 17,20, 17,25, 18,20, 18,21, 18,24,  19,26, 19,25, 19,20, 20,21, 
                              20,23, 21,24,20,26, 21,23, 21,22, 22,26, 22,24, 24,25), n=27, directed=F)
plot( gAfterMating )
str( gAfterMating )

#Degree to compare to Hamede et al
AfterMating.degrees <- as.data.frame(degree(gAfterMating)); colnames(AfterMating.degrees) <- 'Degree'
AfterMating.degrees$Degree <- as.numeric(AfterMating.degrees$Degree)
summary(AfterMating.degrees$Degree)#7.19 in the Hamede et al

# Now, plot it! 
ggplot(AfterMating.degrees, aes(x = Degree )) +
  geom_histogram(breaks=seq(0,30, by=2), col= 'black',fill='grey') + 
  labs(title="Degree Distribution",x="Degree", y="Frequency", y="count") +
  theme_bw()

#After Mating analysis

#convert to a laplacian matrix
M <- laplacian_matrix(gAfterMating,sparse = FALSE)
#not the most efficent way but fine for the networks we are looking at. ARPACK eiganvalue algorithm may be better
#cal eigenvectors
spec <- eigen(M)
#plot all values
plot(spec$values, col="blue")
#take the second smallest eigenvalue
FiedlerValue <- spec$values[length(spec$value)-1]
#if  eigenvalue is greater than 0 graph is connected
FiedlerVector <- spec$vectors[,length(spec$values)-1]
FiedlerVectorA <- spec$vectors[,length(spec$values)-2]
plot(FiedlerVector)

combinedA <- as.data.frame(cbind(FiedlerVector,FiedlerVectorA))

#what is cluster membership
clust <- kmeans(FiedlerVector,2)
clust[1]#cluster membership  across nodes

ggplot(combinedA, aes(FiedlerVector,FiedlerVectorA ))+
  geom_point(aes(colour= factor(kmeans(combinedA, centers=3,nstart=999)$cluster)), size=3)+ 
  geom_text(label=rownames(combinedA), size=3,nudge_x = 0.02)+
  theme_bw()+
  scale_colour_brewer('My groups', palette = 'Set2')


#Mating season analysis

#convert to a laplacian matrix
M1 <- laplacian_matrix(gMating,sparse = FALSE)
#not the most efficent way but fine for the networks we are looking at. ARPACK eiganvalue algorithm may be better
#cal eigenvectors
spec1 <- eigen(M1)
#plot all values
plot(spec1$values, col="blue")
#take the second smallest eigenvalue
FiedlerValue1 <- spec1$values[length(spec1$value)-1]
FiedlerValue1 

#if  eigenvalue is greater than 0 graph is connected
FiedlerVector1 <- spec1$vectors[,length(spec1$values)-1]
FiedlerVector1a <- spec1$vectors[,length(spec1$values)-2]

row.names(combined) <- seq(1:27)

combined <- as.data.frame(cbind(FiedlerVector1,FiedlerVector1a))

clust1 <- kmeans(combined,2, nstart=999)
clust1[1]#cluster membership  across nodes

#what is cluster membership
clust1 <- kmeans(combined,3, nstart=999)
#3 looks right

ggplot(combined, aes(FiedlerVector1,FiedlerVector1a ))+
  geom_point(aes(colour= factor(kmeans(combined, centers=3,nstart=999)$cluster)), size=3)+ 
  geom_text(label=rownames(combined), size=3,nudge_x = 0.02)+
  theme_bw()+
  scale_colour_brewer('My groups', palette = 'Set2')






#--------------------------------------------------------------------
#scg
#--------------------------------------------------------------------

## SCG of real-world network


n <- vcount(gMating)
#interv <- c(100,100,50,25,12,6,3,2,2)
cg <- scg(gMating, ev= n-(1:2),nt=2, mtype="laplacian",
          algo="exact_scg", epairs=TRUE) #ev gives eigenpairs to be oreserved. Dont get the algo bit. nt isn't clear either (number of groups)

## are the eigenvalues well-preserved?
gt <- cg$Xt # this is the spectral graph
nt <- vcount(gt)
Lt <- laplacian_matrix(gt)
evalt <- eigen(Lt, only.values=TRUE)$values[nt-(1:2)]
res <- cbind(cg$values, evalt)
res <- round(res,5) #rounds number
colnames(res) <- c("lambda_i","lambda_tilde_i")
rownames(res) <- c("N-1","N-2")
print(res)

## use SCG to get the communities
com <- scg(laplacian_matrix(gMating), ev=n-c(1,2), nt=2)$groups
col <- rainbow(max(com))
layout <- layout_nicely(gMating)

#plot untransformed network with scg groups
plot(gMating, layout=layout, vertex.size=3, vertex.color=col[com])

## display the coarse-grained graph
gt <- simplify(as.undirected(gt))
#layout.cg <- layout_with_kk(gt) #The Kamada-Kawai layout algorithm
com.cg <- scg(laplacian_matrix(gt), nt-c(1,2), 2)$groups
vsize <- sqrt(as.vector(table(cg$groups)))
#side by side
op <- par(mfrow=c(1,2))
plot(gMating, layout=layout, vertex.size=3, vertex.color=col[com])
plot(gt, layout=layout, vertex.size=3, 
     vertex.color=col[com.cg])
par(op)
#-------------------------

#------------other code------------------------

#calaculate laplacian matrix
embed <- embed_laplacian_matrix(gMating, 5)
#D are the eiganvalues
embed$D
embed$D[length(embed$D)-1]
#spectralTopology package

adjMat <- as.matrix(as_adjacency_matrix(gMating, type = c("both")))
p <- ncol(adjMat)
graph <- learn_k_component_graph(cov(adjMat), w0 = "qp", beta = 10, k = 2,verbose = FALSE)
graph$Laplacian

ei <- eigen(graph$Laplacian)
ei$values[length(ei$value)-1] #gives second smallest eigenvalue of matrix

str(graph)



# plots
net <- graph_from_adjacency_matrix(graph$Laplacian, mode = "undirected", weighted = TRUE)
colors <- brewer.blues(100)
c_scale <- colorRamp(colors)
E(net)$color = apply(c_scale(abs(E(net)$weight) / max(abs(E(net)$weight))), 1,
                     function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(net)$color = "pink"
par(bg = NA)
plot(net, vertex.label = names,
     vertex.size = 3,
     vertex.label.dist = 1,
     vertex.label.family = "Helvetica",
     vertex.label.cex = .8,
     vertex.label.color = "black")


#other stuff of use from Igraph


#spec <- spectrum(gMating)[c("values", "vectors")]
#str( spec )
gp <- scg(gMating,1,10, algo = c("optimum"))
str(gp)
layout <- layout_with_kk(g)
plot(gp$Xt, layout = layout_with_kk)

