#--------------Creating a network from scratch using Igraph-----------------------------
#code by Nick Fountain-Jones 9 Dec 2019

library(igraph)

#create graph object

#mating season (Fig 2a Hamemede et al 2009) - not weighted in this case
g <- graph(edges=c(1,2, 1,5, 1,19, 1,15, 1,21, 1,16, 1,23, 1,22, 1,8, 2,19, 2,10, 2,26, 2,20, 2,12, 2,8, 3,19, 3,25, 3,16, 3,4
                   3,18, 3,7, 3,14, 3,10, 3,11, 3,24, 3,9, 4,13, 4,9, 4,18, 4,16, 4,7, 4,14, 4,25,
                   5,19,  5,10, 5,8, 5,16, 6,12, 7,25,7,19, 7,9, 7,18, 7,23, 7,12, 7,26, 7,22, 7,16, 
                   7,8, 7,21, 7,11, 7,17, 7,14, 7,10, 7,24 ,n=25, directed=F)
plot(g)
