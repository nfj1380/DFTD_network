## SCG of real-world network
library(igraphdata)
data(immuno)
summary(immuno)
n <- vcount(immuno)
interv <- c(100,100,50,25,12,6,3,2,2)
cg <- scg(immuno, ev= n-(1:9), nt=interv, mtype="laplacian",
          algo="interv", epairs=TRUE)

## are the eigenvalues well-preserved?
gt <- cg$Xt
nt <- vcount(gt)
Lt <- laplacian_matrix(gt)
evalt <- eigen(Lt, only.values=TRUE)#$values[nt-(1:9)] #why only 1:9
res <- cbind(interv, cg$values, evalt)
res <- round(res,5)
colnames(res) <- c("interv","lambda_i","lambda_tilde_i")
rownames(res) <- c("N-1","N-2","N-3","N-4","N-5","N-6","N-7","N-8","N-9")
print(res)

## use SCG to get the communities
com <- scg(laplacian_matrix(immuno), ev=n-c(1,2), nt=2)$groups
col <- rainbow(max(com))
layout <- layout_nicely(immuno)

plot(immuno, layout=layout, vertex.size=3, vertex.color=col[com],
     vertex.label=NA)

## display the coarse-grained graph
gt <- simplify(as.undirected(gt))
layout.cg <- layout_with_kk(gt)
com.cg <- scg(laplacian_matrix(gt), nt-c(1,2), 2)$groups
vsize <- sqrt(as.vector(table(cg$groups)))

op <- par(mfrow=c(1,2))
plot(immuno, layout=layout, vertex.size=3, vertex.color=col[com],
     vertex.label=NA)
plot(gt, layout=layout.cg, vertex.size=15*vsize/max(vsize), 
     vertex.color=col[com.cg],vertex.label=NA)
par(op)
