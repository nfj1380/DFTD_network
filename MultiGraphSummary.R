#-------------------------------------------
# Calculating network properties for multiple graphs
#-------------------------------------------
#code by Nick Fountain-Jones May 2020

library(igraph)
library(EpiModel)
library(intergraph)

#create graph object
folder_name = c('Networks')

Network_sum <- function(folder_name){
  
  networks <- list.files(folder_name)
  Global_summary <- data.frame()
  
    for (j in 1:length(networks))
    {
      g <- read_graph(networks[j], format = c("edgelist"))
      g_names <- gsub(".edges","",networks[j])
      

#remove disconnected nodes
    iso <- which(degree(g)==0) #disconencted networks are going to be an issue - very common. Atleast we can easily remove totally disconencted nodes.
    g2 <- delete.vertices(g, iso)
    
    
    #SIR model using inbuild igraph functionality doesn't work on unconnected graphs. ERGMS based models are better
    #epi model
    #intergraph to convert between igraph and statnet

    g3 <- asNetwork(g2)
 
    formation <- ~edges #very basic ERGM at this stage
    target.stats <- 40 #number of nodes - this may help standardise things a bit
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
    coef.diss
    plot(coef.diss)
    #run the ERGM
    
    est <- netest(g3 , formation, target.stats, coef.diss, edapprox = TRUE) #false makes it a STERGM
   
    dx <- netdx(est, nsims = 5, nsteps = 100,
                nwstats.formula = ~edges)
    plot(dx)
    # Epidemic model (SIR)
    param <- param.net(inf.prob = 0.1, act.rate = 5, rec.rate = 0.02)
    init <- init.net(i.num = 1, r.num = 0)
    control <- control.net(type = "SIR", nsteps = 20, nsims = 100, verbose.int = 0) #R0 = 5
    mod1 <- netsim(est, param, init, control)
    plot(mod1)
    plot(mod1, type = "network", at = 20, col.status = TRUE,
         main = "Prevalence at t20")
    
    d <- as.data.frame(mod1)
    avg_inf <- mean(d$i.num)
    sd_inf <- sd(d$i.num)
    
    #convert to a laplacian matrix
    M <- laplacian_matrix(g2,sparse = FALSE)
    
    #cal eigenvectors
    spec <- eigen(M)
    #take the second smallest eigenvalue
    FiedlerValue<- spec$values[length(spec$value)-1]
    
    Global_summary <- rbind(FiedlerValue,Global_summary )
    row.names(Global_summary) <-  g_names
  
    }
  
  colnames(Global_summary) <- c("Fiedler")
  return(Global_summary)
}
#asian elephant network
a <- Network_sum('Networks')



#convert to a laplacian matrix
M <- laplacian_matrix(g3,sparse = FALSE)

#cal eigenvectors
spec <- eigen(M)
#plot all values
plot(spec$values, col="blue")
#take the second smallest eigenvalue
FiedlerValue<- spec$values[length(spec$value)-1]
FiedlerValue
