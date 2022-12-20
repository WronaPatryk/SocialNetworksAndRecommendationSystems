
#install.packages("igraph")

library(igraph)


A1 <- as.matrix(read.csv("competition/D1-K=2.csv", header = FALSE))
A2 <- as.matrix(read.csv("competition/D2-K=7.csv", header = FALSE))
A3 <- as.matrix(read.csv("competition/D3-K=12.csv", header = FALSE))

class(A1) == "matrix"


A_becomes_directed <- function(A){
  # get the components
  cmpts <- components(graph_from_adjacency_matrix(A, mode = 'undirected'))$membership
  G <- split(1:length(cmpts), cmpts)
  numberOfComponents <- length(G)
  
  # verify if the graph is connected
  if ( numberOfComponents > 1){
    # if it is not connected, makes it connected:
    # (!) note that there can be numerous 'connecting' techniques
    # -> I take THE LAST vertex of each component 
    # -> and connect it to THE FIRST vertex of the next component
    for(m in 1:(numberOfComponents - 1)){
      A[[G[[m]][length(G[[m]])], G[[m+1]][1] ]] <- 1
      A[[G[[m+1]][1], G[[m]][length(G[[m]])] ]] <- 1
    }
  }
  return(A)
}

A1_directed <- A_becomes_directed(A1)


# calculates the laplasian of the connected graph
# then returns the matrix of its k eigen vectors as columns
laplacian_eigen <- function(A,k){
  # number of observations should be at least the number of communities
  stopifnot(k < nrow(A) && k > 1)
  # the formula is:  L = D - A, 
  # D is a diagonal matrix such that D(i,i) = degree of i-th vertex
  D <- degree(graph_from_adjacency_matrix(A))
  A <- - A 
  for (i in 1:nrow(A)){
    A[i,i] <- D[i]
  }
  # by default the eigen values & corresp. vectors are in decreasing order
  ev <- eigen(A)
  vectors <- ev$vectors
  # gets the eigen vectors corresponding to the first k smallest eigen values:
  E <- matrix(NA, nrow = nrow(A), ncol = k)
  for (i in 1:k){
    E[,i] <- vectors[,nrow(A)-i]
  }
  return(E)
}

# now, we run kmeans on the result - EXAMPLE ( A1)
E1 <- laplacian_eigen(A1, 2)
kmeans(E1, centers = 2, nstart = 34*34)$cluster





spectral_clustering <- function(A, k){
  connected_A <- A_becomes_directed(A)
  eigen_matrix <- laplacian_eigen(connected_A, k)
  
  # 100000 is a heuristic - high enough to make kmeans convrge
  return(kmeans(eigen_matrix, 
                centers = k, nstart = 100000/nrow(A))$cluster)
}




# ALL RESULTS:



spectral_clustering(A1, 2)
spectral_clustering(A2, 7)
spectral_clustering(A3, 12)



determineK <- function(A){
  ev <- eigen(A)
  
  min_eigval <- min(abs(ev$values))
  second_min_eigval <- min( abs(ev$values)[abs(ev$values)!=min(abs(ev$values))] )
  difference <- second_min_eigval - min_eigval
  
  
  # heuristic '7' made on D3-K=12.csv file
  k <- length(which(abs(ev$values) < min_eigval + difference*7))
  
  return(k)
}


spectral_clustering_without_labels <- function(A){
  connected_A <- A_becomes_directed(A)
  
  # innovative step - determining k
  
  k <- determineK(connected_A)
  
  
  eigen_matrix <- laplacian_eigen(connected_A, k)
  
  # 100000 is also a heuristic - high enough to make kmeans convrge
  return(kmeans(eigen_matrix, 
                centers = k, nstart = 100000/nrow(A))$cluster)
}

A1_UNC <- as.matrix(read.csv("competition/D1-UNC.csv", header = FALSE))
A2_UNC <- as.matrix(read.csv("competition/D2-UNC.csv", header = FALSE))
A3_UNC <- as.matrix(read.csv("competition/D3-UNC.csv", header = FALSE))

spectral_clustering_without_labels(A1_UNC) # 10
spectral_clustering_without_labels(A2_UNC) # 6
spectral_clustering_without_labels(A1_UNC) # 12
