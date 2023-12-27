library("igraph")
library(ggplot2)
library(igraphdata)

set.seed(42)

############################# Spectral clustering algorithms #################################
# This section focuses on creating functions for normalized and absolute spectral clustering 
# and testing them on various simple graph structures

## Part 1: Implement Normalized and Absolute Spectral Clustering

# Function for normalizing a vector
normalize <- function(x) {
  return(x / sqrt(sum(x^2)))
}

spectral_clustering <- function(G, K = 3) {
  n <- vcount(G)
  Lap <- laplacian_matrix(G, norm = TRUE, sparse = FALSE)
  espec <- eigen(Lap)
  
  # Plot eigenvalues
  df_plot <- data.frame(index = 1:n, eigenvalue = espec$values[n:1])
  ggplot(data = df_plot, aes(x = index, y = eigenvalue)) +
    geom_point(col = "red") +
    labs(title = "Eigenvalues", ylab = 'eigenvalues', xlab = 'indices')
  
  eigvec <- espec$vectors[, n:(n - K + 1)]
  eigvec <- t(apply(eigvec, 1, normalize))
  
  res <- kmeans(eigvec, centers = K, nstart = 15)
  plot(G, vertex.color = res$cluster)
  return(res$cluster)
}


# Function to perform absolute spectral clustering and plot clusters
spectral_clustering_abs <- function(G, K = 3) {
  n <- vcount(G)
  Lap <- laplacian_matrix(G, norm = TRUE, sparse = FALSE)
  Labs <- diag(1, n) - Lap
  espec <- eigen(Labs)
  
  # Plot absolute eigenvalues
  df_plot <- data.frame(index = 1:n, eigenvalue = abs(espec$values[n:1]))
  ggplot(data = df_plot, aes(x = index, y = eigenvalue)) +
    geom_point() +
    labs(title = "Absolute Eigenvalues of L_abs", ylab = 'absolute eigenvalues', xlab = 'indices')
  index <- order(abs(espec$values), decreasing = TRUE)[1:K]
  eigvec <- espec$vectors[, index]
  res <- kmeans(eigvec, centers = K, nstart = 15)
  plot(G, vertex.color = res$cluster)
  return(res$cluster)
}


# test 
clusters <- 3
G1 <- sample_gnp(3, 0.85)
G2 <- sample_gnp(4, 0.7)
G3 <- sample_gnp(5, 0.8)
G <- G1 + G2 + G3 
plot(G)

spectral_clustering(G,clusters) 
spectral_clustering_abs(G,clusters)

G <- erdos.renyi.game(20, 0.15, directed = FALSE)

spectral_clustering(G,clusters) 
spectral_clustering_abs(G,clusters)


## Part 2: Analyzing Spectral Clusters on Various Graphs

n <- 30 
pi <- c(0.3, 0.5, 0.2)  # Group proportions
group_sizes <- n * pi
connectivity_matrix <- matrix(c(0.8, 0.2, 0.2,
                                0.2, 0.9, 0.15,
                                0.2, 0.15, 0.85), nrow = clusters)  # Matrix of connectivities per group pairs
G <- sample_sbm(n, pref.matrix = connectivity_matrix, block.sizes = group_sizes)
spectral_clustering(G, clusters)
spectral_clustering_abs(G, clusters)

n <- 30
p <- 0.2
G <- sample_gnp(n, p)
spectral_clustering(G, clusters)
spectral_clustering_abs(G, clusters)

n1<- 10
n2<- 15
G <- sample_bipartite(n1,n2,p=0.7)
spectral_clustering(G, clusters)
spectral_clustering_abs(G, clusters)


clusters <-2 
star1 <- make_star(10, mode = "undirected")
star2 <- make_star(10, mode = "undirected")
G_star <- star1 + star2
spectral_clustering(G_star, clusters)
spectral_clustering_abs(G_star, clusters)

G_stars <- make_star(8, mode = "undirected")
G_stars <- G_stars + make_star(8, mode = "undirected")
G_stars <- add_edges(G_stars, c(8, 16)) 
spectral_clustering(G_stars, clusters)
spectral_clustering_abs(G_stars, clusters)



# Part 3: Practical Application
# Retrieve the "karate" dataset already present in the library
data(karate)

# Extract edge data from the graph
data <- as_data_frame(karate, "edges")

# Create a graph from the data frame
G <- graph_from_data_frame(data, directed = FALSE)  # The original graph is undirected

# Upgrade the graph for compatibility 
G <- upgrade_graph(G)

tab <-spectral_clustering(G, clusters)
tabulate(tab)

res <- spectral_clustering_abs(G, clusters)
tabulate(tab)