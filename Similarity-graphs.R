library("igraph")
library("ggplot2")
library(igraphdata)
library(mlbench)

############################# Similarity graphs #################################

# Setting up the environment and generating a two-spiral dataset
n <- 100
simu <- mlbench.spirals(100, 1, 0.025)
data <- simu$x

# Applying k-means clustering to the dataset for two clusters
kmeans_result <- kmeans(data, centers = 2)

# Visualizing the k-means clustering result on the spiral dataset
plot(data, col = kmeans_result$cluster)

# Defining a function to calculate Gaussian similarity between two data points
gauss <- function(x1, x2, sigma = 1){   
  return( exp(-sum((x1 - x2)^2) /(2 * sigma^2)) )
}

# Constructing a Gaussian similarity matrix for the dataset
sim_gauss <- function(data) {
  n <- nrow(data)
  S <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      S[i, j] <- gauss(data[i, ], data[j, ])
      S[j, i] <- S[i, j] 
    }
  }
  return(S)
}

# Creating a dense graph from the similarity matrix
G_dense <- sim_gauss(data)
G_spiral <- graph_from_adjacency_matrix(G_dense, mode = "undirected", weighted = TRUE)
plot(G_spiral)

# Applying spectral clustering to the dense graph and visualizing the result
spec <- spectral_clustering(G_spiral,clusters)
plot(data, col = spec)

# Computing the normalized Laplacian matrix and its eigenvalues and eigenvectors
lap <- laplacian_matrix(G_spiral, norm = TRUE, sparse = FALSE)
espec <- eigen(lap)
eigvec <- espec$vectors[, n:(n - clusters + 1)]
eigvec <- t(apply(eigvec, 1,normalize))

# Separating the data points based on original classes and plotting them in a transformed space
class1_x <- eigvec[simu$classes == 1, 1]
class1_y <- eigvec[simu$classes == 1, 2]
class2_x <- eigvec[simu$classes == 2, 1]
class2_y <- eigvec[simu$classes == 2, 2]
x_range <- range(eigvec[, 1])
y_range <- range(eigvec[, 2])
plot(class1_x, class1_y, col = 'red', xlim = x_range, ylim = y_range)
points(class2_x, class2_y, col = 'blue')

# Applying absolute spectral clustering and visualizing the result
spec_abs <- spectral_clustering_abs(G_spiral, clusters)
plot(data, col = spec_abs)

# Preparing for similarity graph construction with different epsilon thresholds
svec <- G_dense[upper.tri(G_dense, diag = FALSE)]
eps <- quantile(svec, probs = 0.75)
eps75 <- G_dense
eps75[G_dense < eps] <- 0
eps75[G_dense >= eps] <- 1
Geps75 <- graph_from_adjacency_matrix(eps75, mode = 'undirected')

# Applying normalized and absolute spectral clustering to the ε-neighborhood graph (75% quantile)
spec <- spectral_clustering(Geps75, 2)
plot(data, col = spec)
spec_abs <- spectral_clustering_abs(Geps75, 2)
plot(data, col = spec_abs)

# Repeating the process for a higher epsilon threshold (95% quantile)
eps <- quantile(svec, probs = 0.95)
G_eps95 <- G_dense
G_eps95[G_dense < eps] <- 0
G_eps95[G_dense >= eps] <- 1
Geps95 <- graph_from_adjacency_matrix(G_eps95, mode = 'undirected')
spec <- spectral_clustering(Geps95, 2)
plot(data, col = spec)
spec_abs <- spectral_clustering(Geps95, 2)
plot(data, col = spec_abs)

# Constructing valued graphs based on mutual and simple nearest neighbors
p <- 2 * floor(log(n))
G_dir <- G_dense
for (i in 1:n) {
  ind <- order(G_dense[i, ], decreasing = TRUE)
  G_dir[i, ind[(p + 1):n]] <- 0
}
Gm <- G_dir
Gm[Gm != t(Gm)] <- 0
Gm_graph <- graph_from_adjacency_matrix(Gm, mode = 'undirected', weighted = TRUE)
plot(Gm_graph)

# Applying spectral clustering to the mutual neighbors graph
cl_normm <- spectral_clustering(Gm_graph, 2)
plot(data, col = cl_normm)
cl_absm <- spectral_clustering_abs(Gm_graph, 2)
plot(data, col = cl_absm)

# Repeating the process for simple nearest neighbors
p <- 2
G_dir <- G_dense
for (i in 1:n) {
  ind <- order(G_dense[i, ], decreasing = TRUE)
  G_dir[i, ind[(p + 1):n]] <- 0
}
Gm <- G_dir
Gm[Gm != t(Gm)] <- 0
Gm_graph <- graph_from_adjacency_matrix(Gm, mode = 'undirected', weighted = TRUE)
cl_abs2 <- spectral_clustering_abs(Gm_graph, 2)
plot_data2 <- plot(data, col = cl_abs2)
cl_abs2_4 <- spectral_clustering_abs(Gm_graph, 4)
plot_data2_4 <- plot(data, col = cl_abs2_4)
