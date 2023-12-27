library("igraph")
library("ggplot2")
library(igraphdata)

################################ Configuration models #################################################

# Fit a power-law to the degree distribution of the karate graph
power_law_fit <- fit_power_law(degree(karate))

# Create a dataframe with degree and frequency (excluding zero degrees)
degree_distribution_df <- data.frame(degree = 1:max(degree(karate)), frequency = degree_distribution(karate)[-1])

# Plot the degree distribution (excluding zero degrees) as a histogram
degree_distribution_plot <- ggplot(degree_distribution_df) +
  aes(x = degree, y = frequency) +
  geom_bar(stat = "identity", fill = "blue") +
  ggtitle("Degree Distribution (Excluding Zero Degrees) - Karate Graph")

# Define a function to calculate power-law density for the range of degrees
power_law_density <- function(x, alpha, n) {
  c <- sum((1:(n - 1))^(-alpha))  # Normalization constant for a sum of 1
  return(x^(-alpha)/c)
}

# Superimpose the theoretical power-law distribution using the estimated alpha
degree_distribution_plot +
  geom_line(data = data.frame(degree = 1:(vcount(karate) - 1)),
            aes(x = degree, y = power_law_density(degree, alpha = power_law_fit$alpha, n = vcount(karate))),
            color = "red")

# Function to calculate the expected degree under the Random Degree (RD) model
calculate_esp <- function(data, m) {
  
  # Calculate the degree of each node in the graph
  di <- degree(karate)
  
  # Use outer product to calculate the cross products d_id_j
  pij <- di %*% t(di)
  
  # Choose the constant C based on the parameter m
  if (m == 0)
    C <- max(pij)  # If m=0, C is chosen as the maximum of the degree products
  else  
    C <- 2 * E  # If m=1, C is chosen as twice the number of edges
  
  # Calculate the expected degree under the RD model
  esp_di <- di * (2 * E - di) / C
  
  # Create a dataframe for the plot
  plot_df <- data.frame(di = di, esp_di = esp_di)
  
  # Create the scatter plot of observed and expected degrees
  ggplot(data = plot_df) + aes(x = di, y = esp_di) + geom_point(col = "red") + geom_abline(intercept = 0, slope = 1, color = "red") + labs(x = "Observed Degrees", y = "Mean Degrees under the RD Model")
}

# Call the function for m=0 and m=1
calculate_esp(data,m=0)
calculate_esp(data,m=1)

################################ Tests #################################################

matching_algorithm <- function(degrees){
  # Initialization
  n <- length(degrees)
  total_edges <- sum(degrees)
  edges <- list()
  nodes <- c()
  is_loop <- TRUE
  is_duplicate <- TRUE
  
  # Create extended node list
  for (i in 1:n){
    nodes <- c(nodes, rep(i, degrees[i]))
  }
  
  # Continue until a valid graph is obtained
  while (is_loop || is_duplicate){
    # Generate edges
    edge_index <- 0
    while (total_edges >= 1){
      edge_index <- edge_index + 1
      random_edge <- sample(1:total_edges, 2)
      random_edge <- sort(random_edge)   # Sort the edge indices
      
      if (nodes[random_edge[1]] == nodes[random_edge[2]]) {  # Check for loops
        is_loop <- TRUE
        break
      } else {
        is_loop <- FALSE
        edges[[edge_index]] <- nodes[random_edge]
        nodes <- nodes[-random_edge]
        total_edges <- total_edges - 2
      }
    }
    
    # Check for duplicate edges
    is_duplicate <- as.logical(1 - (length(unique(edges)) == edge_index))
  }
  
  return(edges)
}

# Example 
node_degrees <- c(1, 2, 1)

# Generate edges based on node degrees
generated_edges <- matching_algorithm(node_degrees)

# Display the edges
cat("Edges:", "\n")
for (i in 1:length(generated_edges)) {
  cat(generated_edges[[i]][1], "-", generated_edges[[i]][2], "\n")
}


rewire_graph <- function(Edge_List, Num_iter = 100 * dim(Edge_List)[1]){
  num_edges <- dim(Edge_List)[1]
  Edge_List <- t(apply(Edge_List, 1, sort))
  
  for (i in 1:Num_iter){
    random_edges <- sample(1:num_edges, 2)
    edge1 <- Edge_List[random_edges[1], ]
    edge2 <- Edge_List[random_edges[2], ]
    
    new_edge1 <- sort(c(edge1[1], edge2[2]))
    new_edge2 <- sort(c(edge1[2], edge2[1]))
    
    # Check for 4 different nodes and that new edges don't already exist
    uninodes <- unique(c(edge1, edge2))
    uniedges <- unique(rbind(Edge_List[-c(random_edges[1], random_edges[2]), ], new_edge1, new_edge2))
    
    if ((length(uninodes) == 4) && (length(uniedges) == num_edges)) {
      Edge_List[random_edges[1], ] <- new_edge1
      Edge_List[random_edges[2], ] <- new_edge2
    }
  }
  return(Edge_List)
}

# Example 
edges <- matrix(c(1, 2, 2, 3, 3, 4, 4, 1), ncol = 2, byrow = TRUE)

# Generate rewired edges
rewired_edges <- rewire_graph(edges)

# Display the rewired edges
cat("Rewired Edges:", "\n")
rewired_edges


degseq <- degree(karate)
R <- 1000  # Replications
T <- numeric(R)
T[] <- NA

for (i in 1:R) {
  G <- sample_degseq(degseq, method = "vl")
  T[i] <- sum(count_triangles(G)) / 3
}

obs_tri_count <- sum(count_triangles(karate)) / 3

ggplot(data.frame(count = T)) + 
  aes(count) + 
  geom_histogram() + 
  geom_vline(xintercept = obs_tri_count, col = "red",  linetype = "dashed")

