library("igraph")
library("ggplot2")
library(igraphdata)


# In the G(n,p) model:
  
  # The estimator p^ for p is calculated as the ratio of the number of edges
  # |E| to n(n-1)/2 for undirected graphs, and to n(n-1) for directed graphs.
  # The marginal distribution of the node's degree Di in G(n,p) follows 
  # a binomial distribution Binom(n-1,p)
  
# Define the value of n
n <- 1000

# Define different values of p
p_values <- c(1/n, 3/n, 5/n^2)

# Create a data frame to store the results
results <- data.frame(n = numeric(0), p = numeric(0), rel_error = numeric(0), explanation = character(0))

# Loop through different values of p
for (p in p_values) {
  G <- sample_gnp(n, p)  # Generate the Erdös-Rényi graph
  
  # Calculate the number of edges and vertices in the graph
  E <- ecount(G)
  V <- vcount(G)
  
  # Calculate the estimator for p and the relative error
  p_hat <- 2 * E / V / (V - 1)
  rel_error <- abs(p_hat - p) / p
  
  # Explain the graph behavior based on the values of p * n
  if (n * p < 1) {
    explanation <- "If np < 1, then the graph certainly does not have connected components larger than O(log(n))."
  } else if (n * p == 1) {
    explanation <- "If np = 1, then the graph certainly has a larger component whose size is of the order n^(2/3)."
  } else {
    explanation <- "If np > 1, the graph certainly has a single giant component that contains a positive fraction of the vertices. No other component will have more than O(log(n)) vertices."
  }
  
  # Store the results in the data frame
  results <- rbind(results, data.frame(n = rep(n, 1), p = rep(p, 1), rel_error = rep(rel_error, 1), explanation = rep(explanation, 1)))
  
  while (ecount(G) < 100) {
    new_graph <- sample_gnp(n, p)  # Generate a new graph
    G <- graph.union(G, new_graph)  # Union the new graph with the existing graph
    E <- ecount(G)
    V <- vcount(G)
  }
  # Calculating empirical and theoretical degree distributions
  degree.df <- data.frame(deg = 0:max(degree(G)), freq = degree_distribution(G))  # Empirical degree distribution
  deg_plot <- ggplot(degree.df) + aes(x = deg, y = freq) + geom_bar(stat = "identity", fill = "blue") + ggtitle("Degree Distribution - Erdös-Rényi Graph") 
  
  # Overlapping theoretical Binomial distribution
  deg_plot_binomial <- deg_plot + geom_line(aes(y = dbinom(deg, size = V - 1, prob = p)), colour = "red") + labs(subtitle = "Empirical vs. Binomial")
  
  # Overlapping Poisson approximation
  deg_plot_poisson <- deg_plot + geom_line(aes(y = dpois(deg, lambda = (V - 1) * p)), colour = "red") + labs(subtitle = "Empirical vs. Poisson")
  
  # Show the plots
  print(deg_plot_binomial)
  print(deg_plot_poisson)
  
  # The theoretical curves, represented by Binomial and Poisson distributions, overlay well with the empirical barplots.
  # This alignment indicates the similarity between the expected and observed node degree distributions.
  # Additionally, the absence of high-degree nodes suggests a relatively homogenous network in terms of node connections.
  
}

# Create a plot of relative error for different values of p
ggplot(results, aes(x = p, y = rel_error, color = as.factor(n))) +
  geom_line() +
  scale_x_log10() +
  ylab("Relative Error (|p^ - p| / p)") +
  xlab("p") +
  ggtitle("Relative Error vs. p for Different Values of n")

print(results)
