# install.packages("ggplot2")
# install.packages("ClusterR")
# install.packages('gtools')
# install.packages('fpc')
# install.packages('mvtnorm')
library(ggplot2)
library(cluster)
library(gtools)
library(ClusterR)
library(fpc)
library(mvtnorm)
set.seed(2)             #sets the seed for random number generation.

# Exercise 5: Cluster Evaluation
# In this exercises, we want to implement some cluster validity indices in order to evaluate clustering
# results. For that purpose:
# 1. Create a dataset sampled from Gaussian distributions. Further, create the according gold standard.

# Define the parameters used for our gaussians (sigma and min/max of cluster)
sig = cbind(c(0.07,0),c(0,0.07))
# lets sample the same number of points for each cluster
points_per_cluster = 40
clusters = 9
sigmas <- matrix(c(4,2,2,3), ncol=2)
means <- t(matrix( c(1,1,1,3,1,6,3,1,3,3,3,6,6,1,6,3,6,6), nrow=2))
means  = means * 5
X <- c()
# create distributations and combine them into a matrix X
# use rmvnorm and rbind to form at least 10 clusters


for( i in 1:clusters){
  # center the dataset
  X <- rbind(X,rmvnorm(points_per_cluster, mean=means[i,], sigma=sigmas))
}

#X <- matrix(X, ncol=2)
X <- scale(X)
plot(X)
# Set the gold standard - each row in our matrix should should be assigned a cluster number
golds = c()
for(i in 1:clusters){
  golds <- c(golds, rep(i, points_per_cluster))
}
plot_data = data.frame(
  x = X[,1],
  y = X[,2],
  num = golds
)
  
  # We define a plotting function - that colors points by the gold standard
  plot_dataset <- function(A, gold_standard, title){
    plot_data = data.frame(
      x = A[,1],
      y = A[,2],
      num = gold_standard
    )
    g <- ggplot(data = plot_data, aes(x = x, y = y, color = factor(num))) + ggtitle(title) + guides(color=FALSE) +
      theme_classic() +
      coord_fixed(xlim = c(-3, 3), ylim = c(-3,3)) +
      xlab("x") + ylab("y")+
      geom_point()
    return(g)
  }

plot_dataset(X, golds, "Dataset")


f# Let's test our implementation with kmeans
cluster = #
  cluster_result = cluster$cluster
plot_dataset(X, cluster_result, "Dataset")

# 3. Implement the contingency matrix.

## Built the contingency Matrix - (slides page 30 + 31)
# given a dataset with n objects, C{ C_1, .. C_N } is the clustering obtained by the algorithm
# Gold standard K = { K_1, ... K_N' } 
# The overlap between two clusters C_i, K_i is defined as:
# n_ij = C_i intersection K_i

N <- length(unique(cluster_result))
N_gold <- length(unique(golds))
# Create the matrix of size N x N'
mat = matrix(0, nrow=N, ncol=N_gold) 


# 4. Implement the F-measure using the mapping approach
# For each cluster k_i in K and element a, c_j ist the cluster with the highest number of common elements
# Here: TP - if a in k_i  and  a in c_j  ,  TN 0
# FP if a not in k_i  and  a in c_j    , FP if a in k_i  and  a not in c_j
# Prec = TP / (TP + FP) , Rec = TP / (TP + FN)
# F_beta = ((1+betaÂ²)*(prec*rec)) / (betaÂ²*prec + rec)
beta = 1

total_fm = 0

# 5. Implement the Silhouette Coeficient. Plot the Silhouette Value for various number of clusters
# for your test dataset. Does the correct number of clusters yield to the best Silhouette value?
# cohesion - average within cluster distance of x
# separation - average distance of x to the vectors in other clusters
# silhouette value(x) = (separation(x) - cohesion(x)) / max(separation(x), cohesion(x))

# Calculate the Silhouette value


# Adjustment for Silhouette plot implementation from library cluster
plot_silhouette <- function(X, cluster, title){
  sil = silhouette(cluster, dist(X))
  sv = mean(sil[,3])
  sil = data.frame(
    cluster = sapply(sil[,1], function(x) order(-table(sil[,1]))[x]),
    orig_cluster = sil[,1],
    neighbor = sil[,2],
    sil_width = sil[,3]
  )
  
  sil = sil[with(sil, order(cluster, -sil_width)),]
  ordered_sil = data.frame(
    id = c(nrow(sil):1),
    cluster = sil$cluster,
    orig_cluster = sil$orig_cluster,
    neighbor = sil$neighbor,
    sil_width = sil$sil_width
  )
  g <- ggplot() + ggtitle(title) +
    theme_classic() +
    theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
    ylab("Silhouette Value") + xlab(" ") + 
    geom_bar(data=ordered_sil, aes(id, sil_width, fill=factor(orig_cluster)), width=1 , stat = "identity") + guides(fill=FALSE) + 
    geom_segment(data=ordered_sil, aes(x=0, xend = nrow(X), y= sv, yend = sv), colour="black", size=1, linetype="dashed") + 
    scale_y_continuous(limits = c(-.1,1)) + 
    coord_flip()
  return(g)
}


#Or using existing R functon
plot(silhouette(cluster_result, dist(X)))
plot_silhouette(X,cluster_result,"Silhouette")

# Calculate different values
res = external_validation(cluster_result, golds, summary_stats = T)
# in comparison with own computation total_fm
res
total_fm

# 6. Revisit the exercise "Hierarchical Clustering". Here, you used a gene expression 
# dataset of bone marrow samples and tried different hierarchical clusterings.
# (a) Now, assess the various clusterings with the newly implemented measures

# (b) Which method yields the best result?

# (c) Also use k-means clustering and compare the performance.

# (d) Calculate other validity measure using the R-package ClusterR.