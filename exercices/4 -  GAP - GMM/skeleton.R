library(ggplot2)
library(cluster)
library(fpc)
library(mvtnorm)
# if import not avalailable use 
# install.packages("fpc") # requires an R version > 3.2 - check it by writing version in the console below - need to upgrade 
# install.packages("mvtnorm")

set.seed(2)             #sets the seed for random number generation.
# Exercise Sheet 4 - GAP statistic & GMM

# GAP statistic

# 1. In this exercises, we want to figure out, how to determine the number of clusters. For that,
# a very simple dataset consisting of 600 two-dimensional points provided in a tab-separated file
# should be used. In order to determine the number of clusters in this dataset we want to employ
# the GAP statistic. The GAP statistic basically looks at the total sum-of-squares of clusterings
# for varying ks. These values are compared to the expected sum-of-squares when clustering an 
# appropriate null reference distribution. Normally, such an null reference dataset is created by
# uniformly distributing the same number of data points in a bounding box aligned at the principal
# components of the original data. We will perform this analysis step by step. Perform the following
# steps:
# (a) Download the dataset from the course website.
# (b) Perform a PCA on the dataset and create a dataset R aligned along the two principal components.

# load the dataset
gap_dataset <- read.table("4 -  GAP - GMM/gap_dataset.dat", sep='\t', header = TRUE)
X <- data.matrix(gap_dataset[2:ncol(gap_dataset)])
row.names(X) = gap_dataset$Object

plot(X, xlab="x", ylab="y")
title("Original Dataset")

# Center the dataset - as before you can use the function scale()
X.center <- scale(X, center = TRUE, scale = FALSE)
plot(X.center, xlab="x", ylab="y")
title("Centered Dataset")

# extract the eigenvalues and eigenvectors - as before you can use eigen()
# first two eigenvectors - as in exercise 2
Q <- t(X.center) %*% X.center * 1/(nrow(X.center)-1)
Q.eigen <- eigen(Q)
eigenValues <- Q.eigen$values
eigenVectors <- Q.eigen$vectors

# We will need the first two eigenvectors for our 
e1 <- eigenVectors[,1]
e2 <- eigenVectors[,2]

# now we can draw the two eigenvectors spanning our space with the function lines
lines(X.center[,1], e1[2]/e1[1]*X.center[,1], col=2, lwd=2)
lines(e2[1]/e2[2] * X.center[,2], X.center[,2], col=2, lwd=2)

# Project our dataset by multiplying with the first two eigenvectors
E2 = cbind(eigenVectors[,1],eigenVectors[,2])
R_aligned = X.center %*% matrix(c(e1, e2), ncol=2, nrow=2)

# plot the transformed data
plot_data = as.data.frame(R_aligned)
plot_data$rowname = row.names(plot_data)
qplot(data=plot_data, x=V1, y=V2)

# (c) For the test statistic, create the null reference datasets N_i by sampling 
# the same number of points as the original dataset uniform at random. 
# # Limit the sample range of each dimensions to:
# # min_ râˆˆR (r.x) â‰¤ x â‰¤ max_ râˆˆR  and  min_ râˆˆR (r.y) â‰¤ y â‰¤ max_ râˆˆR r.y .
min_x = min(R_aligned[,1])
max_x = max(R_aligned[,1])
min_y = min(R_aligned[,2])
max_y = max(R_aligned[,2])
# Projected data bounding box - get 4 points by using min and max values from R_aligned[,1] for X - 
# from R_aligned[,2] for Y
bb <- cbind(
  c(min_x, min_y), 
  c(max_x, min_y),
  c(max_x, max_y), 
  c(min_x, max_y),
  c(min_x, min_y))

plot(R_aligned)

lines(bb[1,], bb[2,],col=2, lwd=2)


# Now randomly sample 600 data poitns in this box using runif
# for the minimal x value we want to sample values between min(R_aligned[,1]) and max(R_aligned[,1])
#See ?runif for parameter order
n = nrow(R_aligned)
r.x <- runif(n, min_x, max_x)
r.y <- runif(n, min_y, max_y)

R <- cbind(r.x, r.y)

plot(R)
title("Null-Hypthesis Data")
lines(bb[1,], bb[2,] , col=2, lwd=2)

# Use the inverse transformation from the pca to project the datapoints and the bounding box
# this will lead to the exact same orientation as in the original dataset
R.projected <- R  %*% solve(matrix(c(e1, e2), ncol=2, nrow=2))
bb.projected <- t(bb) %*%  solve(matrix(c(e1, e2), ncol=2, nrow=2)) 
# Plot the projected matrix and the bounding box
plot(R.projected)
lines(bb.projected[,1], bb.projected[,2] , col=2, lwd=2)


# (d) Calculate the total sum-of-squares W_k for varying 1 â‰¤ k â‰¤ 10 
# by applying k-means to the original dataset. 
# Furthermore, calculate the expected total sum-of-squares E[W_k] 
# by performing the same procedure using the null reference datasets N_i with 1 < i < 100

# Lets start clustering
# To get the sum of squares error, we apply kmeans for different ks to the original (centered) data using sapply and kmeans
km = kmeans(X.center, 2, nstart = 100)
# kmeans objects have a fitted method which give us access to different evaluation parameters, eg cluster or size, 
# but also the total within-cluster sum of squares error using tot.withinss
km$cluster
km$size
km$totss
km$tot.withinss
# ?kmeans

# now we apply kmeans for the range of ks and save the tot.withinss
orig.x <- 1:10
orig.y <- sapply(orig.x, function(i){return(kmeans(X.center, centers = i)$tot.withinss)})

plot(orig.x, orig.y, type = "l", col=3, lwd=3, xlab="k", ylab="Within Sum of Squares")

null.x <- orig.x
null.y <- sapply(null.x, function(i){return(kmeans(R.projected, centers = i)$tot.withinss)})
lines(orig.x, null.y, type = "l", col=2, lwd=3, lty=2)

# (e) Calculate and plot the gap statistic g(k) = log E[W_k] - log(W_k) . The optimal k is the k maximizing the gap statistic
gap.y <- log(null.y) - log(orig.y)
plot(orig.x, gap.y, type = "l", lwd=3, xlab="number of clusters k", ylab="GAP")

# we want to select the k, which maximizes the GAP statistic
best_k = which.max(gap.y)

# use best k to cluster dataset and color by cluster_id
cl = kmeans(X.center, best_k, nstart=100)
plot(X.center, xlab="x", ylab="y", col=cl$cluster)
title("Clustered Dataset")

# (f ) Extra:
# Use the bone marrow dataset from the last exercise sheet.
# Perform a principal coordinate analysis 
#and perform the GAP statistic on that dataset ... 
# does the GAP statistic report 3 as the ideal number of clusters?





#######################################################
# EXERCISE 2 - EXPECTATION MAXIMIZATION
#######################################################
# Create Our dataset
x.g1 <- rnorm(10, 1, .8)
x.g2 <- rnorm(10, 4, .8)

X <- sort(c(x.g1,x.g2))
hist(X)

# Guess initial Paramters:
min_ = min(X)
max_ = max(X)

mu.1 = runif(1, min_, max_)
mu.2 =  runif(1, min_, max_)
sig.1 = runif(1, .1, 1)
sig.2 = runif(1, .1, 1)
pi.1 = 0.49
pi.2 = 0.51

gauss <- function( nums, mu, sigma ){
  #f(x) = 1/(√(2 π) σ) e^-((x - μ)^2/(2 σ^2))
  ex_p = exp(-(nums - mu)^2 / (2 * sigma^2))
  y = (1 / (sigma * sqrt(2 * pi))) * ex_p
  return(y)
}

## EXPECTATION STEP ##
# Compute the expected classes for all data points of the two Gaussian distributions
y.1 = pi.1 * gauss(X, mu.1, sig.1)
y.2 = pi.2 * gauss(X, mu.2, sig.2) 
#plot(y.2)

y.1.scaled = y.1 * (1/sum(y.1, y.2))
y.2.scaled = y.2 * (1/sum(y.1, y.2))

# Plot the membership for each point
plot(y.1.scaled)
plot(y.2.scaled)

## MAXIMIZATION STEP ##
# Update mu
den.1 = sum(y.1.scaled)
den.2 = sum(y.2.scaled)

mu.1 = sum(X * y.1.scaled) / den.1
mu.2 = sum(X * y.2.scaled) / den.2

# Update sig
sig.1 = sum( y.1.scaled * ((X - mu.1)^2)) / den.1
sig.2 = sum( y.2.scaled * ((X - mu.2)^2)) / den.2

# Update pi
pi.1 = sum(y.1.scaled) / length(y.1.scaled)
pi.2 = sum(y.2.scaled) / length(y.2.scaled)


# Plot the membership for each point
plot(y.1.scaled)
plot(y.2.scaled)

model <- Mclust(X)
