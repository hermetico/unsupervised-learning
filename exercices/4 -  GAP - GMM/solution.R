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
gap_dataset <- read.table("./gap_dataset.dat", sep='\t', header = TRUE)
X <- data.matrix(gap_dataset[2:ncol(gap_dataset)])
row.names(X) = gap_dataset$Object

plot(X, xlab="x", ylab="y")
title("Original Dataset")

# Center the dataset - as before you can use the function scale()
X.center <- (scale((X), center=TRUE, scale=FALSE))
plot(X.center, xlab="x", ylab="y")
title("Centered Dataset")

# extract the eigenvalues and eigenvectors - as before you can use eigen()
# first two eigenvectors - as in exercise 2
Q <- 1/(nrow(X.center)-1)*((t(X.center) %*% X.center))
eigenValues <- eigen(Q)$values
eigenVectors <- eigen(Q)$vectors

# We will need the first two eigenvectors for our 
e1 <- eigenVectors[,1]
e2 <- eigenVectors[,2]

# now we can draw the two eigenvectors spanning our space with the function lines
lines(X.center[,1], e1[2]/e1[1]*X.center[,1], col=2, lwd=2)
lines(e2[1]/e2[2] * X.center[,2], X.center[,2], col=2, lwd=2)

# Project our dataset by multiplying with the first two eigenvectors
E2 = cbind(eigenVectors[,1],eigenVectors[,2])
R_aligned = t(t(E2) %*% t(X.center))

# plot the transformed data
plot_data = as.data.frame(R_aligned)
plot_data$rowname = row.names(plot_data)
qplot(data=plot_data, x=V1, y=V2)

# (c) For the test statistic, create the null reference datasets N_i by sampling 
# the same number of points as the original dataset uniform at random. 
# # Limit the sample range of each dimensions to:
# # min_ râˆˆR (r.x) â‰¤ x â‰¤ max_ râˆˆR  and  min_ râˆˆR (r.y) â‰¤ y â‰¤ max_ râˆˆR r.y .

# Projected data bounding box - get 4 points by using min and max values from R_aligned[,1] for X - 
# from R_aligned[,2] for Y
bb <- cbind(c(min(R_aligned[,1]), min(R_aligned[,1]), max(R_aligned[,1]), max(R_aligned[,1]), min(R_aligned[,1])), c(max(R_aligned[,2]), min(R_aligned[,2]), min(R_aligned[,2]), max(R_aligned[,2]), max(R_aligned[,2])))
plot(R_aligned)
lines(bb[,1], bb[,2] , col=2, lwd=2)

# Now randomly sample 600 data poitns in this box using runif
# for the minimal x value we want to sample values between min(R_aligned[,1]) and max(R_aligned[,1])
#See ?runif for parameter order
r.x <- runif(n=nrow(X), min = min(R_aligned[,1]), max = max(R_aligned[,1]))
r.y <- runif(n=nrow(X), min = min(R_aligned[,2]), max = max(R_aligned[,2]))
R <- cbind(r.x, r.y)

plot(R)
title("Null-Hypthesis Data")
lines(bb[,1], bb[,2] , col=2, lwd=2)

# Use the inverse transformation from the pca to project the datapoints and the bounding box
# this will lead to the exact same orientation as in the original dataset
R.projected <- t(E2 %*% t(R))
bb.projected <- t(E2 %*% t(bb))
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
orig.y <- sapply(orig.x, function(k) kmeans(X.center, k, nstart=100)$tot.withinss)
plot(orig.x, orig.y, type = "l", col=3, lwd=3, xlab="k", ylab="Within Sum of Squares")

orig.y
null.x <- 1:10
null.y <- sapply(null.x, function(k) kmeans(R.projected, k, nstart=100)$tot.withinss)
lines(orig.x, null.y, type = "l", col=2, lwd=3, lty=2)

# (e) Calculate and plot the gap statistic g(k) = log E[W_k] âˆ’ log(W_k) . The optimal k is the k maximizing the gap statistic
gap.y <- log(null.y)- log(orig.y)
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

bone <- read.table("./bone_marrow.dat", sep='\t')
X <- data.matrix(bone)
row.names(X) = bone$V1
X.center = scale(X, center = TRUE, scale = FALSE)
B <- t(X.center) %*% X.center
eigenValues <- eigen(B)$values
eigenVectors <- eigen(B)$vectors

E2 = cbind(eigenVectors[,1],eigenVectors[,2])

TD = t(E2) %*% t(X)

plot.data = t(TD)
plot.data = as.data.frame(plot.data)
plot.data$oname = row.names(plot.data)
qplot(data=plot.data, x=V1, y=V2)

ccode = c()
for(i in 1:(nrow(plot.data))){
  tmp <- plot.data[i,]$oname
  if(grepl("B.cell", tmp)){
    ccode = c(ccode, "B-Cell")
  }
  if(grepl("T.cell", tmp)){
    ccode = c(ccode, "T-Cell")
  }
  if(grepl("AML", tmp)){
    ccode = c(ccode, "AML")
  }
}
plot.data$color = ccode
qplot(data=plot.data, x=V1, y=V2, label = oname) + geom_point(aes(colour = color))
# With this matrix, you could repeat the GAP statistic analysis)
X = t(TD)


# Center the dataset
X.center <- (scale((X), center=TRUE, scale=FALSE))
plot(X.center, xlab="x", ylab="y")
title("Centered Dataset")

# extract the eigenvectors - as in exercise 2
Q <- 1/(nrow(X.center)-1)*((t(X.center) %*% X.center))
eigenValues <- eigen(Q)$values
eigenVectors <- eigen(Q)$vectors

e1 <- eigenVectors[,1]
e2 <- eigenVectors[,2]

# now we can draw the two eigenvectors spanning our space
lines(X.center[,1], e1[2]/e1[1]*X.center[,1], col=2, lwd=2)
lines(e2[1]/e2[2] * X.center[,2], X.center[,2], col=2, lwd=2)


# Projection our dataset by multiplying with the first two eigenvectors -> Projection
E2 = cbind(eigenVectors[,1],eigenVectors[,2])
# Calculate Projection
R_aligned = t(t(E2) %*% t(X.center))

# plot the transformed data
plot_data = as.data.frame(R_aligned)
plot_data$rowname = row.names(plot_data)
qplot(data=plot_data, x=V1, y=V2)

# create the null reference datasets 
# Projected data bounding box
bb <- cbind(c(min(R_aligned[,1]), max(R_aligned[,1]), max(R_aligned[,1]), min(R_aligned[,1]), min(R_aligned[,1])), c(max(R_aligned[,2]), max(R_aligned[,2]), min(R_aligned[,2]), min(R_aligned[,2]), max(R_aligned[,2])))
plot(R_aligned)
lines(bb[,1], bb[,2] , col=2, lwd=2)

# Now randomly sample data points in this box using runif
r.x <- runif(n=nrow(X), min = min(R_aligned[,1]), max = max(R_aligned[,1]))
r.y <- runif(n=nrow(X), min = min(R_aligned[,2]), max = max(R_aligned[,2]))
R <- cbind(r.x, r.y)

plot(R)
title("Null-Hypthesis Data")
lines(bb[,1], bb[,2] , col=2, lwd=2)

# Use the inverse transformation from the pca to project the datapoints and the bounding box
R.projected <- t(E2 %*% t(R))
bb.projected <- t(E2 %*% t(bb))
# Plot the projected matrix and the bounding box
plot(R.projected)
lines(bb.projected[,1], bb.projected[,2] , col=2, lwd=2)

# now we apply kmeans with a bigger range for k
orig.x <- 1:20
orig.y <- sapply(orig.x, function(k) kmeans(X.center, k, nstart=100)$tot.withinss)
plot(orig.x, orig.y, type = "l", col=3, lwd=3, xlab="k", ylab="Within Sum of Squares")

orig.y
null.x <- orig.x
null.y <- sapply(null.x, function(k) kmeans(R.projected, k, nstart=100)$tot.withinss)
lines(orig.x, null.y, type = "l", col=2, lwd=3, lty=2)

# Calculate and plot the gap statistic g(k) = log E[W_k] âˆ’ log(W_k) . The optimal k is the k maximizing the gap statistic
gap.y <- log(null.y)- log(orig.y)
plot(orig.x, gap.y, type = "l", lwd=3, xlab="number of clusters k", ylab="GAP")

# we want to select the k, which maximizes the GAP statistic
best_k = which.max(gap.y)
best_k
# 16 - so not 3. Sadly this shows, that the GAP statistic is also not infallible


#######################################################
# EXERCISE 2 - EXPECTATION MAXIMIZATION
#######################################################


# Develop an expectation maximization algorithm fitting Gaussians over a given dataset. For the
# sake of simplicity, we will only use 1-D Gaussians.
# 
# (a) Generate an artificial dataset by sampling random values from different Gaussian distributions.

# Create a new dataset by sampling from two normal distributions with rnorm and combining them to matrix X
x.g1 <- rnorm(250, 1, .8)
x.g2 <- rnorm(150, 3, .8)
X <- sort(c(x.g1,x.g2))
hist(X)

# set parameters mu, sigma describing the gaussian and the mixing coefficient pi for two gaussian models
#  Î¼ is the mean of the distribution and Ïƒ the standard deviation.
mu.1 = 8
mu.2 = 12
sig.1 = 0.8
sig.2 = 0.8
pi.1 = .5
pi.2 = .5


## EXPECTATION STEP ##
# Compute the expected classes for all data points of the two Gaussian distributions
y.1 = pi.1 * dnorm(X, mu.1, sd = sqrt(sig.1), log=F)
y.2 = pi.2 * dnorm(X, mu.2, sd = sqrt(sig.2), log=F)

# Scale y.1 and y.2 so all elements with same index numbers sum to 1
y.1.scaled = y.1 * (1/(y.1+y.2))
y.2.scaled = y.2 * (1/(y.1+y.2))

# Plot the membership for each point
plot(y.1.scaled)
plot(y.2.scaled)



## MAXIMIZATION STEP ##
# Update mu
mu.1 = sum(y.1.scaled * X)/sum(y.1.scaled)
mu.2 = sum(y.2.scaled * X)/sum(y.2.scaled)

# Update sig
sig.1 = sum(y.1.scaled * (X-mu.1)*(X-mu.1))/sum(y.1.scaled)
sig.2 = sum(y.2.scaled * (X-mu.2)*(X-mu.2))/sum(y.2.scaled)

# Update pi
pi.1 = sum(y.1.scaled)/length(X)
pi.2 = sum(y.2.scaled)/length(X)


# Lets build three clusters for our evaluation and concatenate them with with rbind into a single matrix X
# use rmvnorm with different parameters
# see ?rmvnorm
x = 0
sigma <- matrix(c(2,0,0,1), ncol=2)
x <- rmvnorm(n=100, mean=c(5,6), sigma=sigma)
sigma <- matrix(c(1,0,0,1), ncol=2)
x <- rbind(x,rmvnorm(n=200, mean=c(-3,4), sigma=sigma))
sigma <- matrix(c(4,2,2,3), ncol=2)
x <- rbind(x,rmvnorm(n=100, mean=c(0,-8), sigma=sigma))
plot(x)

# Initialization
# lets assume we have k = 4 clusters 
k = 4
# then our initial estimation of pi would be (1/4, 1/4, 1/4, 1/4) 
pi = rep(1/k,k)
# and our initial mu
mu = matrix(c(runif(2*k, min=-5,max=5)), nrow = k)
sigma = array(0, c(k, 2, 2));
for(i in 1:k){
  sigma[i,,] = c(1,0,0,1)
}

# Compute expected lables for objects
p = {}
for(i in 1:nrow(x)){
  temp = apply(cbind(1:k), 1, function(j) {pi[j] * dmvnorm(x[i,], mean = mu[j,], sigma = sigma[j,,])})
  temp = temp/sum(temp)
  p = cbind(p,temp)
}

#update the means
for(i in 1:k){
  mu[i,] = c(0,0)
  for(j in 1:nrow(x)){
    mu[i,] = mu[i,] + p[i,j] * x[j,]
  }
  mu[i,] = mu[i,]/sum(p[i,])
}

# Update the covariance
for(i in 1:k){
  sigma[i,,] = c(0,0,0,0)
  for(j in 1:nrow(x)){
    sigma[i,,] = sigma[i,,] + p[i,j] * ((x[j,] - mu[i,]) %*% t(x[j,] - mu[i,]))
  }
  sigma[i,,] = sigma[i,,]/sum(p[i,])
}


# Update mixing coefficient
for(i in 1:k){
  pi[i] = 0
  for(j in 1:nrow(x)){
    pi[i] = pi[i] + p[i,j]/nrow(x)
  }
}


plotdata = as.data.frame(x)
plotdata$num = apply(cbind(1:nrow(x)), 1, function(i){which.max(p[,i])})
pl <- ggplot(plotdata, aes(x= V1, y=V2))
pl + geom_point(aes(colour = plotdata$num))
