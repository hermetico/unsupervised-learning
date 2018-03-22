library(ggplot2)

# Exercise 1 : Implement the following distance measures by yourself and compare to the built-in capabilities of R:  
# (a) Manhattan distance
# (b) Pearson Correlation
# (c) Spearman's Rank Correlation 

# define the measures as measure.manhattan, measure.pearson and measure.spearman
# Remember, the  manhatten distance is given by:
# d(x,y) = sum _i=1 ^n  ( |x_i - y_i| )
measure.manhattan <- function(x,y){
  return(sum(abs(x - y)))
}

# The pearson correlation is calculated by:
# phi(x,y) = sum_i=1 ^n ((x_i-x_i_mean)*(y_i-y_i_mean)) / sqrt( sum_i=1 ^n (x_i - x_i_mean)^2  *  sum_i=1 ^n (y_i-y_i_mean)^2 )
measure.pearson <- function(x,y){
  mean.x <- mean(x)
  mean.y <- mean(y)
  return(
    sum((x - mean.x) * (y - mean.y)) / 
      sqrt( sum((x - mean.x)^2) * sum((y - mean.y)^2) )
  )
}

# The spearman correlation is defined as the Pearson correlation of the ranked variables.
# make use of the builtin function rank()
measure.spearman <- function(x,y){
  measure.pearson(rank(x), rank(y))
}

# 2. Download the dataset for the website. The file is tab-separated. Each row represents one patient. 
# The row starts wuth the row name followed by 999 expression values.
# All patients with a name like "ALL_19769_B.cell" are of typ e B-ALL, 
# all patient with a name like "ALL_19881_T.cell_2" are of type T 
# and all patients with a name like "AML_16" are of type ALL.


# Read in the datatable and assign row names 
inp <- read.table("3 - hierarchical clustering/bone_marrow.dat", sep='\t')
X <- data.matrix(inp)
row.names(X) = inp$V1

# 3. Calculate a 38 Ã— 38 proximity matrix P of all pair-wise proximities by using
# (a) Manhattan distance -> dist(x, method="manhattan")
# (b) Pearson Correlation -> cor(x, method="pearson")
# (c) Spearman's Rank Correlation -> cor(x, method="spearman")


# Implement a function, that applies a measure m to all pairs of row-vectors of s and returns a distance matrix t
generate.matrix <- function(s, m, init){
  t <- matrix(init, nrow(s), nrow(s))
  # compute distance matrix here
  for(i in 1:(nrow(s)-1)){
    for(j in (i+1):nrow(s)){
        r = m(s[i,], s[j,])
        t[i,j] = r
        t[j,i] = r
    }
  }
  return(t)
}


# Initialize the matrix for measures, use dist() for computing the already implemented measure. 
# Make sure to use the correct initialization parameter - 0 for similarity measures and 1 for distances
dist.manhattan <- generate.matrix(X, measure.manhattan, 0)
dist.manhattan.r <- dist(X, method="manhattan")

sim.pearson <- generate.matrix(X, measure.pearson, 1)
sim.pearson.r <- cor(t(X), method="pearson")

sim.spearman <- generate.matrix(X, measure.spearman, 1)
sim.spearman.r <- cor(t(X), method="spearman")


# 4. Calculate a hierarchical clustering (set the cutoff such that you receive 3 clusters) based on the calculated matrices above.
# Be aware that you might have to convert between similarities and distances. Use hierarchical clustering with the following linkage functions:
# (a) Single linkage clustering
# (b) Complete linkage clustering
# (c) Average linkage clustering

# We can use hclust() - which conveniently provides builtin methods for single, average and complete Linkage
# But: hclust wants to have distance measures stored as dist-objects therefore we need to convert similarity measures to distance measures
dist.manhattan <- as.dist(dist.manhattan)
dist.pearson.r <- as.dist(1 - sim.pearson.r) 
dist.spearman.r <- as.dist(1 - sim.spearman.r)


# Let's cluster all 9 combinations using hclust:
hc.man.avg <- hclust(dist.manhattan.r, method="average")
hc.man.max <- hclust(dist.manhattan, method="complete")
hc.man.min <- hclust(dist.manhattan, method="single")

hc.pear.avg <- hclust(dist.pearson.r, method="average")
hc.pear.max <- hclust(dist.pearson, method="complete")
hc.pear.min <- hclust(dist.pearson, method="single")

hc.spear.avg <- hclust(dist.spearman.r, method="average")
hc.spear.max <- hclust(dist.spearman, method="complete")
hc.spear.min <- hclust(dist.spearman, method="single")

# Lets have a look at the dendrograms
plot(hc.pear.avg, hang = -1)
plot(hc.pear.max, hang = -1)
plot(hc.pear.min, hang = -1)

plot(hc.spear.avg, hang = -1)
plot(hc.spear.max, hang = -1)
plot(hc.spear.min, hang = -1)

plot(hc.man.avg, hang = -1)
plot(hc.man.max, hang = -1)
plot(hc.man.min, hang = -1)


# 5. Compare the results. What is the best method when comparing to the gold-standard?
# Find an appropriate way of displaying the results.

# convert to dendrogram for comparison
den.pear.avg <- as.dendrogram(hc.pear.avg)
den.pear.max <- as.dendrogram(hc.pear.max)
den.pear.min <- as.dendrogram(hc.pear.min)

den.spear.avg <- as.dendrogram(hc.spear.avg)
den.spear.max <- as.dendrogram(hc.spear.max)
den.spear.min <- as.dendrogram(hc.spear.min)

den.man.avg <- as.dendrogram(hc.man.avg)
den.man.max <- as.dendrogram(hc.man.max)
den.man.min <- as.dendrogram(hc.man.min)

# define vector of colorbliend friendly colors for classes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# barplot(seq(1,length(cbPalette)), col=cbPalette)

colorLabelsByLeafClassEdgesByMembership <- function(n, membership){
  if(is.leaf(n)){
    # Find the attributes of current node
    # TODO pattern match by class labels with grepl and assign color
    # TODO determine edge color by label
    a <- attributes(n)
    if(grep("AML", a$label)){
      attr(n, "nodePar") <- c(a$nodePar, list(col=cbPalette[1], lwd=4))
    }
    

    edgeCol <- cbPalette[membership[which(names(membership) == a$label)]]
    attr(n, "edgePar") <- c(a$edgePar, list(col=edgeCol, lwd=4))
  }
  
}

# Visualize clustering with manhattan distance


# Visualize clustering with Pearson correlation


# Visualize clustering with Spearman correlation



# How to compare it to the gold standard - what is the gold standard?