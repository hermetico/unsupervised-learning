library(readr)
bone_marrow_df <- read_delim("3 - hierarchical clustering/bone_marrow.dat", 
                          "\t", escape_double = FALSE, col_names = FALSE, 
                          trim_ws = TRUE)

# 1 
dummy = matrix(rnorm(15), nrow=3)
## a) Manhattan distance

manhattan <- function(a, b){
  return(sum(abs(a - b)))
}

dist_manhattan <- function(data){
  size <- nrow(data)
  R <- matrix(0, ncol=size, nrow=size)
  for(i in 1:nrow(data)){
    for(j in i:nrow(data)){
      if(!i == j){
        r = manhattan(data[i,], data[j,])
        R[i,j] = r
        R[j,i] = r
      }
    }
  }
  return (R)
}

dist(dummy, method="manhattan", upper=TRUE, diag=TRUE) 
dist_manhattan(dummy)
## b) Pearson Correlation
pearson <- function(a, b){
  a_ = mean(a)
  b_ = mean(b)
  return(
    sum((a - a_) * (b - b_)) / 
      sqrt( sum((a - a_)^2) * sum((b - b_)^2) )
    )
  
}

dist_pearson <- function(data){
  size <- nrow(data)
  R <- matrix(0, ncol=size, nrow=size)
  for(i in 1:nrow(data)){
    for(j in i:nrow(data)){
      
      if(!i == j){
        r = cov(dummy[i,],dummy[j,]) / sd(dummy[i,]) * sd(dummy[j,])
        #r = pearson(data[i,], data[j,])
        R[i,j] = r
        R[j,i] = r
      }
    }
  }
  return (R)
}
dist_pearson(dummy)

cov(dummy[1,],dummy[2,]) / sd(dummy[1,]) * sd(dummy[2,])
## c) Spearman's Rank Correlation