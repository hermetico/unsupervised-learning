## 1
library(readr)
## A)
dataset <- read_delim("4 -  GAP - GMM/gap_dataset.dat", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
dataset <- dataset[,-1]

## B)
pca <- prcomp(dataset, rank=2)
R <- predict(pca, newdata=dataset)

## C)
null_reference <- function(base, pca){
  min_x = min(base[,1])
  max_x = max(base[,1])
  min_y = min(base[,2])
  max_y = max(base[,2])
  n = nrow(base)
  
  x = runif(n, min_x, max_x)
  y = runif(n, min_y, max_y)
  
  return(data.frame('x' = x, 'y' = y))
}

reference <- null_reference(R, pca)
plot(R)
