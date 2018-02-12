## Exercice 2
library(readr)
dataframe <- read_csv("dummy_data.csv", col_names = FALSE)
matrix_a = data.matrix(dataframe, rownames.force = NA)
matrix_b = t(data.matrix(dataframe, rownames.force = NA))

## Exercice 3
dummy_matmult <- function(A, B){
  C  = matrix(0, nrow = nrow(A), ncol = ncol(B))
  for(i in 1:nrow(A))
  {
    for(j in 1:ncol(B))
    {
      acc = 0
      for(k in 1:ncol(A)){
        acc = acc + A[i,k] * B[k,j]
      }
      C[i,j] = acc
    }
  }
  return(C)
}

dummy_matmult(matrix_a, matrix_b)
matrix_a %*% matrix_b

## Exercice 4
### square matrix
matrix_a = cbind(matrix_a, matrix_a[,5])
### kernel computation
dummy_kernel <- function(A){
  return(A[1,1] * A[2,2] - A[2,1] * A[1,2])
}
### the algorithm
dummy_determinant <-function(A){
  
  if(nrow(A) == 2){
    return (dummy_kernel(A))
  }
  
  acc = 0 
  sign = 1
  for(i in 1:nrow(A)){
    A_ = A[-i,-1]
    num = A[i,1]
    
    #print(sign * num)
    #print(A_)
    
    acc = acc + sign * (num * dummy_determinant(A_))
    sign = sign * -1
  }
  return(acc)
}

dummy_determinant(matrix_a)
as.integer(det(matrix_a))

## Exercice 5
dummy_random_matrix <-function(ncol=1, nrow=1){
  # numbers between 0 and 1
  M = matrix(0, ncol=ncol, nrow=nrow)
  for(i in 1:nrow){
    for(j in 1:ncol){
      M[i,j] = runif(1,0,1)
    }
  }
  return (M)
}

# benchmark
tests = 1000

size = 1:tests
dummy = c()
builtin = c()

for(i in 1:tests){
  dummy = c(dummy, system.time(dummy_random_matrix(i,i))[1])
  builtin = c(builtin, system.time(matrix(rnorm(i*i),i,i))[1])
}

bench <- data.frame(size=size, dummy=dummy, builtin=builtin)
plot(bench$dummy, type="o", col="blue")
lines(bench$builtin,type="o", col="red")
legend(1,1, c("dummy","builtin"), cex=0.8, col=c("blue","red"), pch=21:22, lty=1:2)
