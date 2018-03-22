# Exercise Sheet 2 - Unsupervised Learning
# Exercise solution skeleton

# For plotting we will make use of the ggplot2 library
library(ggplot2)
library(readr)

# to inspect function description use ?func




# 1. PCA
# (a) Calculate a centered version of the dataset X
# first column as rownames
dataframe <- read_delim("2 - pca/uk.dat", "\t", escape_double = FALSE, trim_ws = TRUE)
dataframe = as.matrix(data.frame(dataframe[,-1], row.names=dataframe$X1))
# read in ./uk.dat as numercial matrix and assign row names ~ 3 line of code

# to center a matrix, use scale() - it requries a - numeric - matrix as imput (see ?scale ), but centers columns wise - we have rows of values!
# ~1 line of code
centered <- t(scale(t(dataframe), scale = FALSE))
M <- centered
# (b) Calculate the Covariance matrix C
# The covariance matrix Q of a centered matrix M is defined as Q = 1/(n-1)*(M,M^T) 
# where n is the number of columns in M and M^T is the transposed matrix M
# ~ 1 line of code
Q <- centered %*% t(centered) * (1/ncol(centered) - 1)
# (c) Calculate the Eigenvalues and Eigenvectors of C
# use builtin function eigen() to compute the eigenvalues and vectors. You can separate them by using the accessor $values and $vectors
# ~ 3 lines of code
eig <- eigen(Q)
eigenvalues = eig$values
eigenvectors = eig$vectors

# (d) Calculate and plot the coordinates of our 4 objects using 
# i  - The first principal component 
# ii - The first and second principal component


# For this we need to extract the first and second eigenvector, use cbind() to form them back into a matrix
# ~ 2 lines of code
f_first <- cbind(eigenvectors[,1])
f_second <- cbind(eigenvectors[,1], eigenvectors[,2])
# Multiply the eigenvector with the centerd matrix  to get the values of the pricipal components
# ~ 2 lines of code
E <- t(f_first) %*% M 
E2 <- t(f_second) %*% M 
# plot by country
# plot the coordinates of the first principal component - use qplot - quickplot from ggplot2
# ~ 5 lines of code
plot.data <- data.frame('x' = t(E), 'y' = 0)
ggplot(plot.data, aes(x=x, y=y))  + geom_point(aes(shape=rownames(principal)))


# plot the coordinates of the first two principal components
# ~ 4 lines of code
plot.data <- data.frame('x' = t(E2)[,1], 'y' = t(E2)[,2])
ggplot(plot.data, aes(x=x, y=y))  + geom_point(aes(shape=rownames(principal)))

# plot by cosumable
# First PC

# First 2 PC
dataframe <- read_delim("2 - pca/denmark.dat", "\t", escape_double = FALSE, trim_ws = TRUE)
dataframe = as.matrix(data.frame(dataframe[,-1], row.names=dataframe$X1))
# ===================================================================

# 2. Multidimensional Scaling
# (a) Calculated the Matrix B as explained in the lecture
# B = X^T x X
X = dataframe
# We provide two datasets, denmark.dat and germany.dat

# drop the X column for the matrix again and assign row names

# convert the distances matrix into B as described in the slides (p. 57)
# A = (a_ij) = -1/2* d_ij ^2
# ~ 1 line of code
A <- -0.5 * X * X
A_row_means <- rowMeans(A)
A_col_means <- colMeans(A)
meanA <- mean(A)
B <- matrix(0, nrow = nrow(A), ncol = ncol(A))
for( i in 1:nrow(A)){
  for(j in 1:ncol(A)){
    B[i,j] <- A[i,j]  - A_row_means[i] - A_col_means[j] + meanA
  }
    
}
colnames(B) <- colnames(A)
row.names(B) <- rownames(A)
# Calculate the average of the row in A, the average of the columns in A and the average of the matrix
# You can use the function mean in combination with apply () - see ?apply for rowwise and column wise application 
# ~ 3 lines of code

# define B = (b_ij) as b_ij = a_ij - a_i[row] - a_j[column] + mean(A)
# where a_i[row] is the average of row i, a_j[column] is the average of column j and mean(A) the average of A
# First copy dimensions of A to B, then overwrite all values
# your for loop here

# (b) Calculate the Eigenvalues and Eigenvectors of B

# as before 
# ~ 2 lines of code
eig <- eigen(B)
eigenvalues = eig$values
eigenvectors = eig$vectors

# (c) Use the first two Eigenvectors of B to calculate the tranformed coordinates of the cities
# As before use cbind to assign and multiply with B
# ~ 2 lines of code
f_first <- cbind(eigenvectors[,1])
f_second <- cbind(eigenvectors[,1], eigenvectors[,2])

E <- t(f_first) %*% B 
E2 <- t(f_second) %*% B
# plot the original data

plot.data <- data.frame('x' = t(E2)[,1], 'y' = t(E2)[,2])
ggplot(plot.data, aes(x=x, y=y))  + geom_point(aes(shape=rownames(plot.data)))
# fimd a rotation facor alpha_denmark and alpha_germany that rotates the points in the regular projection (e.g. 0.3*pi)
#Rotation Denmark

#Rotation Germany

# To apply rotation, we multiply with a 2x2 matrix:
# 2x2 matrix:
# cos(alpha), sin(alpha)
# -sin(alpha), cos(alpha)
alpha = 0.9 * pi
R = matrix( c(cos(alpha), sin(alpha),-sin(alpha), cos(alpha)),nrow=2, ncol=2)
rotated_data <-  as.matrix(plot.data) %*%R
plot.data <- data.frame('x' = rotated_data[,1], 'y' = rotated_data[,2])
ggplot(plot.data, aes(x=x, y=y))  + geom_point(aes(shape=rownames(plot.data)))
# plot the result again with qplot()
# alpha = alpha_denmark
# alpha = alpha_germany
# rotation.M = matrix(c(cos(alpha_denmark), sin(alpha_denmark), -sin(alpha_denmark), cos(alpha_denmark)), ncol=2)