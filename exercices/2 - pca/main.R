library(readr)
library(ggplot2)
# Exercice 1
dataframe <- read_delim("2 - pca/uk.dat", "\t", escape_double = FALSE, trim_ws = TRUE)
# first column as rownames
dataframe = as.matrix(data.frame(dataframe[,-1], row.names=dataframe$X1))

# a)
substract_mean <- function( data ){
  means <- rowMeans(data)
  return(data - means)
}

data <-substract_mean(dataframe)

# b)
#Q <- var(data) # it is the same
#Q <- cov(data) # is equal to (data %*% t(data)) / 3
Q <- (data %*% t(data)) / (nrow(data) - 1)

# c)
E <- eigen(Q)
e_vectors <- E$vectors
e_values <- E$values

# d
######## Do we have to add the mean??
rotated <- ( t(data) %*% e_vectors ) #+ rowMeans(dataframe)
colnames(rotated) <- rownames(data)
# .1)
argument <- order(e_values, decreasing=TRUE)[1]
#components = E$vectors[argument,]
principal <- data.frame(rotated[,argument],'y'=0)
colnames(principal)[1] <- colnames(rotated)[argument]

ggplot(principal, aes(x=Confectionery, y=y))  + geom_point(aes(shape=rownames(principal)))


#legend(1,1, rownames(principal), cex=0.8, col=c("blue","red"), pch=21:22, lty=1:2)
# .2)
arguments <- order(e_values, decreasing=TRUE)[1:2]
#components = E$vectors[argument]
principals <- data.frame(rotated[,arguments])
ggplot(principals, aes(x=Confectionery, y=Alcoholic.drinks))  + geom_point(aes(shape=rownames(principal)))


# Exercice 2
dataframe <- read_delim("2 - pca/denmark.dat", "\t", escape_double = FALSE, trim_ws = TRUE)
dataframe = as.matrix(data.frame(dataframe[,-1], row.names=dataframe$X1))
centered = substract_mean(dataframe)
# a)
B <- t(centered) %*% centered

# b)
BE <- eigen(B)
