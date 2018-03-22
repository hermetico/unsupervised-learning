# Class notes

c = c(1,2,3,4,5,6)

mean(c) # the mean
rnorm(c) # normal distribution

seq(0, 10, 0.33) # creates a sequence with init, end step

for(i in seq(0, 10, 0.33)){
  print(i)
}

c = c(1,2,3)
c[6] = 5 # it fills the gaps
mean(c) # fails because na
mean(c, na.rm=TRUE) # ignores na
table <- read.table("dummy_data.csv", sep=',') # reads a table without dataframe
runif(5) # uniform distribution of 5 elements
