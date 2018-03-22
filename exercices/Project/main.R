library(readr)
## Step 1 data preprocessing
all_vs_all <- read_delim("Project/all-vs-all.tsv", 
                         "\t", escape_double = FALSE, col_names = c('query', 'target', 3, 4, 5, 6, 7, 8, 9, 10, 'expectation', 12), 
                         trim_ws = TRUE)
# now we create a matrix
sequences = unique(all_vs_all$query)
N = length(unique(all_vs_all$query))
M = matrix(0, nrow=N, ncol=N)
## create a dataframe out of the 0 matrix
M.frame = data.frame(M)
colnames(M.frame) <- sequences
rownames(M.frame) <- sequences

# adds pseudocounts
# applies -log(X) to the column to convert the expectation in similarity
all_vs_all$similarity <- all_vs_all$expectation + 0.0000000000000000000001
all_vs_all$similarity <- -log(all_vs_all$similarity)

# loop over all the all_vs_all keeping the best results
for(i in 1:nrow(all_vs_all)){
  query = all_vs_all$query[i]
  target = all_vs_all$target[i]
  similarity = all_vs_all$similarity[i]
  M.frame[query, target] = max(M.frame[query, target], similarity)
}

# Write CSV in R
write.csv(M.frame, file = "Project/similarity.csv")
