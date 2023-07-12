## BUILDING VIRUS-RECEPTOR NETWORK

# Load matrix
my_rec <- read.delim ("receptor_virus_matrix.tab", sep = "\t", header=TRUE, na.strings="")
dim(my_rec)
## [1] 246 230

# Boolean matrix           
x <- as.matrix(my_rec[,-1])
x[is.na(x)] <- 0
x <- apply(x, 2,  function(x) as.numeric(x > 0))  #recode as 0/1
v <- x %*% t(x)                                   #the magic matrix 
diag(v) <- 0                                      #repalce diagonal
dimnames(v) <- list(my_rec[, 1], my_rec[,1])                #name the dimensions
v

# Export adjacency matrix to cytoscape format
library(WGCNA)

# convert connection values to [-1,1] range
v_df <- as.data.frame(v)
max(v_df)
## [1] 23
my_fun <- function(x) {x/23}
v_df[] <- lapply(v_df, my_fun)

my_cytosc <- exportNetworkToCytoscape(
  v_df,
  edgeFile = NULL,
  nodeFile = NULL,
  weighted = TRUE,
  threshold = 0.0, # to export ALL connections
  nodeNames = NULL,
  altNodeNames = NULL,
  nodeAttr = NULL,
  includeColNames = TRUE)

edge_df <- as.data.frame(my_cytosc[["edgeData"]])
write.table(file="Cytoscape_RECEPTOR_edgeData_rec.txt", edge_df, sep="\t", na="")

