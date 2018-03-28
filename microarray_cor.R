##Charles C.N. Wang 0228
#example:GSE48213_RAW

library(dplyr)
nameList <- list.files("/Users/charleswang/Desktop/GSE48213_RAW/")
matrix <- read.table(paste0("/Users/charleswang/Desktop/GSE48213_RAW/",nameList[1]),header = T)
for (i in 2:length(nameList)){
  matrix <- inner_join(matrix,
                       read.table(paste0("/Users/charleswang/Desktop/GSE48213_RAW/",nameList[i]),header = T),
                       by="EnsEMBL_Gene_ID")
}

#save(matrix,file = "56_cell_expression.Rda")
write.table(matrix,"/Users/charleswang/Desktop/GSE48213_RAW/56_cell_expression.txt", sep="\t")

############

#pathfile <- system.file("extdata","expression_example.txt",package = "coexnet")
pathfile <- "/Users/charleswang/Desktop/GSE48213_RAW/56_cell_expression.txt"
cell_data <- read.table(pathfile, sep = "\t")
cell_data$EnsEMBL_Gene_ID

row.names(cell_data) <- cell_data$EnsEMBL_Gene_ID
cell_data$EnsEMBL_Gene_ID = NULL

#Impute Missing Value
library(missForest)
complete.cases(cell_data)
rm_data <- cell_data[complete.cases(cell_data), ]

#cor(x, y, method = c("pearson", "kendall", "spearman"))
#mc_data <- rm_data[,2:length(rm_data)]
#mc_data <- rm_data[,2:10]
#t_mc_data <- t(rm_data)

# Finding threshold value
#library(coexnet)
#cor_pearson <- findThreshold(expData = rm_data,method = "correlation")
#cor_pearson
#threshold=0.8

# Finding PCC value
simil<- round(cor(rm_data, method="pearson"),2)

#plot
library(GGally)
library(ggplot2)
ggcorr(simil, label = TRUE)

#Create an empty matrix

Ad <- matrix(0,ncol = nrow(simil),nrow = nrow(simil))

# Transform to adjacency matrix
threshold=0.8
for(i in seq_len(nrow(simil))){
  Ad[which(simil[,i]>=threshold),i]<-1
  Ad[which(simil[,i]<threshold),i]<-0
}

#for (i in seq_len(nrow(simil))){
#  for (j in seq_len(nrow(simil))){
#    if (simil[i,j]>=threshold){
#      Ad[i,j] <- simil[i,j] 
#    }
#  }
#}

# Change the names of adjacency matrix to the names of the genes
colnames(Ad)<-rownames(Ad)<-rownames(simil)

# Diagonal equal to zero
diag(Ad)<-0

# Creating the network from the adjacency matrix
#library(igraph)
#Gr=graph.adjacency(Ad,mode="undirected",add.colnames=NULL,diag=FALSE)

# Obtaining the degree value for each node in the network

#de <- degree(Gr,loops = FALSE)

# Delete the nodes without any edge

#Ad <- Ad[which(de > 0), which(de > 0)]

# Create the final network from the adjacency matrix filtered

#net <- graph.adjacency(Ad,mode="undirected",add.colnames=NULL,diag=FALSE)
#l <- layout_on_sphere(net)
#plot(net,layout=l)

##############
index <- which(Ad==1, arr.ind=TRUE)
index=as.data.frame(index)

newAd= as.data.frame(paste(rownames(Ad)[index[,1]], colnames(Ad)[index[,2]], sep="-"))

te1=function(x,j){
  b<-simil[index$row[x],index$col[j]]
  return(b)
}

newval=mapply(te1,1:length(index$row),1:length(index$col))
newAd$value = newval
colnames(newAd)=c("name","value")

write.csv(newAd,"/Users/charleswang/Desktop/newAD.csv",row.names = F)

