# Installing package
install.packages(c("caret", "ggplot2", "tidyverse", "tsne", "umap", "Rtsne" , "factoextra", "mclust", "matrixStats", "plotly"))

# Loading package
library(caret)
library(ggplot2)
library(tidyverse)
library(tsne)
library(umap)
library(Rtsne)
library(factoextra)
library(mclust)
library(matrixStats)
library(plotly)
library("RColorBrewer")

#Path of the data set files
datatsets_path = "~/MASTER/2nd Semester/ML in Computational Biology/Assignments/2/Assignment2_Datasets"
#Change the working directory
setwd(datatsets_path)

#Load data set
df <- read.csv("dataset1.csv")
#Turn first column into row index
df <- data.frame(df[,-1], row.names = df[,1])

par(mfrow=c(1,1))
boxplot(df, main="Barplot of original Dataset 1", 
        xlab="Genes")

#####################Remove highly correlated variables (greater than 0.7)
data <- df
Corrmatrix <- cor(data, method = "pearson")
highlyCorrelated <- findCorrelation(Corrmatrix, cutoff = .7, verbose = FALSE, names = TRUE)
highlyCorrelated
# Remove highly correlated variables and create a new dataset
df <- data[, ! names(data) %in% highlyCorrelated]
data <- df

#####################Data preprocessing
#1.MinMax normalization
mean_norm_minmax <- function(x){
  (x- mean(x)) /(max(x)-min(x))
}
minmax_df <- as.data.frame(lapply(data, mean_norm_minmax))
boxplot(minmax_df, main="Barplot of Min-Max normalization Dataset 1", 
        xlab="Genes")

#2. Log2 normalization
logNorm_df <- log2(data+1) #+1 to make up for the 0.0 values
boxplot(logNorm_df, main="Barplot of logNorm Dataset 1", 
        xlab="Genes")

#3. Centered and Scaled
standarized <- preProcess(data, method=c("center", "scale"))
standarized_df <- predict(standarized, data)
boxplot(standarized_df, main="Barplot of Centered and Scaled Dataset 1", 
        xlab="Genes")

#4. Centered and Scaled and LogNorm
standarized_log <- preProcess(logNorm_df, method=c("center", "scale"))
standarized_log_df <- predict(standarized_log, logNorm_df)
boxplot(standarized_log_df, 
        main="Barplot of Centered and Scaled and LogNorm Dataset 1", 
        xlab="Genes")

#5. Centered and LogNorm
centered_log <- preProcess(logNorm_df, method=c("center"))
centered_log_df <- predict(centered_log, logNorm_df)
boxplot(centered_log_df, 
        main="Barplot of Centered and LogNorm Dataset 1", 
        xlab="Genes")


data <- standarized_df

##################### PCA
plot.new()

#Perform PCA
df_pca <- prcomp(data)

fviz_pca_ind(df_pca, col.ind = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

#Find the explained variance per df_pca
expl_var <- df_pca$sdev^2/sum(df_pca$sdev^2)

#How many principal components do we need to reach the explained variance threshold (=0.8)
count <- 0
threshold <- 0.8
sum <-0
for(k in expl_var){
   count <- count + 1
   sum <- sum + k
   if (sum >= threshold){
     break
   }
}
print(count)


N_perm <- 10
expl_var_perm <- matrix(NA, ncol = length(df_pca$sdev), nrow = N_perm)
for(k in 1:N_perm)
{
  expr_perm <- apply(df,2,sample)
  PC_perm <- prcomp(t(expr_perm), center=TRUE, scale=FALSE)
  expl_var_perm[k,] <- PC_perm$sdev^2/sum(PC_perm$sdev^2)
}

#Optimal number of principal components
plot(expl_var[1:50]~seq(1:50), ylab="EXPLAINED VARIANCE",
     col="darkgreen", type='o', xlab="PRINCIPAL COMPONENTS")
lines(colMeans(expl_var_perm)[1:50]~seq(1:50),col="red")
legend("topright", c("Explained by PCs", "Explained by chance"),
       fill=c("darkgreen","red"), inset=0.02)

pval <- apply(t(expl_var_perm) >= expl_var,1,sum) / N_perm
optPC<-head(which(pval>=0.05),1)-1
mtext(paste0("OPTIMAL NUMBER OF PRINCIPAL COMPONENTS = ",optPC))

#3D scatter plot (3 PCs)
components <- data.frame(df_pca[["x"]])
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3

fig <- plot_ly(components, 
               x = ~PC1, y = ~PC2, z = ~PC3)
fig <- fig %>%
  layout(
    title = "3D scatter plot",
    scene = list(bgcolor = "#e5ecf6")
  )
fig



##################### t-SNE
data <- standarized_df


#Test for the optimum number of components (=2 for dataset1)
N_cells<-dim(data)[2]
par(mfrow=c(3,3))
perp_range<-vector(length=9)
perp_range[1]<-3; perp_range[9] <- N_cells/3-1
optPerp <- round(sqrt(N_cells),0); 
perp_step<-(optPerp-3)/4
for(s in 2:8){perp_range[s]<-perp_range[s-1]+perp_step}
for(j in round(perp_range,0))
{
  tsne_perp_iter<-Rtsne(data,perplexity=j,
                        initial_dims=optPC,max_iter=5000)
  plot(tsne_perp_iter$Y,col="blue",xlab="tSNE1",ylab="tSNE2", bg='tomato2', pch=21, cex=1, lwd=0.2)
  mtext(paste0("perplexity = ",j))
  
}
mtext(paste0("tSNE PC=", optPC), side = 3, line = -2, outer = TRUE,font=2)



#Test for 30 initial components
N_cells<-dim(data)[2]
par(mfrow=c(3,3))
perp_range<-vector(length=9)
perp_range[1]<-3; perp_range[9] <- N_cells/3-1
optPerp <- round(sqrt(N_cells),0); 
perp_step<-(optPerp-3)/4
for(s in 2:8){perp_range[s]<-perp_range[s-1]+perp_step}
for(j in round(perp_range,0))
{
  tsne_perp_iter<-Rtsne(data,perplexity=j,
                        initial_dims=30,max_iter=5000)
  plot(tsne_perp_iter$Y,col="blue",xlab="tSNE1",ylab="tSNE2", bg='tomato2', pch=21, cex=1, lwd=0.2)
  mtext(paste0("perplexity = ",j))

}
mtext(paste0("tSNE PC=", 30), side = 3, line = -2, outer = TRUE,font=2)



#Test for number of components equal to 80% of the dataset's variance (=59 for dataset1)
N_cells<-dim(data)[2]
par(mfrow=c(3,3))
perp_range<-vector(length=9)
perp_range[1]<-3; perp_range[9] <- N_cells/3-1
optPerp <- round(sqrt(N_cells),0); perp_step<-(optPerp-3)/4
for(s in 2:8){perp_range[s]<-perp_range[s-1]+perp_step}
for(j in round(perp_range,0))
{
  tsne_perp_iter<-Rtsne(data,perplexity=j,
                        initial_dims=15,max_iter=2000)
  plot(tsne_perp_iter$Y,col="blue",xlab="tSNE1",ylab="tSNE2", bg='tomato2', pch=21, cex=1, lwd=0.2)
  mtext(paste0("perplexity = ",j))
  
}
mtext(paste0("tSNE PC=", 15), side = 3, line = -2, outer = TRUE,font=2)


#Test number of iterations
par(mfrow=c(3,3))
optPerp <- round(sqrt(N_cells),0)
for(i in c(50, 100, 200, 500, 1000, 2000, 3000, 5000, 10000))
{
  tsne <- Rtsne(data, initial_dims=30,
                perplexity=15, max_iter=i)
  plot(tsne$Y, col="blue", xlab="tSNE1", ylab="tSNE2", bg='tomato2', pch=21, cex=1, lwd=0.2)
  mtext(paste0("max_iter = ", i))
}
mtext(paste0("tSNE PC=30 Perplexity=15"), side = 3, line = -2, outer = TRUE,font=2)


###Final tSNE results
tsne_results <- Rtsne(data, perplexity = 15, 
                      initial_dims=30, 
                      max_iter=5000)                     
par(mfrow=c(1,1))
plot(tsne_results$Y,col="blue",xlab="tSNE1",ylab="tSNE2", bg='tomato2', pch=21, cex=1, lwd=0.2)
mtext(paste0("tSNE\n Initial components=30, Perplexity=15, Max. Iter.=5000"),
      font=2)


##################### UMAP
data <- standarized_df
df_pca <- prcomp(data)
par(mfrow=c(4,5))
for(j in c(7, 15, 25, 35))
{
  for(i in c(0.05, 0.1, 0.2, 0.3, 0.4))
  {
    df.umap <- umap(df_pca$x[,1:50], n_epochs=2000, n_neighbors=j, min_dist=i)
    plot(df.umap$layout, xlab="UMAP1", ylab="UMAP2", col="blue", bg='tomato2', pch=21, cex=1, lwd=0.2)
    mtext(paste0("neigh.=", j, " | dist=", i))
  }
}
mtext(paste0("UMAP PC=", 50), side = 3, line = -2, outer = TRUE,font=2)

###Final UMAP results
par(mfrow=c(1,1))
df.umap <- umap(df_pca$x[,1:50], n_epochs=2000, n_neighbors=25, min_dist=0.1)
plot(df.umap$layout, xlab="UMAP1", ylab="UMAP2", col="blue", bg='tomato2', pch=21, cex=1, lwd=0.2)
mtext(paste0("UMAP\nn_neighbors=", 25, " | dist=", 0.1), font=2)




##################### MLCLUST


########### PCA
df.pca <- data.frame(x = df_pca$x[,1],
                     y = df_pca$x[,2])

par(mfrow=c(1,3))

BIC1 <- mclustBIC(df.pca)

cluster.data <- Mclust(df.pca, x = BIC1)
#df.PCA.MM1 <- as.data.frame(df_pca$rotation[,1:3])

## preview plots
par(mfrow=c(1,3), oma=c(0,0,2,0))

plot(df.pca, col="blue", xlab="PCA1", ylab="PCA2", 
     bg='tomato2', pch=21,  lwd=0.2)
mtext("PC1 vs PC2", font=2)

plot(BIC1)
mtext("BIC criterion", font=2)

plot(cluster.data, bg=cluster.data$classification, pch=21,
     xlab='PC1', ylab='PC2', what="classification")
mtext("PCA Clustering", font=2)


title(paste0("PCA"), cex.main = 2, outer=TRUE)

#Plots
post <- as.data.frame(cluster.data[["z"]])
post[post < 0.01] <- 0

names = c()
for (i in 1:ncol(post)){
  names = c(names, paste0("Cluster ",i))
}
colnames(post) <- names
pal <- brewer.pal(ncol(post), "Dark2")

par(mfrow=c(1,1))
post1 <- post[do.call(order, post),]
matplot(y=post1, type =  "n",  xlab = "Cell index", ylab = "Posterior",
        main = paste("PCA Trajectory plot"), lty = "solid")
matlines(post1,   col = pal, lty = "solid")
legend("right", legend=names, col = pal,lty=1:2, cex=1) 





############ t-SNE
df_tsne <- data.frame(x = tsne_results$Y[,1],
                      y = tsne_results$Y[,2])


BIC2 <- mclustBIC(df_tsne)
cluster.data <- Mclust(df_tsne, x = BIC2)
df.PCA.MM1 <- as.data.frame(df_tsne)

## preview plots
par(mfrow=c(1,3), oma=c(0,0,3,0))

plot(tsne_results$Y, col="blue", xlab="tSNE1", ylab="tSNE2", 
           bg='tomato2', pch=21,  lwd=0.2)
mtext("tSNE1 vs tSNE2", font=2)

plot(BIC2)
mtext("BIC criterion", font=2)

plot(cluster.data, bg=cluster.data$classification, pch=21,
     xlab='tSNE1', ylab='tSNE2', what="classification")
mtext("tSNE Clustering", font=2)


title(paste0("tSNE\n Perplexity=", 15, " Initial Dimenstions=", 30), cex.main = 2, outer=TRUE)

#Plots
post <- as.data.frame(cluster.data[["z"]])
post[post < 0.01] <- 0

names = c()
for (i in 1:ncol(post)){
  names = c(names, paste0("Cluster ",i))
}
colnames(post) <- names
pal <- brewer.pal(ncol(post), "Dark2")

par(mfrow=c(1,1))
post1 <- post[do.call(order, post),]
matplot(y=post1, type =  "n",  xlab = "Cell index", ylab = "Posterior",
        main = paste("tSNE Trajectory plot"), lty = "solid")
matlines(post1,   col = pal, lty = "solid")
legend("right", legend=names, col = pal,lty=1:2, cex=0.8) 



########### UMAP
df_umap <- data.frame(x = df.umap$layout[,1],
                      y = df.umap$layout[,2])

par(mfrow=c(1,1))
plot(df.umap$layout, xlab="UMAP1", ylab="UMAP2", col="blue", bg='tomato2', pch=21, cex=1, lwd=0.2)
mtext(paste0("UMAP\nn_neighbors = 25, min_dist = 0.1,	n_epochs = 2000"), font=2)

BIC3 <- mclustBIC(df_umap)
cluster.data <- Mclust(df_umap, x = BIC3)

## preview plots

par(mfrow=c(1,3), oma=c(0,0,3,0))

plot(df.umap$layout,xlab="UMAP1",ylab="UMAP2",
     col="blue", bg='tomato2', pch=21, lwd=0.2)
mtext("UMAP1 vs UMAP2", font=2)

plot(BIC3)
mtext("BIC criterion", font=2)

plot(cluster.data, bg=cluster.data$classification,
     pch=21,xlab='UMAP1', ylab='UMAP2', what = "classification")

mtext("UMAP Clustering", font=2)


title(paste0("UMAP\nn_neighbors = 25, min_dist = 0.1,	n_epochs = 2000"), cex.main = 2, outer=TRUE)



#Plots


post <- as.data.frame(cluster.data[["z"]])
post[post < 0.01] <- 0

names = c()
for (i in 1:ncol(post)){
  names = c(names, paste0("Cluster ",i))
}
colnames(post) <- names
pal <- brewer.pal(ncol(post), "Dark2")

par(mfrow=c(1,1))
post1 <- post[do.call(order, post),]
matplot(y=post1, type =  "n",  xlab = "Cell index", ylab = "Posterior",
        main = paste("UMAP Trajectory plot"), lty = "solid")
matlines(post1,   col = pal, lty = "solid")
legend("right", legend=names, col = pal,lty=1:2, cex=0.8) 







