setwd("C:/Users/dan91/Rstudio/statistic-in-ecology/stat_data/mock_dataTS")
cop_compo <- read.table("Copepod_composition.txt", header = T)
cop_den <- read.table("Cop_density.txt")
cop_sp <- read.table("copepodSPlist.txt", fill = T, sep = "\n")

sp_den <- matrix(0, nrow = length(cop_sp$V1), ncol = length(cop_den$V1))
for(i in 1:length(cop_sp$V1)){
  for(j in 1:length(cop_den$V1)){
    sp_den[i, j] <- cop_compo[i, j]/100 * cop_den[j,]
  }
}
sp_den.df <- as.data.frame(sp_den)
colnames(sp_den.df) <- colnames(cop_compo)
rownames(sp_den.df) <- cop_sp$V1
head(sp_den.df)

dominant_sp.df <- data.frame(ifelse(cop_compo >= 5, 1, 0))
dominant_species <- rownames(dominant_sp.df)[rowSums(dominant_sp.df) >= 1]
as.numeric(dominant_species)

dom_sp_compo <- t(cop_compo[dominant_species, ])
colnames(dom_sp_compo) <- cop_sp[dominant_species, ]
rownames(dom_sp_compo) <- colnames(cop_compo)
dom_sp_compo.df <- as.data.frame(dom_sp_compo)

#Apply PCA or MDS using matrix algebra to extract major gradients of the dominant species. 
#Make a bi-plot to show the relationships among species and sites. 
#Then, check your results with those from build in functions of PCA, MDS, and NMDS.

#PCA (manually)
dom_sp_center <- scale(dom_sp_compo.df, center = T, scale = F)

cov_matrix <- cov(dom_sp_center)

eigen_decompo <- eigen(cov_matrix)
eigen_vector <- eigen_decompo$vectors
eigen_value <- eigen_decompo$values

pca_score <- dom_sp_center %*% eigen_vector

plot(pca_score[, 1], pca_score[, 2], 
     xlab = "PC1", 
     ylab = "PC2", 
     main = "Biplot")

text(pca_score[,1], pca_score[,2], 
     labels = rownames(dom_sp_compo.df), 
     pos = 4, cex = 0.8)

for (i in 1:ncol(dom_sp_compo.df)) {
  arrows(0, 0, 
         eigen_vector[i, 1] * sqrt(eigen_value[1]) * 2,
         eigen_vector[i, 2] * sqrt(eigen_value[2]) * 2,
         col = "blue", length = 0.1)
  
  text(eigen_vector[i, 1] * sqrt(eigen_value[1]) * 1.8 * 1.2, 
       eigen_vector[i, 2] * sqrt(eigen_value[2]) * 1.8 * 1.2, 
       labels = colnames(dom_sp_compo.df)[i], 
       col = "blue", cex = 0.8)
}

#PCA (automatically)
pca_dom_copepod <- prcomp(dom_sp_center)
biplot(pca_dom_copepod, scale = 0)



dom_sp_dist <- dist(dom_sp_compo.df, method = "euclidean")

D2 <- as.matrix(dom_sp_dist)^2

n <- nrow(D2)
H <- diag(n) - matrix(1, n, n)/n  # centering matrix
B <- -0.5 * H %*% D2 %*% H

eigen_B <- eigen(B)
eigen_values_mds <- eigen_B$values
eigen_vectors_mds <- eigen_B$vectors

# 只取特徵值>0的部分
positive_indices <- which(eigen_values_mds > 0)
Lambda_sqrt <- diag(sqrt(eigen_values_mds[positive_indices]))
V <- eigen_vectors_mds[, positive_indices]

# MDS座標
mds_scores <- V %*% Lambda_sqrt

plot(mds_scores[,1], mds_scores[,2], 
     xlab = "MDS1", ylab = "MDS2", main = "Classical MDS (manual)", pch = 19)


mds_builtin <- cmdscale(dom_sp_dist, k = 2)  # k是要幾維，通常k=2
plot(mds_builtin[,1], mds_builtin[,2], 
     xlab = "MDS1", ylab = "MDS2", main = "Classical MDS (cmdscale)", pch = 19)



