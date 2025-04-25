# Dominant copepod species

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

#Use the dominant copepod species data (from HW 1). 
#Perform cluster analysis of stations based on percent composition data of the dominant species and tell your story about these copepod data. 
#Compare the results based on different distance measures. 
#Compare the results based on different cluster algorithms. 
#Determine final number of clusters and describe the differences among them.

library(vegan)
library(cluster)

#different distance method (hierarchical)
dist_eu <- dist(dom_sp_compo, "euclidean")

hc_eu_average <- hclust(dist_eu, "average")
hc_eu_ward <- hclust(dist_eu, "ward.D2")

par(mfrow=c(1,2))

plot(hc_eu_average, main = "Euclidean + Average")
cutree(hc_eu_average, k = 8)
rect.hclust(hc_eu_average, k = 8)

plot(hc_eu_ward, main = "Euclidean + Ward")
cutree(hc_eu_ward, k = 3)
rect.hclust(hc_eu_ward, k = 3)


dist_man <- dist(dom_sp_compo, "manhattan")

hc_man_average <- hclust(dist_man, "average")
hc_man_ward <- hclust(dist_man, "ward.D2")

par(mfrow=c(1,2))

plot(hc_man_average, main = "Manhattan + Average")
cutree(hc_man_average, k = 3)
rect.hclust(hc_man_average, k = 3)

plot(hc_man_ward, main = "Manhattan + Ward")
cutree(hc_man_ward, k = 3)
rect.hclust(hc_man_ward, k = 3)


dist_bc <- vegdist(dom_sp_compo, "bray")

hc_bc_average <- hclust(dist_bc, "average")
hc_bc_ward <- hclust(dist_bc, "ward.D2")

par(mfrow=c(1,2))

plot(hc_bc_average, main = "Bray-Curtis + Average")
cutree(hc_bc_average, k = 3)
rect.hclust(hc_bc_average, k = 3)

plot(hc_bc_ward, main = "Bray-Curtis + Ward")
cutree(hc_bc_ward, k = 3)
rect.hclust(hc_bc_ward, k = 3)


dist_ja <- vegdist(dom_sp_compo, "jaccard", binary = T)

hc_ja_average <- hclust(dist_ja, "average")
hc_ja_complete <- hclust(dist_ja, "complete")

par(mfrow=c(1,2))

plot(hc_ja_average, main = "Jaccard + Average")
cutree(hc_ja_average, k = 4)
rect.hclust(hc_ja_average, k = 4)

plot(hc_ja_complete, main = "Jaccard + Complete")
cutree(hc_ja_complete, k = 4)
rect.hclust(hc_ja_complete, k = 4)

distance <- c("Euclidean", "Manhattan", "Bray-Curtis", "Jaccard")
k_average <- c(8, 3, 3, 4)
k_ward <- c(3, 3, 3, NA)
k_complete <- c(NA, NA, NA, 4)
k_compare <- data.frame(distance, k_average, k_ward, k_complete)
k_compare

#cophenetic correlation

coph_eu_average <- round(cor(cophenetic(hc_eu_average), dist_eu), 3)
coph_eu_ward <- round(cor(cophenetic(hc_eu_ward), dist_eu), 3)
coph_man_average <- round(cor(cophenetic(hc_man_average), dist_man), 3)
coph_man_ward <- round(cor(cophenetic(hc_man_ward), dist_man), 3)
coph_bc_average <- round(cor(cophenetic(hc_bc_average), dist_bc), 3)
coph_bc_ward <- round(cor(cophenetic(hc_bc_ward), dist_bc), 3)
coph_ja_average <- round(cor(cophenetic(hc_ja_average), dist_ja), 3)
coph_ja_complete <- round(cor(cophenetic(hc_ja_complete), dist_ja), 3)

distance <- c("Euclidean", "Manhattan", "Bray-Curtis", "Jaccard")
coph_average <- c(coph_eu_average, coph_man_average, coph_bc_average, coph_ja_average)
coph_ward <- c(coph_eu_ward, coph_man_ward, coph_bc_ward, NA)
coph_complete <- c(NA, NA, NA, coph_ja_complete)
cop_compare <- data.frame(distance, coph_average, coph_ward, coph_complete)
cop_compare

#non-hierarchical

sil_eu <- sapply(2:8, function(k){
  pam_eu <- pam(dist_eu, k)
  pam_eu$silinfo$avg.width
})

par(mfrow=c(1,1))
plot(2:8, sil_eu, type = "b",
     xlab = "k",
     ylab = "Average Silhouette Width",
     main = "Euclidean")

pam_eu_best <- pam(dist_eu, k = 3)
pam_eu_best$clustering

plot(silhouette(pam_eu_best), border = NA, main = "k = 3")

#man

sil_man <- sapply(2:8, function(k){
  pam_man <- pam(dist_man, k)
  pam_man$silinfo$avg.width
})

par(mfrow=c(1,1))
plot(2:8, sil_man, type = "b",
     xlab = "k",
     ylab = "Average Silhouette Width",
     main = "Manhattan")

pam_man_best <- pam(dist_man, k = 3)
pam_man_best$clustering

plot(silhouette(pam_man_best), border = NA, main = "k = 3")

#bc

sil_bc <- sapply(2:8, function(k){
  pam_bc <- pam(dist_bc, k)
  pam_bc$silinfo$avg.width
})

par(mfrow=c(1,1))
plot(2:8, sil_bc, type = "b",
     xlab = "k",
     ylab = "Average Silhouette Width",
     main = "Bray-Curtis")

pam_bc_best <- pam(dist_bc, k = 3)
pam_bc_best$clustering

plot(silhouette(pam_bc_best), border = NA, main = "k = 3")

#ja

sil_ja <- sapply(2:8, function(k){
  pam_ja <- pam(dist_ja, k)
  pam_ja$silinfo$avg.width
})

par(mfrow=c(1,1))
plot(2:8, sil_ja, type = "b",
     xlab = "k",
     ylab = "Average Silhouette Width",
     main = "Jaccard")

pam_ja_best <- pam(dist_ja, k = 2)
pam_ja_best$clustering

plot(silhouette(pam_ja_best), border = NA, main = "k = 2")

distance <- c("Euclidean", "Manhattan", "Bray-Curtis", "Jaccard")
best_k <- c(3, 3, 3, 2)
average_silhouette <- c(0.27, 0.28, 0.29, 0.44)
compare <- data.frame(distance, best_k, average_silhouette)
compare
