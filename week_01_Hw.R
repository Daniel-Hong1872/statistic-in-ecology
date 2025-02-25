setwd("C:/Users/dan91/Rstudio/stat_data/mock_dataTS")
cop_compo <- read.table("Copepod_composition.txt", header = T)
cop_den <- read.table("Cop_density.txt")
cop_sp <- read.table("copepodSPlist.txt", fill = T, sep = "\n") #seperate by different row


# 1.Calculate the copepod density for each species for each cruise-station

colnames(cop_compo) <- c("1", "3", "4", "6", "13", "16", "19a", "21", "23a", "25a", "18", "19b",
                         "20", "22a", "23b", "25b", "27a", "29a", "Aa", "Ba", "Ca", "Da", "E",
                         "F", "G", "22b", "23c", "25c", "27b", "29b", "Ab", "Bb", "Cb", "Db")
head(cop_compo)

rownames(cop_den) <- colnames(cop_compo)
head(cop_den)

sp_den_cr <- matrix(0, nrow = length(cop_sp$V1), ncol = length(cop_den$V1))
for(i in 1:length(cop_sp$V1)){
  for(j in 1:length(cop_den$V1)){
    sp_den_cr[i, j] <- cop_compo[i, j]/100 * cop_den[j,]
  }
}
sp_den_cr.df <- as.data.frame(sp_den_cr)
colnames(sp_den_cr.df) <- colnames(cop_compo)
head(sp_den_cr.df)

sp_den_cr.df$`19a` <- sp_den_cr.df$`19a` + sp_den_cr.df$`19b`

sp_den_cr.df$`23a` <- sp_den_cr.df$`23a` + sp_den_cr.df$`23b` + sp_den_cr.df$`23c`

sp_den_cr.df$`25a` <- sp_den_cr.df$`25a` + sp_den_cr.df$`25b` + sp_den_cr.df$`25c`

sp_den_cr.df$`22a` <- sp_den_cr.df$`22a` + sp_den_cr.df$`22b`

sp_den_cr.df$`27a` <- sp_den_cr.df$`27a` + sp_den_cr.df$`27b`

sp_den_cr.df$`29a` <- sp_den_cr.df$`29a` + sp_den_cr.df$`29b`

sp_den_cr.df$Aa <- sp_den_cr.df$Aa + sp_den_cr.df$Ab

sp_den_cr.df$Ba <- sp_den_cr.df$Ba + sp_den_cr.df$Bb

sp_den_cr.df$Ca <- sp_den_cr.df$Ca + sp_den_cr.df$Cb

sp_den_cr.df$Da <- sp_den_cr.df$Da + sp_den_cr.df$Db

sp_den_cr.df <- sp_den_cr.df[, -c(12, 15, 16, 26, 27, 28, 29, 30, 31,32,33,34)]
colnames(sp_den_cr.df) <- c("1", "3", "4", "6", "13", "16", "19", "21", "23", "25", "18", 
                            "20", "22", "27", "29", "A", "B", "C", "D", "E", "F", "G")
rownames(sp_den_cr.df) <- cop_sp$V1
head(sp_den_cr.df)

# 2.For each cruise-station, calculate the species richness (number of species) and Shannon diversity index
# 2(1). species richness 

sp_rich <- ifelse(sp_den_cr.df > 0, 1, 0) 
sp_rich.df <- as.data.frame(sp_rich)

richness <- data.frame(colSums(sp_rich.df))
head(richness)

# 2(2). Shannon diversity index

cop_compo$`19a` <- cop_compo$`19a` + cop_compo$`19b`

cop_compo$`23a` <- cop_compo$`23a` + cop_compo$`23b` + cop_compo$`23c`

cop_compo$`25a` <- cop_compo$`25a` + cop_compo$`25b` + cop_compo$`25c`

cop_compo$`22a` <- cop_compo$`22a` + cop_compo$`22b`

cop_compo$`27a` <- cop_compo$`27a` + cop_compo$`27b`

cop_compo$`29a` <- cop_compo$`29a` + cop_compo$`29b`

cop_compo$Aa <- cop_compo$Aa + cop_compo$Ab

cop_compo$Ba <- cop_compo$Ba + cop_compo$Bb

cop_compo$Ca <- cop_compo$Ca + cop_compo$Cb

cop_compo$Da <- cop_compo$Da + cop_compo$Db

cop_compo <- cop_compo[, -c(12, 15, 16, 26, 27, 28, 29, 30, 31,32,33,34)]
colnames(cop_compo) <- c("1", "3", "4", "6", "13", "16", "19", "21", "23", "25", "18", 
                            "20", "22", "27", "29", "A", "B", "C", "D", "E", "F", "G")
rownames(cop_compo) <- cop_sp$V1
head(cop_compo)

shannon_index <- data.frame(-colSums(cop_compo * log(cop_compo / 100), na.rm = T))
colnames(shannon_index) <- "Shannon index"
head(shannon_index)

# 3. For each of the dominant species (species >=5% of total composition in any cruise-station regardless of the seasons), calculate the average density by seasons.

cop_compo_season <- read.table("Copepod_composition.txt", header = T)
sp_den_cr_season <- matrix(0, nrow = length(cop_sp$V1), ncol = length(cop_den$V1))
for(i in 1:length(cop_sp$V1)){
  for(j in 1:length(cop_den$V1)){
    sp_den_cr_season[i, j] <- cop_compo_season[i, j]/100 * cop_den[j,]
  }
}
sp_den_cr_season.df <- as.data.frame(sp_den_cr_season)
colnames(sp_den_cr_season.df) <- colnames(cop_compo_season)
head(sp_den_cr_season.df)

dominant_sp.df <- data.frame(ifelse(cop_compo_season >= 5, 1, 0))
dominant_species <- rownames(dominant_sp.df)[rowSums(dominant_sp.df) >= 1]
dominant_species
 
dominant_sp_den <- sp_den_cr_season.df[dominant_species, ]
dominant_sp_den$spring_mean <- rowMeans(dominant_sp_den[, 1:10])
dominant_sp_den$summer_mean <- rowMeans(dominant_sp_den[, 11:25])
dominant_sp_den$winter_mean <- rowMeans(dominant_sp_den[, 26:34])

dominant_sp_den_season <- dominant_sp_den[, 35:37]
