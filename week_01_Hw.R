cop_compo <- read.table("Copepod_composition.txt", header = T)
cop_den <- read.table("Cop_density.txt")
cop_sp <- read.table("copepodSPlist.txt", fill = T, sep = "\n") #seperate by different row


# 1.Calculate the copepod density for each species for each cruise-station
sp_den_cr <- matrix(0, nrow = length(cop_sp$V1), ncol = length(cop_den$V1))
for(i in 1:length(cop_sp$V1)){
  for(j in 1:length(cop_den$V1)){
    sp_den_cr[i, j] <- cop_compo[i, j]/100 * cop_den[j,]
  }
}
sp_den_cr.df <- as.data.frame(sp_den_cr)
sp_den_cr.df

# 2.For each cruise-station, calculate the species richness (number of species) and Shannon diversity index
 