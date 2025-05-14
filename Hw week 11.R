#dominant species
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

dom_sp_compo <- t(cop_compo[dominant_species, ])
colnames(dom_sp_compo) <- cop_sp[dominant_species, ]
rownames(dom_sp_compo) <- colnames(cop_compo)
dom_sp_compo.df <- as.data.frame(dom_sp_compo)
rownames(dom_sp_compo.df)[1] <- "p1"

#environmental data
library(readxl)
envi_data <- read_excel("enviANDdensity.xls")
envi.df <- as.data.frame(envi_data)
rownames(envi.df) <- envi_data$station
envi.df$station <- NULL
rownames(envi.df)[11:34] <- c("s18", "s19", "s20", "s22", "s23", "s25", "s27", "s29", "sA",
                              "sB", "sC", "sD", "sE", "sF", "sG", "w22", "w23", "w25", "w27",
                              "w29", "wA", "wB", "wC", "wD")

#Apply constrained ordination
library(vegan)
dca_y <- decorana(dom_sp_compo.df)
summary(dca_y)

cop_hel <- decostand(dom_sp_compo.df, method = "hellinger")
rda_y <- rda(cop_hel ~ Depth + Temperature + Salinity + Fluorescence +
               DissolvedOxygen, data = envi.df)
cca_y <- cca(dom_sp_compo.df ~ Depth + Temperature + Salinity + Fluorescence +
               DissolvedOxygen, data = envi.df)

summary(rda_y)
summary(cca_y)

vif.cca(rda_y)
vif.cca(cca_y)

anova(rda_y, by = "term", permutations = 999)
anova(cca_y, by = "term", permutations = 999)

rda_noDO <- rda(cop_hel ~ Depth + Temperature + Salinity + Fluorescence
                , data = envi.df)
cca_noDO <- cca(dom_sp_compo.df ~ Depth + Temperature + Salinity + Fluorescence
                , data = envi.df)

anova(rda_noDO, by = "term", permutations = 999)
anova(cca_noDO, by = "term", permutations = 999)

plot(rda_noDO, scaling = 1, main = "RDA (scaling 1)")
plot(rda_noDO, scaling = 2, main = "RDA (scaling 2)")
plot(cca_noDO, scaling = 1, main = "CCA (scaling 1)")
plot(cca_noDO, scaling = 2, main = "CCA (scaling 2)")


