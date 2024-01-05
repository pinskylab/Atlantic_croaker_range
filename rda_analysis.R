##RDA Analysis##

#Using environmental data from Copernicus and filtered genomic data with ~23K SNPs

##Individual frequencies using all filters applied dataset (~23K SNPs)

#Load in genomic data
indfreq_filters <- np$load("pcangsd_minmaf0.05_bam_list_realigned_mindp132_maxdp4000_minind0_minmaf0.indf.npy")
#400 individuals, 22,548 SNPs

indfreq_filters <- as.data.frame(indfreq_filters)
indfreq_filters[1:5,1:5] #each individual has a row and each SNP has a column
#Need to transfer over the contig/position numbers for each SNP
sum(is.na(indfreq_filters)) #No NAs, so no imputation needed

#Load in mafs file from ANGSD to get contig/position numbers for each SNP
mafs_filters <- read.table("bam_list_realigned_mindp132_maxdp4000_minind0_minmaf0.05.mafs.txt", header = T)

#Combine chromo and position columns
mafs_filters$chromo_position <- paste(mafs_filters$chromo, mafs_filters$position, sep="_")
snp_position_filters <- as.vector(mafs_filters$chromo_position)

colnames(indfreq_filters) <- snp_position_filters
indfreq_filters[1:5, 1:5] #Now each column has the contig and position number as the header


#pRDA

## Running a PCA on neutral genetic markers
pca_filters <- rda(indfreq_filters, scale=T) # PCA in vegan uses the rda() call without any predictors

#Screeplot of the PCA eigenvalues:

screeplot(pca_filters, type = "barplot", npcs=10, main="PCA Eigenvalues")

#Based on the screeplot, one PC would be a reasonable set to retain as a proxy for
#neutral population structure in downstream analyses.

## Neutral population structure table
PC_filters <- scores(pca_filters, choices=1, display="sites", scaling=0)
PopStruct_filters <- data.frame(Individual = seq(1:400), PC_filters)
colnames(PopStruct_filters) <- c("Individual", "PC1")
PopStruct_filters[1:5,]

#Add PC1 to the environmental variables data frame
env_ocean_filters <- subset(env_oct, select = c(lat, depth, oxygen, temp, salinity))
env_ocean_filters$pc1 <- PopStruct_filters$PC1

#Make sure variables are numeric
env_ocean_filters$oxygen <- as.numeric(env_ocean_filters$oxygen)
env_ocean_filters$temp <- as.numeric(env_ocean_filters$temp)
env_ocean_filters$salinity <- as.numeric(env_ocean_filters$salinity)

#pRDA accounting for population structure using PC1
#Using latitude, depth, oxygen, temperature, and salinity as explanatory variables and PC1 as a conditioning variable
#to account for neutral population structure.
library(vegan) #version 2.6-4
rda_ocean_filters <- rda(indfreq_filters ~ lat + depth + oxygen + temp + salinity + Condition(pc1),  data=env_ocean_filters, scale=T)
rda_ocean_filters
#                Inertia Proportion Rank
#Total         2.252e+04  1.000e+00
#Conditional   2.250e+04  9.993e-01    1
#Constrained   3.450e-01  1.532e-05    5
#Unconstrained 1.543e+01  6.853e-04  163
#Inertia is correlations

#Calculate adjusted R-squared
RsquareAdj(rda_ocean_filters)
#r.squared=1.532103e-05 adj.r.squared=6.618678e-06

summary(eigenvals(rda_ocean_filters, model = "constrained"))
#Most variance explained by RDA1 and RDA2, some on RDA3
#Proportion explained RDA1: 0.64 RDA2: 0.21319 RDA3: 0.10983 RDA4:0.020442 RDA5:0.016529

screeplot(rda_ocean_filters) #Looks reasonable to include RDA1,2,and 3 for SNP candidates

plot(rda_ocean_filters)

load.rda.ocean.filters <- scores(rda_ocean_filters, choices=c(1:3), display="species")

hist(load.rda.ocean.filters[,1], main="Loadings on RDA1")
hist(load.rda.ocean.filters[,2], main="Loadings on RDA2")
hist(load.rda.ocean.filters[,3], main="Loadings on RDA3")

cand1_ocean_filters <- outliers(load.rda.ocean.filters[,1],3) #447
cand2_ocean_filters <- outliers(load.rda.ocean.filters[,2],3) #20
cand3_ocean_filters <- outliers(load.rda.ocean.filters[,3],3) #331

ncand_ocean_filters <- length(cand1_ocean_filters) + length(cand2_ocean_filters) + length(cand3_ocean_filters)
ncand_ocean_filters #798

#Organize results by making one data frame with the axis, SNP name, loading, and correlation with each predictor

cand1_ocean_filters <- cbind.data.frame(rep(1,times=length(cand1_ocean_filters)), names(cand1_ocean_filters), unname(cand1_ocean_filters))
cand2_ocean_filters <- cbind.data.frame(rep(2,times=length(cand2_ocean_filters)), names(cand2_ocean_filters), unname(cand2_ocean_filters))
cand3_ocean_filters <- cbind.data.frame(rep(3,times=length(cand3_ocean_filters)), names(cand3_ocean_filters), unname(cand3_ocean_filters))


colnames(cand1_ocean_filters) <- colnames(cand2_ocean_filters) <- colnames(cand3_ocean_filters)  <- c("axis","snp","loading")

cand_ocean_filters <- rbind(cand1_ocean_filters, cand2_ocean_filters, cand3_ocean_filters)
cand_ocean_filters$snp <- as.character(cand_ocean_filters$snp)

#Add in correlations of each candidate SNP with the four environmental predictors:

foo <- matrix(nrow=(ncand_ocean_filters), ncol=5)  # 5 columns for 5 predictors
colnames(foo) <- c("lat","depth","oxygen","temp", "salinity")


env_ocean_filters_cand <- subset(env_ocean_filters, select = c(lat, depth, oxygen, temp, salinity))

for (i in 1:length(cand_ocean_filters$snp)) {
  nam <- cand_ocean_filters[i,2]
  snp.gen <- indfreq_filters[,nam]
  foo[i,] <- apply(env_ocean_filters_cand,2,function(x) cor(x,snp.gen))
}

cand_ocean_filters <- cbind.data.frame(cand_ocean_filters,foo)
head(cand_ocean_filters)

#Investigate the SNP candidates
#Look for duplicate detections (SNPs identified as candidates on more than one RDA axis)

length(cand_ocean_filters$snp[duplicated(cand_ocean_filters$snp)]) #202 duplicate detections

foo <- cbind(cand_ocean_filters$axis, duplicated(cand_ocean_filters$snp))
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) #13 duplicates on axis 2
table(foo[foo[,1]==3,2]) #189 duplicates on axis 3

cand_ocean_filters <- cand_ocean_filters[!duplicated(cand_ocean_filters$snp),] # remove duplicate detections
#now have 596 unique candidate SNPs

#See which of the predictors each SNP is most strongly correlated with
for (i in 1:length(cand_ocean_filters$snp)) {
  bar <- cand_ocean_filters[i,]
  cand_ocean_filters[i,9] <- names(which.max(abs(bar[4:8]))) # gives the variable
  cand_ocean_filters[i,10] <- max(abs(bar[4:8]))              # gives the correlation
}

colnames(cand_ocean_filters)[9] <- "predictor"
colnames(cand_ocean_filters)[10] <- "correlation"

table(cand_ocean_filters$predictor) #depth=8, oxygen=583, salinity=2, temp=3
#oxygen is closely followed by salinity and temperature (but with inverse correlations,
#ie oxygen is negative and then salinity and temperature are positive)

#Identify how many unique contigs the SNPs are on
#Separate the snp column into contig and position
library(dplyr)
library(tidyr)
cand_ocean_filters_contig <- cand_ocean_filters %>% separate_wider_delim(snp, "_", names=c('contig',  'position'))

cand_ocean_filters$contig <- cand_ocean_filters_contig$contig
cand_ocean_filters$position <- cand_ocean_filters_contig$position

unique_ocean_filters_contigs <- unique(cand_ocean_filters[c("contig")]) #130 unique contigs

##Plotting SNPs

#Formatting table for ggplot
locus_scores_ocean_filters <- scores(rda_ocean_filters, choices=c(1:2), display="species", scaling="none")
TAB_loci_ocean_filters <- data.frame(names = row.names(locus_scores_ocean_filters), locus_scores_ocean_filters)
TAB_loci_ocean_filters$type <- "Neutral"
TAB_loci_ocean_filters$type[TAB_loci_ocean_filters$names%in%cand_ocean_filters$snp] <- "Outlier"
TAB_loci_ocean_filters$type <- factor(TAB_loci_ocean_filters$type, levels = c("Neutral", "Outlier"))
TAB_loci_ocean_filters <- TAB_loci_ocean_filters[order(TAB_loci_ocean_filters$type),]
TAB_var_ocean_filters <- as.data.frame(scores(rda_ocean_filters, choices=c(1,2), display="bp")) # pull the biplot scores

## Biplot of RDA loci and variables scores
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci_ocean_filters, aes(x=RDA1, y=RDA2, colour = type), shape=19, size = 2) +
  scale_color_manual(values = c("gray60", "#F9A242FF")) +
  #geom_segment(data = TAB_var_oct_filters, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.2, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  #geom_text(data = TAB_var_oct_filters, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var_oct_filters)), size = 3.75, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  #facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(axis.text.x= element_text(size=11), axis.text.y=element_text(size=12)) +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=11))

#Plot the SNPs by predictor
#Subset based on predictor
temp_ocean_filters <- subset(cand_ocean_filters, predictor=='temp', select=c(snp, predictor))
salinity_ocean_filters <- subset(cand_ocean_filters, predictor=='salinity', select=c(snp, predictor))
depth_ocean_filters <- subset(cand_ocean_filters, predictor=='depth', select=c(snp, predictor))
oxygen_ocean_filters <- subset(cand_ocean_filters, predictor=='oxygen', select=c(snp, predictor))

TAB_loci_ocean_filters$predictor <- "Neutral"
TAB_loci_ocean_filters$predictor[TAB_loci_ocean_filters$names%in%temp_ocean_filters$snp] <- "Temperature"
TAB_loci_ocean_filters$predictor[TAB_loci_ocean_filters$names%in%salinity_ocean_filters$snp] <- "Salinity"
TAB_loci_ocean_filters$predictor[TAB_loci_ocean_filters$names%in%depth_ocean_filters$snp] <- "Depth"
TAB_loci_ocean_filters$predictor[TAB_loci_ocean_filters$names%in%oxygen_ocean_filters$snp] <- "Oxygen"
TAB_loci_ocean_filters$predictor <- factor(TAB_loci_ocean_filters$predictor, levels = c("Neutral", "Temperature", "Salinity", "Depth", "Oxygen"))
TAB_loci_ocean_filters <- TAB_loci_ocean_filters[order(TAB_loci_ocean_filters$predictor),]

ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci_ocean_filters, aes(x=RDA1, y=RDA2, colour = predictor), shape=19, size = 1.75) +
  scale_color_manual(values = c("gray60", '#e31a1c', '#1f78b4','#6a3d9a', '#F9A242FF')) +
  geom_segment(data = TAB_var_ocean_filters, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.2, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var_ocean_filters, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var_ocean_filters)), size = 3, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  #facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(axis.text.x= element_text(size=11), axis.text.y=element_text(size=12)) +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=11))
