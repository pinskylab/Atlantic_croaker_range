##RDA Analysis##

#Using environmental data from Copernicus and filtered genomic data with ~23K SNPs
#Following tutorial by Brenna Forester (https://popgen.nescent.org/2018-03-27_RDA_GEA.html
#and https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html)

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

#Get summary data for manuscript table
sum(env_oct$population=='GA')

env_oct %>% group_by(population) %>% summarize(mean_value = mean(lat)) %>% as.data.frame()

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

##Checking how this worked
cand_ocean_filters_trial <- rbind(cand1_ocean_filters, cand2_ocean_filters, cand3_ocean_filters)
cand_ocean_filters_trial$snp <- as.character(cand_ocean_filters_trial$snp)

foo <- matrix(nrow=(ncand_ocean_filters), ncol=5)  # 5 columns for 5 predictors
colnames(foo) <- c("lat","depth","oxygen","temp", "salinity")

env_ocean_filters_cand <- subset(env_ocean_filters, select = c(lat, depth, oxygen, temp, salinity))

for (i in 1:length(cand_ocean_filters_trial$snp)) {
  nam <- cand_ocean_filters_trial[i,2]
  snp.gen <- indfreq_filters[,nam]
  foo[i,] <- apply(env_ocean_filters_cand,2,function(x) cor(x,snp.gen))
}

cand_ocean_filters <- cbind.data.frame(cand_ocean_filters,foo)
head(cand_ocean_filters)

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

#Formatting table for ggplot
TAB_var_ocean_filters <- as.data.frame(scores(rda_ocean_filters, choices=c(1,2), display="bp")) # pull the biplot scores

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


##Identifying outliers using the rdaadapt function

## Function rdadapt
source("./rdadapt.R")

library(robust)
library(qvalue)

## Running the function with K = 2
rdadapt_ocean <- rdadapt(rda_ocean_filters, 2)

##Add in q-values
rdadapt_ocean$q.values <- p.adjust(rdadapt_ocean$p.values, method="fdr")
sum(rdadapt_ocean$q.values<0.05) #2320 SNPs

## Identifying the loci that are below the q-value threshold
outliers_ocean <- data.frame(Loci = colnames(indfreq_filters)[which(rdadapt_ocean$q.values<0.05)], q.value =
                             rdadapt_ocean$q.values[which(rdadapt_ocean$q.values<0.05)], contig =
                             unlist(lapply(strsplit(colnames(indfreq_filters)[which(rdadapt_ocean$q.values<0.05)],
                                                    split = "_"), function(x) x[1])))
## Top hit outlier per contig
outliers_ocean <- outliers_ocean[order(outliers_ocean$contig, outliers_ocean$q.value),]

## List of outlier names
outliers_rdadapt_ocean <- as.character(outliers_ocean$Loci[!duplicated(outliers_ocean$contig)]) #228 contigs

#Add in correlations of each candidate SNP with the four environmental predictors:

foo <- matrix(nrow=(2320), ncol=5)  # 5 columns for 5 predictors
colnames(foo) <- c("lat","depth","oxygen","temp", "salinity")

env_ocean_filters_cand <- subset(env_ocean_filters, select = c(lat, depth, oxygen, temp, salinity))
outliers_ocean2 <- outliers_ocean #copy just in case

for (i in 1:length(outliers_ocean$Loci)) {
  nam2 <- outliers_ocean[i,1]
  snp.gen2 <- indfreq_filters[,nam2]
  foo[i,] <- apply(env_ocean_filters_cand,2,function(x) cor(x,snp.gen2))
}

full_outliers_ocean <- cbind.data.frame(outliers_ocean,foo)
head(full_outliers_ocean)

#See which of the predictors each SNP is most strongly correlated with
for (i in 1:length(full_outliers_ocean$Loci)) {
  bar <- full_outliers_ocean[i,]
  full_outliers_ocean[i,9] <- names(which.max(abs(bar[4:8]))) # gives the variable
  full_outliers_ocean[i,10] <- max(abs(bar[4:8]))              # gives the correlation
}

colnames(full_outliers_ocean)[9] <- "predictor"
colnames(full_outliers_ocean)[10] <- "correlation"

table(full_outliers_ocean$predictor) #depth 8, oxygen 2307, salinity 2, temp 3
#oxygen is closely followed by salinity and temperature (but with inverse correlations,
#ie oxygen is negative and then salinity and temperature are positive)
#Using only the top hit per contig, all 228 are oxygen


#Once the outliers have been identified, it can be useful to visualize their distribution in
#comparison with neutral loci using either an RDA biplot or a Manhattan plot.

## Formatting table for ggplot
locus_scores_ocean <- scores(rda_ocean_filters, choices=c(1:2), display="species", scaling="none")
locus_scores_ocean <- scores(rdadapt_ocean, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci_ocean <- data.frame(names = row.names(locus_scores_ocean), locus_scores_ocean)
TAB_loci_ocean$type <- "Neutral"
TAB_loci_ocean$type[TAB_loci_ocean$names%in%outliers_ocean$Loci] <- "All outliers"
TAB_loci_ocean$type[TAB_loci_ocean$names%in%outliers_rdadapt_ocean] <- "Top outliers"
TAB_loci_ocean$type <- factor(TAB_loci_ocean$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci_ocean <- TAB_loci_ocean[order(TAB_loci_ocean$type),]
TAB_var_ocean <- as.data.frame(scores(rdadapt_ocean, choices=c(1:2), display="species")) # pull the biplot scores
TAB_var_ocean <- as.data.frame(scores(rda_ocean_filters, choices=c(1,2), display="bp"))

rownames(TAB_var_ocean_filters) <- c("Latitude", "Depth", "Oxygen", "Temperature", "Salinity")
## Biplot of RDA loci and variables scores
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci_ocean, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var_ocean_filters, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var_ocean_filters, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var_ocean_filters)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

## Manhattan plot
Outliers_ocean <- rep("Neutral", length(colnames(indfreq_filters)))
Outliers_ocean[colnames(indfreq_filters)%in%outliers_ocean$Loci] <- "All outliers"
Outliers_ocean[colnames(indfreq_filters)%in%outliers_rdadapt_ocean] <- "Top outliers"
Outliers_ocean <- factor(Outliers_ocean, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_manhattan_ocean <- data.frame(pos = 1:length(colnames(indfreq_filters)),
                                  Loci=colnames(indfreq_filters),
                               qvalues = rdadapt_ocean$q.values,
                               Outliers = Outliers_ocean)
TAB_manhattan_ocean <- TAB_manhattan_ocean[order(TAB_manhattan_ocean$Outliers),]

ggplot(data = TAB_manhattan_ocean) +
  geom_point(aes(x=pos, y=-log10(qvalues), col = Outliers), size=1.4, shape=19, alpha=0.5) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  xlab("Loci") + ylab("-log10(q.values)") +
  #geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  #facet_wrap(~"Manhattan plot", nrow = 3) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(),
        legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(),
        legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11)) +
  ylim(0, 1000)

#Replace 0 q-value with 1e-300 to avoid the -log10 being infinity
TAB_manhattan_ocean2 <- TAB_manhattan_ocean
TAB_manhattan_ocean2[TAB_manhattan_ocean2==0.000000e+00] <- 1e-300
TAB_manhattan_ocean2[TAB_manhattan_ocean2==1e-300] <- 1e-320

TAB_manhattan_ocean3 <- TAB_manhattan_ocean
TAB_manhattan_ocean3[TAB_manhattan_ocean3==0.000000e+00] <- 1e-320
TAB_manhattan_ocean3$Outliers <- as.character(TAB_manhattan_ocean3$Outliers)
TAB_manhattan_ocean3[TAB_manhattan_ocean3=="All outliers"] <- "Outlier"
TAB_manhattan_ocean3[TAB_manhattan_ocean3=="Top outliers"] <- "Outlier"


ggplot(data = TAB_manhattan_ocean3) +
  geom_point(aes(x=pos, y=-log10(qvalues), col = Outliers), size=1.4, shape=19, alpha=0.75) +
  scale_color_manual(values = c("gray90", "#F9A242FF")) +
  xlab("Loci") + ylab("-log10(q.values)") +
  #geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  #facet_wrap(~"Manhattan plot", nrow = 3) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 12, base_family = "Times") +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(),
        legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(),
        legend.text=element_text(size=rel(.8)), strip.text = element_text(size=12))

TAB_manhattan_ocean2$Outliers <- as.character(TAB_manhattan_ocean2$Outliers)
TAB_manhattan_ocean2[TAB_manhattan_ocean2=="All outliers"] <- "Outlier"
TAB_manhattan_ocean2[TAB_manhattan_ocean2==NA] <- "Outlier"
TAB_manhattan_ocean2[TAB_manhattan_ocean2=="Top outliers"] <- "Outlier"

#Plot the SNPs by predictor
#Subset based on predictor
temp_ocean_outliers <- subset(full_outliers_ocean, predictor=='temp', select=c(Loci, predictor))
salinity_ocean_outliers <- subset(full_outliers_ocean, predictor=='salinity', select=c(Loci, predictor))
depth_ocean_outliers <- subset(full_outliers_ocean, predictor=='depth', select=c(Loci, predictor))
oxygen_ocean_outliers <- subset(full_outliers_ocean, predictor=='oxygen', select=c(Loci, predictor))

TAB_loci_ocean$predictor <- "Neutral"
TAB_loci_ocean$predictor[TAB_loci_ocean$names%in%temp_ocean_outliers$Loci] <- "Temperature"
TAB_loci_ocean$predictor[TAB_loci_ocean$names%in%salinity_ocean_outliers$Loci] <- "Salinity"
TAB_loci_ocean$predictor[TAB_loci_ocean$names%in%depth_ocean_outliers$Loci] <- "Depth"
TAB_loci_ocean$predictor[TAB_loci_ocean$names%in%oxygen_ocean_outliers$Loci] <- "Oxygen"
TAB_loci_ocean$predictor <- factor(TAB_loci_ocean$predictor, levels = c("Neutral", "Temperature", "Salinity", "Depth", "Oxygen"))
TAB_loci_ocean <- TAB_loci_ocean[order(TAB_loci_ocean$predictor),]

ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci_ocean, aes(x=RDA1, y=RDA2, colour = predictor), shape=1, size = 1.75) +
  scale_color_manual(values = c("gray60", '#e31a1c', '#1f78b4','#6a3d9a', '#F9A242FF')) +
  geom_segment(data = TAB_var_ocean, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.2, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var_ocean, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var_ocean_filters)), size = 3, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  #facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(axis.text.x= element_text(size=11), axis.text.y=element_text(size=12)) +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=11))

#Manhattan plot with SNPs by predictor
TAB_manhattan_ocean$predictor <- "Neutral"
TAB_manhattan_ocean$predictor[TAB_manhattan_ocean$Loci%in%temp_ocean_outliers$Loci] <- "Temperature"
TAB_manhattan_ocean$predictor[TAB_manhattan_ocean$Loci%in%salinity_ocean_outliers$Loci] <- "Salinity"
TAB_manhattan_ocean$predictor[TAB_manhattan_ocean$Loci%in%depth_ocean_outliers$Loci] <- "Depth"
TAB_manhattan_ocean$predictor[TAB_manhattan_ocean$Loci%in%oxygen_ocean_outliers$Loci] <- "Oxygen"
TAB_manhattan_ocean$predictor <- factor(TAB_manhattan_ocean$predictor, levels = c("Neutral", "Oxygen", "Temperature", "Salinity", "Depth"))
TAB_manhattan_ocean <- TAB_manhattan_ocean[order(TAB_manhattan_ocean$predictor),]

ggplot(data = TAB_manhattan_ocean) +
  geom_point(aes(x=pos, y=-log10(qvalues), col = predictor), size=1.4, shape=19) +
  scale_color_manual(values = c("gray60", '#F9A242FF', '#e31a1c', '#1f78b4','#6a3d9a')) +
  xlab("Loci") + ylab("-log10(q.values)") +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  #facet_wrap(~"Manhattan plot", nrow = 3) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(),
        legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(),
        legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

##Check variance inflation factors
library(vegan)
vif.cca(rda_ocean_filters)
#pc1       lat     depth    oxygen      temp  salinity
#1.022314 30.264763  2.619818 22.474251 27.576771  9.727815
#a little high for lat, oxygen, and temp but to be expected due to their physical properties/correlations

##Create vcf file for SnpEff

vcf_2320 <- as.data.frame(full_outliers_ocean$contig)
colnames(vcf_2320)[1] <- "#CHROM"
vcf_2320$POS <- full_outliers_ocean$Loci
vcf_2320$ID <- full_outliers_ocean$Loci
vcf_2320$REF <- "NA"
vcf_2320$ALT <- "NA"
vcf_2320$QUAL <- "."
vcf_2320$FILTER <- "PASS"
vcf_2320$INFO <- "."
vcf_2320<- vcf_2320 %>% separate(POS, c('contig', 'pos'))
colnames(vcf_2320)[3] <- 'POS'
vcf_2320 <- vcf_2320[,-2]

#Add ref and alt information
mafs_filters$POS <- mafs_filters$position
mafs_filters$REF <- mafs_filters$major
mafs_filters$ALT <- mafs_filters$minor
mafs_filters_tojoin <- matrix(ncol=3, nrow=22548)
mafs_filters_tojoin <- as.data.frame(mafs_filters_tojoin)
mafs_filters_tojoin$V1 <- mafs_filters$POS
mafs_filters_tojoin$V2 <- mafs_filters$REF
mafs_filters_tojoin$V3 <- mafs_filters$ALT
colnames(mafs_filters_tojoin) <- c("POS", "REF", "ALT")
vcf_2320$POS <- as.integer(vcf_2320$POS)
vcf_2320_join <- full_join(vcf_2320, mafs_filters_tojoin, by="POS")

#Issue with position column in the join, so will go with ID to make sure it is joined correctly
mafs_filters_tojoin$POS <- mafs_filters$chromo_position
colnames(mafs_filters_tojoin) <- c("ID", "REF", "ALT")
vcf_2320_join <- left_join(vcf_2320, mafs_filters_tojoin, by="ID")
vcf_2320_join <- vcf_2320_join[,-4]
vcf_2320_join <- vcf_2320_join[,-4]
vcf_2320_join <- vcf_2320_join[,c(1:3, 7:8, 4:6)]
colnames(vcf_2320_join) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

#Export to text file
write.table(vcf_2320_join, "2320snps.txt", quote=FALSE, sep="\t", row.names=FALSE)
#Will change the ending from .txt to .vcf and will add the vcf header information in text editor

##Trying to get the allele frequency trends for our candidate SNPs

#Flip the allele frequency dataframe
indfreq_filters_flipped <- indfreq_filters
indfreq_filters_flipped <- t(indfreq_filters_flipped)
colnames(indfreq_filters_flipped) <- seq(1:400)

indfreq_qvalues <- cbind(TAB_manhattan_ocean, indfreq_filters_flipped)
indfreq_qvalues_cand <- subset(indfreq_qvalues, qvalues<0.05)

indfreq_cand <- indfreq_qvalues_cand
indfreq_cand <- indfreq_cand[,-1]
indfreq_cand <- indfreq_cand[,-2:-4]
rownames(indfreq_cand) <- indfreq_cand[,1]
indfreq_cand <- indfreq_cand[,-1]
indfreq_cand_flip <- t(indfreq_cand)
indfreq_cand_flip <- as.data.frame(indfreq_cand_flip)

#Order by latitude
indfreq_cand_flip$latitude <- env_ocean_filters$lat
indfreq_cand_flip <- indfreq_cand_flip %>% relocate(latitude, .before = ptg000092l_596944)
indfreq_cand_flip_order <- indfreq_cand_flip[order(indfreq_cand_flip$latitude, decreasing = TRUE),]

#Subset to the 14 candidates on genes
indfreq_cand_flip_order_gene <- as.data.frame(indfreq_cand_flip_order$latitude)
indfreq_cand_flip_order_gene$ptg000011l_21814 <- indfreq_cand_flip_order$ptg000011l_21814
indfreq_cand_flip_order_gene$ptg000092l_1178693 <- indfreq_cand_flip_order$ptg000092l_1178693

indfreq_cand_flip_order_gene$ptg000092l_1178706 <- indfreq_cand_flip_order$ptg000092l_1178706

indfreq_cand_flip_order_gene$ptg000318l_265757 <- indfreq_cand_flip_order$ptg000318l_265757

indfreq_cand_flip_order_gene$ptg000008l_15485912 <- indfreq_cand_flip_order$ptg000008l_15485912

indfreq_cand_flip_order_gene$ptg000029l_20218499 <- indfreq_cand_flip_order$ptg000029l_20218499

indfreq_cand_flip_order_gene$ptg000008l_15485909 <- indfreq_cand_flip_order$ptg000008l_15485909

indfreq_cand_flip_order_gene$ptg000013l_15215686 <- indfreq_cand_flip_order$ptg000013l_15215686

indfreq_cand_flip_order_gene$ptg000001l_3598709 <- indfreq_cand_flip_order$ptg000001l_3598709

indfreq_cand_flip_order_gene$ptg000038l_1789331 <- indfreq_cand_flip_order$ptg000038l_1789331

indfreq_cand_flip_order_gene$ptg000038l_14775337 <- indfreq_cand_flip_order$ptg000038l_14775337

indfreq_cand_flip_order_gene$ptg000065l_12653405 <- indfreq_cand_flip_order$ptg000065l_12653405

indfreq_cand_flip_order_gene$ptg000074l_2906217 <- indfreq_cand_flip_order$ptg000074l_2906217

indfreq_cand_flip_order_gene$ptg000038l_1789904 <- indfreq_cand_flip_order$ptg000038l_1789904

#Create vector of differences in allele frequencies

freq_range <- as.vector(c(0.0178, 0.0141, 0.13022, 0.3763, 0.1253, 0.028, 0.1329, 0.0884,
                          0.0242, 0.0859, 0, 0.0546, 0.1374, 0.0176))
library(plotrix)
mean(freq_range) #0.08805143
std.error(freq_range) #0.02595717

##Import depth information

seqdepth <- read_table("bam_list_realigned_mindp132_maxdp4000_minind0.50_minq20.pos.gz")

mean(seqdepth$totDepth)

