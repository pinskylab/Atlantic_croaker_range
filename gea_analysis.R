## GEA Analysis Using LFMM and RDA ##

##Following worked example by Brenna Forester https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html#WE_11

#Import the genetic data
#Genotypes called in ANGSD and saved to PLINK format (.tped) then converted to .raw (-recodeA in plink) in PLINK

#Read in .raw data using adegenet

library(adegenet) #version 2.1.10

genotypes <- read.PLINK("genotypes_mindepth1_a.raw")

genotypes_matrix <- as.matrix(genotypes)

genotypes_mindepth132 <- read.PLINK("genotypes_mindepth132_a.raw")

genotypes_132_matrix <- as.matrix(genotypes_mindepth132)

##Missing genotype imputation using the LEA package

library(LEA)


##Using population-level allele frequencies instead

#Load in .mafs allele frequency files from ANGSD (one file per population)

de_freq <- read.table("DE_global_snp_list_bam_list_realigned_mindp132_maxdp4000_minind0_snppopminind1.mafs.txt", header = T)
fl_freq <- read.table("FL_global_snp_list_bam_list_realigned_mindp132_maxdp4000_minind0_snppopminind1.mafs.txt", header = T)
ga_freq <- read.table("GA_global_snp_list_bam_list_realigned_mindp132_maxdp4000_minind0_snppopminind1.mafs.txt", header = T)
lb_freq <- read.table("LB_global_snp_list_bam_list_realigned_mindp132_maxdp4000_minind0_snppopminind1.mafs.txt", header = T)
md_freq <- read.table("MD_global_snp_list_bam_list_realigned_mindp132_maxdp4000_minind0_snppopminind1.mafs.txt", header = T)
nc_freq <- read.table("NC_global_snp_list_bam_list_realigned_mindp132_maxdp4000_minind0_snppopminind1.mafs.txt", header = T)
ob_freq <- read.table("OB_global_snp_list_bam_list_realigned_mindp132_maxdp4000_minind0_snppopminind1.mafs.txt", header = T)
sc_freq <- read.table("SC_global_snp_list_bam_list_realigned_mindp132_maxdp4000_minind0_snppopminind1.mafs.txt", header = T)
va_freq <- read.table("VA_global_snp_list_bam_list_realigned_mindp132_maxdp4000_minind0_snppopminind1.mafs.txt", header = T)

#Combine chromo and position columns

de_freq$chromo_position <- paste(de_freq$chromo, de_freq$position, sep="_")
fl_freq$chromo_position <- paste(fl_freq$chromo, fl_freq$position, sep="_")
ga_freq$chromo_position <- paste(ga_freq$chromo, ga_freq$position, sep="_")
lb_freq$chromo_position <- paste(lb_freq$chromo, lb_freq$position, sep="_")
md_freq$chromo_position <- paste(md_freq$chromo, md_freq$position, sep="_")
nc_freq$chromo_position <- paste(nc_freq$chromo, nc_freq$position, sep="_")
ob_freq$chromo_position <- paste(ob_freq$chromo, ob_freq$position, sep="_")
sc_freq$chromo_position <- paste(sc_freq$chromo, sc_freq$position, sep="_")
va_freq$chromo_position <- paste(va_freq$chromo, va_freq$position, sep="_")

#Subset data frame to chromo_position column and frequency (knownEM)

de_freq_subset <- subset(de_freq, select = c(chromo_position, knownEM))
fl_freq_subset <- subset(fl_freq, select = c(chromo_position, knownEM))
ga_freq_subset <- subset(ga_freq, select = c(chromo_position, knownEM))
lb_freq_subset <- subset(lb_freq, select = c(chromo_position, knownEM))
md_freq_subset <- subset(md_freq, select = c(chromo_position, knownEM))
nc_freq_subset <- subset(nc_freq, select = c(chromo_position, knownEM))
ob_freq_subset <- subset(ob_freq, select = c(chromo_position, knownEM))
sc_freq_subset <- subset(sc_freq, select = c(chromo_position, knownEM))
va_freq_subset <- subset(va_freq, select = c(chromo_position, knownEM))

#Rename knownEM column to specify the population it came from

colnames(de_freq_subset) <- c("chromo_position", "de_knownEM")
colnames(fl_freq_subset) <- c("chromo_position", "fl_knownEM")
colnames(ga_freq_subset) <- c("chromo_position", "ga_knownEM")
colnames(lb_freq_subset) <- c("chromo_position", "lb_knownEM")
colnames(md_freq_subset) <- c("chromo_position", "md_knownEM")
colnames(nc_freq_subset) <- c("chromo_position", "nc_knownEM")
colnames(ob_freq_subset) <- c("chromo_position", "ob_knownEM")
colnames(sc_freq_subset) <- c("chromo_position", "sc_knownEM")
colnames(va_freq_subset) <- c("chromo_position", "va_knownEM")

#Join population dataframes

freq_list <- list(de_freq_subset, fl_freq_subset, ga_freq_subset, lb_freq_subset, md_freq_subset, nc_freq_subset, ob_freq_subset, sc_freq_subset, va_freq_subset)

freq_combined <- Reduce(function(x, y) merge(x, y, all=TRUE), freq_list) #all SNPs but there are NAs for some populations (380,773 SNPs)

freq_combined_full <- na.omit(freq_combined) #SNPs that have an allele frequency for each population (35,305 SNPs)

row.names(freq_combined_full) <- as.character(freq_combined_full$chromo_position)

#Flip rows and columns to match Capblancq and Forester RDA example code

freq_combined_full <- t(freq_combined_full)

freq_combined_full <- freq_combined_full[-1,] #Looks good

freq_combined_full <- as.data.frame(freq_combined_full)

freq_combined_full <- sapply(freq_combined_full, as.numeric)

##Load in environmental data

env_pop <- read.csv("Env_Data_byPop.csv", header = T)
#interested in latitude, depth, temperature, and salinity

#Convert degrees F to degrees C

env_pop$noaa_oct_temp <- (env_pop$noaa_oct_temp_f - 32)*(5/9)

row.names(env_pop) <- as.character(env_pop$population)

#Subset to lat, depth, bottom_temp, and bottom_salinity for analysis
env_pop_pred <- subset(env_pop, select = c(avg_lat, avg_depth, avg_temp, avg_salinity))

##RDA

library(vegan) #version 2.6-2

rda_pop <- rda(freq_combined_full ~ ., data=env_pop_pred, scale=T)
rda_pop
#                Inertia Proportion Rank
#Total         3.528e+04  1.000e+00
#Constrained   1.739e+04  4.929e-01    4
#Unconstrained 1.789e+04  5.071e-01    4
#Inertia is correlations

#Calculate adjusted R-squared
RsquareAdj(rda_pop)
#r.squared=0.492937 adj.r.squared=-0.0141259

summary(eigenvals(rda_pop, model = "constrained"))
#Most variance explained by RDA1

screeplot(rda_pop)

plot(rda_pop)

#Test for significance
signif.full <- anova.cca(rda_pop, parallel=getOption("mc.cores"))

signif.axis <- anova.cca(rda_pop, by="axis", parallel=getOption("mc.cores"))

##Identify candidate SNPs involved in local adaptation

load.rda.pop <- scores(rda_pop, choices=c(1:4), display="species")

hist(load.rda.pop[,1], main="Loadings on RDA1")
hist(load.rda.pop[,2], main="Loadings on RDA2")
hist(load.rda.pop[,3], main="Loadings on RDA3")
hist(load.rda.pop[,4], main="Loadings on RDA4")

cand1_pop <- outliers(load.rda.pop[,1],2.5) #15
cand2_pop <- outliers(load.rda.pop[,2],2.5) #25
cand3_pop <- outliers(load.rda.pop[,3],2.5) #65
cand4_pop <- outliers(load.rda.pop[,4],2.5) #244

ncand_pop <- length(cand1_pop) + length(cand2_pop) + length(cand3_pop) + length(cand4_pop)
ncand_pop #349

## Function rdadapt
source("./rdadapt.R")

## Running the function with K = 2
rdadapt_pop <-rdadapt(rda_pop, 2)

#One critical step when conducting a genome scan is to set a pertinent p-value threshold to
#identify the outlier loci. Here, we used a Bonferroni correction to account for multiple testing.

#**Note:** the rdadapt function returns both p-values and q-values, which means it is possible
#to use a FDR (False Discovery Rate) approach instead of a p-value threshold to identify outliers;
#see [Francois et al. 2015](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.13513).

## P-values threshold after Bonferroni correction
thres_pop <- 0.01/length(rdadapt_pop$p.values) #2.83e-7

## Identifying the loci that are below the p-value threshold
outliers_pop <- data.frame(Loci = colnames(freq_combined_full)[which(rdadapt_pop$p.values<0.05)], p.value =
                            rdadapt_pop$p.values[which(rdadapt_pop$p.values<0.05)], contig =
                            unlist(lapply(strsplit(colnames(freq_combined_full)[which(rdadapt_pop$p.values<0.05)],
                                                   split = "_"), function(x) x[1])))
#yields 6 possible SNPs

#To avoid redundancy among loci that are in strong physical linkage along the genome
#we only kept the outlier locus with the lowest p-values for each genomic contig.

## Top hit outlier per contig
outliers_pop <- outliers_pop[order(outliers_pop$contig, outliers_pop$p.value),]

## List of outlier names
outliers_rdadapt_pop <- as.character(outliers_pop$Loci[!duplicated(outliers_pop$contig)])
#5 contigs

#Once the outliers have been identified, it can be useful to visualize their distribution in
#comparison with neutral loci using either an RDA biplot or a Manhattan plot.

## Formatting table for ggplot
locus_scores_pop <- scores(rda_pop, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci_pop <- data.frame(names = row.names(locus_scores_pop), locus_scores_pop)
TAB_loci_pop$type <- "Neutral"
TAB_loci_pop$type[TAB_loci_pop$names%in%outliers_pop$Loci] <- "All outliers"
TAB_loci_pop$type[TAB_loci_pop$names%in%outliers_rdadapt_pop] <- "Top outliers"
TAB_loci_pop$type <- factor(TAB_loci_pop$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci_pop <- TAB_loci_pop[order(TAB_loci_pop$type),]
TAB_var_pop <- as.data.frame(scores(rda_pop, choices=c(1,2), display="bp")) # pull the biplot scores

## Biplot of RDA loci and variables scores
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci_pop, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var_pop, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var_pop, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var_pop)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

## Manhattan plot
Outliers_pop <- rep("Neutral", length(colnames(freq_combined_full)))
Outliers_pop[colnames(freq_combined_full)%in%outliers_pop$Loci] <- "All outliers"
Outliers_pop[colnames(freq_combined_full)%in%outliers_rdadapt_pop] <- "Top outliers"
Outliers_pop <- factor(Outliers_pop, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_manhatan_pop <- data.frame(pos = 1:length(colnames(freq_combined_full)),
                           pvalues = rdadapt_pop$p.values,
                           Outliers = Outliers_pop)
TAB_manhatan_pop <- TAB_manhatan_pop[order(TAB_manhatan_pop$Outliers),]
ggplot(data = TAB_manhatan_pop) +
  geom_point(aes(x=pos, y=-log10(pvalues), col = Outliers), size=1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  xlab("Loci") + ylab("-log10(p.values)") +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  facet_wrap(~"Manhattan plot", nrow = 3) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(),
        legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(),
        legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))


##Population-level RDA using October mean temperature

env_pop_oct <- subset(env_pop, select = c(avg_lat, avg_depth, noaa_oct_temp, avg_salinity))

rda_pop_oct <- rda(freq_combined_full ~ ., data=env_pop_oct, scale=T)
rda_pop_oct
#                Inertia Proportion Rank
#Total         3.528e+04  1.000e+00
#Constrained   1.726e+04  4.892e-01    4
#Unconstrained 1.802e+04  5.108e-01    4
#Inertia is correlations

#Calculate adjusted R-squared
RsquareAdj(rda_pop_oct)
#r.squared=0.4892302 adj.r.squared=-0.02153966

summary(eigenvals(rda_pop_oct, model = "constrained"))

screeplot(rda_pop_oct)

plot(rda_pop_oct)

#Test for significance
signif.full <- anova.cca(rda_pop_oct, parallel=getOption("mc.cores"))

signif.axis <- anova.cca(rda_pop_oct, by="axis", parallel=getOption("mc.cores"))

##Identify candidate SNPs involved in local adaptation

load.rda.pop.oct <- scores(rda_pop_oct, choices=c(1:4), display="species")

hist(load.rda.pop.oct[,1], main="Loadings on RDA1")
hist(load.rda.pop.oct[,2], main="Loadings on RDA2")
hist(load.rda.pop.oct[,3], main="Loadings on RDA3")
hist(load.rda.pop.oct[,4], main="Loadings on RDA4")

cand1_pop_oct <- outliers(load.rda.pop.oct[,1],2.5) #5
cand2_pop_oct <- outliers(load.rda.pop.oct[,2],2.5) #33
cand3_pop_oct <- outliers(load.rda.pop.oct[,3],2.5) #73
cand4_pop_oct <- outliers(load.rda.pop.oct[,4],2.5) #348

ncand_pop_oct <- length(cand1_pop_oct) + length(cand2_pop_oct) + length(cand3_pop_oct) + length(cand4_pop_oct)
ncand_pop_oct #459

##Run without the scale=T in the RDA call

rda_pop_oct_not <- rda(freq_combined_full ~ ., data=env_pop_oct)
rda_pop_oct_not

#Calculate adjusted R-squared
RsquareAdj(rda_pop_oct_not)
#r.squared=0.4920667 adj.r.squared=-0.01586664

summary(eigenvals(rda_pop_oct_not, model = "constrained"))

screeplot(rda_pop_oct_not)

plot(rda_pop_oct_not)

##Looks about the same as

##Trying out individual allele frequencies from PCAngsd as the genetic input
#individual frequencies with SNPs found at 200+ individuals (50%+)

library(reticulate)
np <- import("numpy")

indfreq <- np$load("pcangsd_freq_filters_snp200.indf.npy")
#400 individuals, 121,156 SNPs

indfreq <- as.data.frame(indfreq)
indfreq[1:5,1:5] #each individual has a row and each SNP has a column
#Need to transfer over the contig/position numbers for each SNP
sum(is.na(indfreq)) #No NAs, so no imputation needed

#Load in mafs file from ANGSD to get contig/position numbers for each SNP
mafs <- read.table("freq_filters_snp200.mafs.txt", header = T)

#Combine chromo and position columns
mafs$chromo_position <- paste(mafs$chromo, mafs$position, sep="_")
snp_position <- as.vector(mafs$chromo_position)

colnames(indfreq) <- snp_position
indfreq[1:5, 1:5] #Now each column has the contig and position number as the header

#Check that the entire data frame looks okay by looking at a subset
indfreq_subset <- indfreq[-c(5001:121156)]
indfreq_subset #Looks good, each individual has a row, and each SNP has a column with the right header

##Load in environmental data

env <- read.csv("Env_Data.csv", header = T)
#interested in latitude, depth, temperature, and salinity

#Subset to lat, depth, bottom_temp, and bottom_salinity for analysis
env_pred <- subset(env, select = c(lat, depth, bottom_temp, bottom_salinity))

##RDA

library(vegan) #version 2.6-4

first_rda <- rda(indfreq ~ ., data=env_pred, scale=T)
first_rda
#                Inertia Proportion Rank
#Total         1.212e+05  1.000e+00
#Constrained   4.402e+03  3.633e-02    4
#Unconstrained 1.168e+05  9.637e-01   18
#Inertia is correlations

#Calculate adjusted R-squared
RsquareAdj(first_rda)
#r.squared=0.03633412 adj.r.squared=0.02657548

summary(eigenvals(first_rda, model = "constrained"))
#Most variance explained by RDA1

screeplot(first_rda)

plot(first_rda)

#Test for significance
signif.full <- anova.cca(first_rda, parallel=getOption("mc.cores"))

signif.axis <- anova.cca(first_rda, by="axis", parallel=getOption("mc.cores"))

##Identify candidate SNPs involved in local adaptation

load.rda <- scores(first_rda, choices=c(1:4), display="species")

hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3")
hist(load.rda[,4], main="Loadings on RDA4")

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3) #0, but seems like all SNPs loading on the edges
cand2 <- outliers(load.rda[,2],3) #3128
cand3 <- outliers(load.rda[,3],3) #1123
cand4 <- outliers(load.rda[,4],3) #1510

ncand <- length(cand1) + length(cand2) + length(cand3) + length(cand4)
ncand #5761

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(4,times=length(cand3)), names(cand4), unname(cand4))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

#Alternate method for identifying candidate SNPs: procedure described in Capblancq et al. (2018) where
#outliers are identified based on their extremeness along a distribution of Mahalanobis distances estimated between each locus
#and the center of the RDA space using a certain number of axes (K).

#Load the [rdadapt] function, described in Capblancq et al. (2018) and use it to conduct the genome scan.

library(robust)
library(qvalue)

## Function rdadapt
source("./rdadapt.R")

## Running the function with K = 2
rdadapt_env<-rdadapt(first_rda, 2)

#One critical step when conducting a genome scan is to set a pertinent p-value threshold to
#identify the outlier loci. Here, we used a Bonferroni correction to account for multiple testing.

#**Note:** the rdadapt function returns both p-values and q-values, which means it is possible
#to use a FDR (False Discovery Rate) approach instead of a p-value threshold to identify outliers;
#see [Francois et al. 2015](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.13513).

## P-values threshold after Bonferroni correction
thres_env <- 0.01/length(rdadapt_env$p.values) #8.23e-8

## Identifying the loci that are below the p-value threshold
outliers_k2 <- data.frame(Loci = colnames(indfreq)[which(rdadapt_env$p.values<thres_env)], p.value =
                         rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig =
                         unlist(lapply(strsplit(colnames(indfreq)[which(rdadapt_env$p.values<thres_env)],
                                                split = "_"), function(x) x[1])))
#yields 59,799 possible SNPs

#To avoid redundancy among loci that are in strong physical linkage along the genome
#we only kept the outlier locus with the lowest p-values for each genomic contig.

## Top hit outlier per contig
outliers_k2 <- outliers_k2[order(outliers_k2$contig, outliers_k2$p.value),]

## List of outlier names
outliers_rdadapt_env <- as.character(outliers_k2$Loci[!duplicated(outliers_k2$contig)])
#186 contigs

#Once the outliers have been identified, it can be useful to visualize their distribution in
#comparison with neutral loci using either an RDA biplot or a Manhattan plot.

## Formatting table for ggplot
locus_scores <- scores(first_rda, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers_k2$Loci] <- "All outliers"
TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(first_rda, choices=c(1,2), display="bp")) # pull the biplot scores

## Biplot of RDA loci and variables scores
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

## Manhattan plot
Outliers <- rep("Neutral", length(colnames(indfreq)))
Outliers[colnames(indfreq)%in%outliers_k2$Loci] <- "All outliers"
Outliers[colnames(indfreq)%in%outliers_rdadapt_env] <- "Top outliers"
Outliers <- factor(Outliers, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_manhatan <- data.frame(pos = 1:length(colnames(indfreq)),
                           pvalues = rdadapt_env$p.values,
                           Outliers = Outliers)
TAB_manhatan <- TAB_manhatan[order(TAB_manhatan$Outliers),]
ggplot(data = TAB_manhatan) +
  geom_point(aes(x=pos, y=-log10(pvalues), col = Outliers), size=1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  xlab("Loci") + ylab("-log10(p.values)") +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  facet_wrap(~"Manhattan plot", nrow = 3) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(),
        legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(),
        legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

##Using approach that accounts for population structure and runs a pRDA

## Running a PCA on neutral genetic markers
pca_200 <- rda(indfreq, scale=T) # PCA in vegan uses the rda() call without any predictors

#Screeplot of the PCA eigenvalues:

screeplot(pca_200, type = "barplot", npcs=10, main="PCA Eigenvalues")

#Based on the screeplot, one PC would be a reasonable set to retain as a proxy for
#neutral population structure in downstream analyses.

## Neutral population structure table
PCs <- scores(pca_200, choices=1, display="sites", scaling=0)
PopStruct <- data.frame(Individual = seq(1:400), PCs)
colnames(PopStruct) <- c("Individual", "PC1")
PopStruct[1:5,]

#Add PC1 to the environmental variables data frame

env_pred$pc1 <- PopStruct$PC1

#pRDA accounting for population structure using PC1
#Using the same 4 environmental variables as explanatory variables and PC1 as a conditioning variable
#to account for neutral population structure.

RDA_pc1 <- rda(indfreq ~ lat + depth + bottom_temp + bottom_salinity + Condition(pc1),  data=env_pred, scale=T)
RDA_pc1
#Inertia Proportion Rank
#Total         1.212e+05  1.000e+00
#Conditional   1.211e+05  9.994e-01    1
#Constrained   9.528e-01  7.864e-06    4
#Unconstrained 7.602e+01  6.275e-04   48
#Inertia is correlations

#Calculate adjusted R-squared
RsquareAdj(RDA_pc1)
#r.squared=7.863891e-06 adj.r.squared=1.49754e-06

summary(eigenvals(RDA_pc1, model = "constrained"))
#Most variance explained by RDA1

screeplot(RDA_pc1)

plot(RDA_pc1)

load.rda.pc1 <- scores(RDA_pc1, choices=c(1:3), display="species")

hist(load.rda.pc1[,1], main="Loadings on RDA1")
hist(load.rda.pc1[,2], main="Loadings on RDA2")
hist(load.rda.pc1[,3], main="Loadings on RDA3")

cand1_pc1 <- outliers(load.rda.pc1[,1],3) #2670
cand2_pc1 <- outliers(load.rda.pc1[,2],3) #3726
cand3_pc1 <- outliers(load.rda.pc1[,3],3) #417

ncand_pc1 <- length(cand1_pc1) + length(cand2_pc1) + length(cand3_pc1)
ncand_pc1 #6813

#Organize results by making one data frame with the axis, SNP name, loading, and correlation with each predictor

cand1_pc1 <- cbind.data.frame(rep(1,times=length(cand1_pc1)), names(cand1_pc1), unname(cand1_pc1))
cand2_pc1 <- cbind.data.frame(rep(2,times=length(cand2_pc1)), names(cand2_pc1), unname(cand2_pc1))
cand3_pc1 <- cbind.data.frame(rep(3,times=length(cand3_pc1)), names(cand3_pc1), unname(cand3_pc1))

colnames(cand1_pc1) <- colnames(cand2_pc1) <- colnames(cand3_pc1) <- c("axis","snp","loading")

cand_pc1 <- rbind(cand1_pc1, cand2_pc1, cand3_pc1)
cand_pc1$snp <- as.character(cand_pc1$snp)

#Add in correlations of each candidate SNP with the four environmental predictors:

foo <- matrix(nrow=(ncand_pc1), ncol=4)  # 4 columns for 4 predictors
colnames(foo) <- c("lat","depth","bottom_temp","bottom_salinity")

env_pred_cand <- subset(env_pred, select = c(lat, depth, bottom_temp, bottom_salinity))

for (i in 1:length(cand_pc1$snp)) {
  nam <- cand_pc1[i,2]
  snp.gen <- indfreq[,nam]
  foo[i,] <- apply(env_pred_cand,2,function(x) cor(x,snp.gen))
}

cand_pc1 <- cbind.data.frame(cand_pc1,foo)
head(cand_pc1)

#Investigate the SNP candidates
#Look for duplicate detections (SNPs identified as candidates on more than one RDA axis)

length(cand_pc1$snp[duplicated(cand_pc1$snp)]) #1229 duplicate detections

foo <- cbind(cand_pc1$axis, duplicated(cand_pc1$snp))
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) #812 duplicates on axis 2
table(foo[foo[,1]==3,2]) #no duplicates on axis 3

cand_pc1 <- cand_pc1[!duplicated(cand_pc1$snp),] # remove duplicate detections
#now have 5,584 unique candidate SNPs

#See which of the predictors each SNP is most strongly correlated with
for (i in 1:length(cand_pc1$snp)) {
  bar <- cand_pc1[i,]
  cand_pc1[i,8] <- names(which.max(abs(bar[4:7]))) # gives the variable
  cand_pc1[i,9] <- max(abs(bar[4:7]))              # gives the correlation
}

colnames(cand_pc1)[8] <- "predictor"
colnames(cand_pc1)[9] <- "correlation"

table(cand_pc1$predictor) #bottom_salinity = 5476, depth = 108
#temperature seems to be close behind salinity as the second highest predictor for individual SNPs

#Plot the SNPs

sel_pc1 <- cand_pc1$snp
env_pc1 <- cand_pc1$predictor
env_pc1[env_pc1=="lat"] <- '#1f78b4'
env_pc1[env_pc1=="depth"] <- '#a6cee3'
env_pc1[env_pc1=="bottom_temp"] <- '#6a3d9a'
env_pc1[env_pc1=="bottom_salinity"] <- '#e31a1c'

# color by predictor:
col.pred <- rownames(RDA_pc1$CCA$v) # pull the SNP names

for (i in 1:length(sel_pc1)) {           # color code candidate SNPs, something going wrong here
  foo <- match(sel_pc1[i],col.pred)
  col.pred[foo] <- env_pc1[i]
}

col.pred[grep("chr",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c')

# axes 1 & 2
plot(RDA_pc1, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(RDA_pc1, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(RDA_pc1, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(RDA_pc1, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("latitude","depth","bottom temperature","bottom salinity"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

#pRDA using the October mean temperature

#Add PC1 to the environmental variables data frame

env_oct_pred$pc1 <- PopStruct$PC1

#pRDA accounting for population structure using PC1
#Using the same 4 environmental variables as explanatory variables and PC1 as a conditioning variable
#to account for neutral population structure.

rda_oct_pc1 <- rda(indfreq ~ lat + depth + noaa_oct_temp + bottom_salinity + Condition(pc1),  data=env_oct_pred, scale=T)
rda_oct_pc1
#                Inertia Proportion Rank
#Total         1.212e+05  1.000e+00
#Conditional   1.211e+05  9.994e-01    1
#Constrained   7.089e-01  5.851e-06    4
#Unconstrained 7.626e+01  6.295e-04   48
#Inertia is correlations

#Calculate adjusted R-squared
RsquareAdj(rda_oct_pc1)
#r.squared=5.851264e-06 adj.r.squared=-5.406292e-07

summary(eigenvals(rda_oct_pc1, model = "constrained"))
#Most variance explained by RDA1, but some explained by RDA2 and RDA3
#Proportion explained RDA1: 0.86 RDA2: 0.11 RDA3: 0.02

screeplot(rda_oct_pc1)

plot(rda_oct_pc1)

load.rda.oct.pc1 <- scores(rda_oct_pc1, choices=c(1:3), display="species")

hist(load.rda.oct.pc1[,1], main="Loadings on RDA1")
hist(load.rda.oct.pc1[,2], main="Loadings on RDA2")
hist(load.rda.oct.pc1[,3], main="Loadings on RDA3")

cand1_oct_pc1 <- outliers(load.rda.oct.pc1[,1],3) #2802
cand2_oct_pc1 <- outliers(load.rda.oct.pc1[,2],3) #3797
cand3_oct_pc1 <- outliers(load.rda.oct.pc1[,3],3) #284

ncand_oct_pc1 <- length(cand1_oct_pc1) + length(cand2_oct_pc1) + length(cand3_oct_pc1)
ncand_oct_pc1 #6883

#Organize results by making one data frame with the axis, SNP name, loading, and correlation with each predictor

cand1_oct_pc1 <- cbind.data.frame(rep(1,times=length(cand1_oct_pc1)), names(cand1_oct_pc1), unname(cand1_oct_pc1))
cand2_oct_pc1 <- cbind.data.frame(rep(2,times=length(cand2_oct_pc1)), names(cand2_oct_pc1), unname(cand2_oct_pc1))
cand3_oct_pc1 <- cbind.data.frame(rep(3,times=length(cand3_oct_pc1)), names(cand3_oct_pc1), unname(cand3_oct_pc1))

colnames(cand1_oct_pc1) <- colnames(cand2_oct_pc1) <- colnames(cand3_oct_pc1) <- c("axis","snp","loading")

cand_oct_pc1 <- rbind(cand1_oct_pc1, cand2_oct_pc1, cand3_oct_pc1)
cand_oct_pc1$snp <- as.character(cand_oct_pc1$snp)

#Add in correlations of each candidate SNP with the four environmental predictors:

foo <- matrix(nrow=(ncand_oct_pc1), ncol=4)  # 4 columns for 4 predictors
colnames(foo) <- c("lat","depth","noaa_oct_temp","bottom_salinity")

env_oct_pred_cand <- subset(env_oct_pred, select = c(lat, depth, noaa_oct_temp, bottom_salinity))

for (i in 1:length(cand_oct_pc1$snp)) {
  nam <- cand_oct_pc1[i,2]
  snp.gen <- indfreq[,nam]
  foo[i,] <- apply(env_oct_pred_cand,2,function(x) cor(x,snp.gen))
}

cand_oct_pc1 <- cbind.data.frame(cand_oct_pc1,foo)
head(cand_oct_pc1)

#Investigate the SNP candidates
#Look for duplicate detections (SNPs identified as candidates on more than one RDA axis)

length(cand_oct_pc1$snp[duplicated(cand_oct_pc1$snp)]) #1177 duplicate detections

foo <- cbind(cand_oct_pc1$axis, duplicated(cand_oct_pc1$snp))
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) #998 duplicates on axis 2
table(foo[foo[,1]==3,2]) #179 duplicates on axis 3

cand_oct_pc1 <- cand_oct_pc1[!duplicated(cand_oct_pc1$snp),] # remove duplicate detections
#now have 5,706 unique candidate SNPs

#See which of the predictors each SNP is most strongly correlated with
for (i in 1:length(cand_oct_pc1$snp)) {
  bar <- cand_oct_pc1[i,]
  cand_oct_pc1[i,8] <- names(which.max(abs(bar[4:7]))) # gives the variable
  cand_oct_pc1[i,9] <- max(abs(bar[4:7]))              # gives the correlation
}

colnames(cand_oct_pc1)[8] <- "predictor"
colnames(cand_oct_pc1)[9] <- "correlation"

table(cand_oct_pc1$predictor) #bottom_salinity = 4383, depth = 22, noaa_oct_temp = 1301
#temperature seems to be close behind salinity as the second highest predictor for individual SNPs

#Identify how many unique contigs the SNPs are on
#Separate the snp column into contig and position
library(dplyr)
library(tidyr)
cand_oct_pc1_contig <- cand_oct_pc1 %>% separate_wider_delim(snp, "_", names=c('contig',  'position'))

cand_oct_pc1$contig <- cand_oct_pc1_contig$contig
cand_oct_pc1$position <- cand_oct_pc1_contig$position

unique_oct_contigs <- unique(cand_oct_pc1[c("contig")]) #146 unique contigs

##Plotting SNPs

#Formatting table for ggplot
locus_scores_oct <- scores(rda_oct_pc1, choices=c(1:2), display="species", scaling="none")
TAB_loci_oct <- data.frame(names = row.names(locus_scores_oct), locus_scores_oct)
TAB_loci_oct$type <- "Neutral"
TAB_loci_oct$type[TAB_loci_oct$names%in%cand_oct_pc1$snp] <- "Outlier"
TAB_loci_oct$type <- factor(TAB_loci_oct$type, levels = c("Neutral", "Outlier"))
TAB_loci_oct <- TAB_loci_oct[order(TAB_loci_oct$type),]
TAB_var_oct <- as.data.frame(scores(rda_oct_pc1, choices=c(1,2), display="bp")) # pull the biplot scores

## Biplot of RDA loci and variables scores
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci_oct, aes(x=RDA1*20, y=RDA2*20, colour = type), shape=19, size = 2) +
  scale_color_manual(values = c("gray60", "#F9A242FF")) +
  geom_segment(data = TAB_var_oct, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.2, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var_oct, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var_oct)), size = 3.75, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  #facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(axis.text.x= element_text(size=11), axis.text.y=element_text(size=12)) +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=11))

#Plot the SNPs by predictor
#Subset based on predictor
temp_oct_pc1 <- subset(cand_oct_pc1, predictor=='noaa_oct_temp', select=c(snp, predictor))
salinity_oct_pc1 <- subset(cand_oct_pc1, predictor=='bottom_salinity', select=c(snp, predictor))
depth_oct_pc1 <- subset(cand_oct_pc1, predictor=='depth', select=c(snp, predictor))

TAB_loci_oct$predictor <- "Neutral"
TAB_loci_oct$predictor[TAB_loci_oct$names%in%temp_oct_pc1$snp] <- "Temperature"
TAB_loci_oct$predictor[TAB_loci_oct$names%in%salinity_oct_pc1$snp] <- "Salinity"
TAB_loci_oct$predictor[TAB_loci_oct$names%in%depth_oct_pc1$snp] <- "Depth"
TAB_loci_oct$predictor <- factor(TAB_loci_oct$predictor, levels = c("Neutral", "Temperature", "Salinity", "Depth"))
TAB_loci_oct <- TAB_loci_oct[order(TAB_loci_oct$predictor),]

ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci_oct, aes(x=RDA1*20, y=RDA2*20, colour = predictor), shape=19, size = 2) +
  scale_color_manual(values = c("gray60", '#e31a1c', '#1f78b4','#6a3d9a')) +
  geom_segment(data = TAB_var_oct, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.2, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var_oct, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var_oct)), size = 3.75, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  #facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(axis.text.x= element_text(size=11), axis.text.y=element_text(size=12)) +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=11))

#Base r plot version
col.pred <- TAB_loci_oct

col.pred$color[col.pred$predictor=="Neutral"] <- "gray60"
col.pred$color[col.pred$predictor=="Temperature"] <- '#e31a1c'
col.pred$color[col.pred$predictor=="Salinity"] <- '#1f78b4'
col.pred$color[col.pred$predictor=="Depth"] <- '#6a3d9a'

bg <- c("gray60",'#e31a1c','#1f78b4','#6a3d9a')

col.pred[grep("chr",col.pred)] <- '#f1eef6' # non-candidate SNPs #not using this code chunk for now
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c')

# axes 1 & 2
plot(rda_oct_pc1, type="n", scaling=3, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
points(rda_oct_pc1, display="species", pch=21, cex=1, col="gray32", bg=col.pred$color, scaling=3)
#points(rda_oct_pc1, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(rda_oct_pc1, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("Neutral", "Temperature","Salinity","Depth"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

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


#pRDA using the October mean temperature

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
env_oct_filters <- env_oct_pred
env_oct_filters$pc1 <- PopStruct_filters$PC1

#pRDA accounting for population structure using PC1
#Using the same 4 environmental variables as explanatory variables and PC1 as a conditioning variable
#to account for neutral population structure.

rda_oct_filters <- rda(indfreq_filters ~ lat + depth + noaa_oct_temp + bottom_salinity + Condition(pc1),  data=env_oct_filters, scale=T)
rda_oct_filters
#                Inertia Proportion Rank
#Total         2.252e+04  1.000e+00
#Conditional   2.250e+04  9.993e-01    1
#Constrained   2.693e-01  1.196e-05    4
#Unconstrained 1.551e+01  6.887e-04  164

#Calculate adjusted R-squared
RsquareAdj(rda_oct_filters)
#r.squared=1.195716e-05 adj.r.squared=4.978032e-06

summary(eigenvals(rda_oct_filters, model = "constrained"))
#Most variance explained by RDA1 and RDA2, some on RDA3 and RDA4
#Proportion explained RDA1: 0.5306 RDA2: 0.29483 RDA3: 0.10983

screeplot(rda_oct_filters)

plot(rda_oct_filters)

load.rda.oct.filters <- scores(rda_oct_filters, choices=c(1:3), display="species")

hist(load.rda.oct.filters[,1], main="Loadings on RDA1")
hist(load.rda.oct.filters[,2], main="Loadings on RDA2")
hist(load.rda.oct.filters[,3], main="Loadings on RDA3")

cand1_oct_filters <- outliers(load.rda.oct.filters[,1],3) #479
cand2_oct_filters <- outliers(load.rda.oct.filters[,2],3) #288
cand3_oct_filters <- outliers(load.rda.oct.filters[,3],3) #245

ncand_oct_filters <- length(cand1_oct_filters) + length(cand2_oct_filters) + length(cand3_oct_filters)
ncand_oct_filters #1012

#Organize results by making one data frame with the axis, SNP name, loading, and correlation with each predictor

cand1_oct_filters <- cbind.data.frame(rep(1,times=length(cand1_oct_filters)), names(cand1_oct_filters), unname(cand1_oct_filters))
cand2_oct_filters <- cbind.data.frame(rep(2,times=length(cand2_oct_filters)), names(cand2_oct_filters), unname(cand2_oct_filters))
cand3_oct_filters <- cbind.data.frame(rep(3,times=length(cand3_oct_filters)), names(cand3_oct_filters), unname(cand3_oct_filters))

colnames(cand1_oct_filters) <- colnames(cand2_oct_filters) <- colnames(cand3_oct_filters) <- c("axis","snp","loading")

cand_oct_filters <- rbind(cand1_oct_filters, cand2_oct_filters, cand3_oct_filters)
cand_oct_filters$snp <- as.character(cand_oct_filters$snp)

#Add in correlations of each candidate SNP with the four environmental predictors:

foo <- matrix(nrow=(ncand_oct_filters), ncol=4)  # 4 columns for 4 predictors
colnames(foo) <- c("lat","depth","noaa_oct_temp","bottom_salinity")

env_oct_filters_cand <- subset(env_oct_filters, select = c(lat, depth, noaa_oct_temp, bottom_salinity))

for (i in 1:length(cand_oct_filters$snp)) {
  nam <- cand_oct_filters[i,2]
  snp.gen <- indfreq_filters[,nam]
  foo[i,] <- apply(env_oct_filters_cand,2,function(x) cor(x,snp.gen))
}

cand_oct_filters <- cbind.data.frame(cand_oct_filters,foo)
head(cand_oct_filters)

#Investigate the SNP candidates
#Look for duplicate detections (SNPs identified as candidates on more than one RDA axis)

length(cand_oct_filters$snp[duplicated(cand_oct_filters$snp)]) #441 duplicate detections

foo <- cbind(cand_oct_filters$axis, duplicated(cand_oct_filters$snp))
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) #201 duplicates on axis 2
table(foo[foo[,1]==3,2]) #240 duplicates on axis 3

cand_oct_filters <- cand_oct_filters[!duplicated(cand_oct_filters$snp),] # remove duplicate detections
#now have 571 unique candidate SNPs

#See which of the predictors each SNP is most strongly correlated with
for (i in 1:length(cand_oct_filters$snp)) {
  bar <- cand_oct_filters[i,]
  cand_oct_filters[i,8] <- names(which.max(abs(bar[4:7]))) # gives the variable
  cand_oct_filters[i,9] <- max(abs(bar[4:7]))              # gives the correlation
}

colnames(cand_oct_filters)[8] <- "predictor"
colnames(cand_oct_filters)[9] <- "correlation"

table(cand_oct_filters$predictor) #bottom_salinity = 5, depth = 6, noaa_oct_temp = 560
#salinity seems to be close behind temperature as the second highest predictor for individual SNPs

#Identify how many unique contigs the SNPs are on
#Separate the snp column into contig and position
library(dplyr)
library(tidyr)
cand_oct_filters_contig <- cand_oct_filters %>% separate_wider_delim(snp, "_", names=c('contig',  'position'))

cand_oct_filters$contig <- cand_oct_filters_contig$contig
cand_oct_filters$position <- cand_oct_filters_contig$position

unique_oct_filters_contigs <- unique(cand_oct_filters[c("contig")]) #126 unique contigs

##Plotting SNPs

#Formatting table for ggplot
locus_scores_oct_filters <- scores(rda_oct_filters, choices=c(1:2), display="species", scaling="none")
TAB_loci_oct_filters <- data.frame(names = row.names(locus_scores_oct_filters), locus_scores_oct_filters)
TAB_loci_oct_filters$type <- "Neutral"
TAB_loci_oct_filters$type[TAB_loci_oct_filters$names%in%cand_oct_filters$snp] <- "Outlier"
TAB_loci_oct_filters$type <- factor(TAB_loci_oct_filters$type, levels = c("Neutral", "Outlier"))
TAB_loci_oct_filters <- TAB_loci_oct_filters[order(TAB_loci_oct_filters$type),]
TAB_var_oct_filters <- as.data.frame(scores(rda_oct_filters, choices=c(1,2), display="bp")) # pull the biplot scores

## Biplot of RDA loci and variables scores
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci_oct_filters, aes(x=RDA1, y=RDA2, colour = type), shape=19, size = 2) +
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
temp_oct_filters <- subset(cand_oct_filters, predictor=='noaa_oct_temp', select=c(snp, predictor))
salinity_oct_filters <- subset(cand_oct_filters, predictor=='bottom_salinity', select=c(snp, predictor))
depth_oct_filters <- subset(cand_oct_filters, predictor=='depth', select=c(snp, predictor))

TAB_loci_oct_filters$predictor <- "Neutral"
TAB_loci_oct_filters$predictor[TAB_loci_oct_filters$names%in%temp_oct_filters$snp] <- "Temperature"
TAB_loci_oct_filters$predictor[TAB_loci_oct_filters$names%in%salinity_oct_filters$snp] <- "Salinity"
TAB_loci_oct_filters$predictor[TAB_loci_oct_filters$names%in%depth_oct_filters$snp] <- "Depth"
TAB_loci_oct_filters$predictor <- factor(TAB_loci_oct_filters$predictor, levels = c("Neutral", "Temperature", "Salinity", "Depth"))
TAB_loci_oct_filters <- TAB_loci_oct_filters[order(TAB_loci_oct_filters$predictor),]

ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci_oct_filters, aes(x=RDA1, y=RDA2, colour = predictor), shape=19, size = 2) +
  scale_color_manual(values = c("gray60", '#e31a1c', '#1f78b4','#6a3d9a')) +
  #geom_segment(data = TAB_var_oct_filters, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.2, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  #geom_text(data = TAB_var_oct_filters, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var_oct_filters)), size = 3, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  #facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(axis.text.x= element_text(size=11), axis.text.y=element_text(size=12)) +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=11))

#Base r plot version
col.pred <- TAB_loci_oct_filters

col.pred$color[col.pred$predictor=="Neutral"] <- "gray60"
col.pred$color[col.pred$predictor=="Temperature"] <- '#e31a1c'
col.pred$color[col.pred$predictor=="Salinity"] <- '#1f78b4'
col.pred$color[col.pred$predictor=="Depth"] <- '#6a3d9a'

bg <- c("gray60",'#e31a1c','#1f78b4','#6a3d9a')

col.pred[grep("chr",col.pred)] <- '#f1eef6' # non-candidate SNPs #not using this code chunk for now
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c')

# axes 1 & 2
plot(rda_oct_filters, type="n", scaling=3, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
points(rda_oct_filters, display="species", pch=21, cex=1, col="gray32", bg=col.pred$color, scaling=3)
#points(rda_oct_pc1, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(rda_oct_filters, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("Neutral", "Temperature","Salinity","Depth"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

##Using allele frequency matrix with only SNPs found in 300+ individuals (75%+)

indfreq300 <- np$load("pcangsd_freq_filters_snp300.indf.npy")
#400 individuals, 63,952 SNPs

indfreq300 <- as.data.frame(indfreq300)
indfreq300[1:5,1:5] #each individual has a row and each SNP has a column
#Need to transfer over the contig/position numbers for each SNP
sum(is.na(indfreq300)) #No NAs, so no imputation needed

#Load in mafs file from ANGSD to get contig/position numbers for each SNP
mafs300 <- read.table("freq_filters_snp300.mafs.txt", header = T)

#Combine chromo and position columns
mafs300$chromo_position <- paste(mafs300$chromo, mafs300$position, sep="_")
snp300_position <- as.vector(mafs300$chromo_position)

colnames(indfreq300) <- snp300_position
indfreq300[1:5, 1:5] #Now each column has the contig and position number as the header

#Check that the entire data frame looks okay by looking at a subset
indfreq300_subset <- indfreq300[-c(5001:63952)]
indfreq300_subset #Looks good

##RDA

rda_300 <- rda(indfreq300 ~ ., data=env_pred, scale=T)
rda_300
#                Inertia Proportion Rank
#Total         6.395e+04  1.000e+00
#Constrained   2.947e+03  4.608e-02    4
#Unconstrained 6.101e+04  9.539e-01   15
#Inertia is correlations

#Calculate adjusted R-squared
RsquareAdj(rda_300)
#r.squared=0.04608071 adj.r.squared=0.03642077

summary(eigenvals(rda_300, model = "constrained"))
#Most variance explained by RDA1

screeplot(rda_300)

plot(rda_300)

#Test for significance
signif.full <- anova.cca(rda_300, parallel=getOption("mc.cores"))

signif.axis <- anova.cca(rda_300, by="axis", parallel=getOption("mc.cores"))

##Identify candidate SNPs involved in local adaptation

load.rda300 <- scores(rda_300, choices=c(1:4), display="species")

hist(load.rda300[,1], main="Loadings on RDA1")
hist(load.rda300[,2], main="Loadings on RDA2")
hist(load.rda300[,3], main="Loadings on RDA3")
hist(load.rda300[,4], main="Loadings on RDA4")

cand1_300 <- outliers(load.rda300[,1],3) #0, but seems like all SNPs loading on the edges
cand2_300 <- outliers(load.rda300[,2],3) #2985
cand3_300 <- outliers(load.rda300[,3],3) #1005
cand4_300 <- outliers(load.rda300[,4],3) #296

ncand300 <- length(cand1_300) + length(cand2_300) + length(cand3_300) + length(cand4_300)
ncand300 #4286


##Run with temperature from NOAA October monthly means

#Load in environmental data

env_oct <- read.csv("Env_Data_Noaa_Oct.csv", header = T)

#Convert degrees F to degrees C

env_oct$noaa_oct_temp <- (env_oct$noaa_oct_temp_f - 32)*(5/9)

#Subset to lat, depth, noaa_oct_temp, and bottom_salinity for analysis
env_oct_pred <- subset(env_oct, select = c(lat, depth, noaa_oct_temp, bottom_salinity))

rda_oct_noaa <- rda(indfreq ~ ., data=env_oct_pred, scale=T)
rda_oct_noaa
#Inertia Proportion Rank
#Total         1.212e+05  1.000e+00
#Constrained   4.379e+03  3.614e-02    4
#Unconstrained 1.168e+05  9.639e-01   18
#Inertia is correlations

#Calculate adjujsted R-squared
RsquareAdj(rda_oct_noaa)
#r.squared=0.03614457 adj.r.squared=0.026384

summary(eigenvals(rda_oct_noaa, model = "constrained"))
#Most variance explained by RDA1

screeplot(rda_oct_noaa)

#Run RDA without latitude as a predictor

env_oct_nolat <- subset(env_oct_pred, select = c(depth, noaa_oct_temp, bottom_salinity))

rda_oct_nolat <- rda(indfreq ~ ., data=env_oct_nolat, scale=T)
rda_oct_nolat
#Inertia Proportion Rank
#Total         1.212e+05  1.000e+00
#Constrained   4.221e+03  3.484e-02    3
#Unconstrained 1.169e+05  9.652e-01   18
#Inertia is correlations
#Says RDA1 is still high, so maybe RDA1 is not latitude

plot(rda_oct_nolat)

##Individual allele frequencies with SNPs left after all filters are applied (22,600 SNPs)

indfreq_filters <- np$load("pcangsd_bam_list_realigned_mindp132_maxdp4000_minind0.indf.npy")
#400 individuals, 22,600 SNPs

indfreq_filters <- as.data.frame(indfreq_filters)
indfreq_filters[1:5,1:5] #each individual has a row and each SNP has a column
#Need to transfer over the contig/position numbers for each SNP
sum(is.na(indfreq_filters)) #No NAs, so no imputation needed

#Load in mafs file from ANGSD to get contig/position numbers for each SNP
mafs_filters <- read.table("bam_list_realigned_mindp132_maxdp4000_minind0.mafs.txt", header = T)
summary(mafs_filters)

#Remove SNPs with MAF<0.05 (this was already done during the indfreq matrix calculation, so need to do this for headers to match)
mafs_filters05 <- mafs_filters[-which(mafs_filters$knownEM>0.95 | mafs_filters$knownEM<0.05),]
#22,548 SNPs, so doesn't perfectly match the pcangsd output
#Will proceed without contig names for now

##RDA

library(vegan) #version 2.6-2

rda_filters <- rda(indfreq_filters ~ ., data=env_pred, scale=T)
rda_filters
#                Inertia Proportion Rank
#Total         2.257e+04  1.000e+00
#Constrained   4.498e+02  1.993e-02    4
#Unconstrained 2.212e+04  9.801e-01   44
#Inertia is correlations

#Calculate adjusted R-squared
RsquareAdj(rda_filters)
#r.squared=0.01992956 adj.r.squared=0.0100048

summary(eigenvals(rda_filters, model = "constrained"))
#Most variance explained by RDA1, lower eigenvalue than previously though

screeplot(rda_filters)

plot(rda_filters)

#Test for significance
signif.full <- anova.cca(rda_filters, parallel=getOption("mc.cores"))

signif.axis <- anova.cca(rda_filters, by="axis", parallel=getOption("mc.cores"))

##Identify candidate SNPs involved in local adaptation

load.rda.filters <- scores(rda_filters, choices=c(1:4), display="species")

hist(load.rda.filters[,1], main="Loadings on RDA1")
hist(load.rda.filters[,2], main="Loadings on RDA2")
hist(load.rda.filters[,3], main="Loadings on RDA3")
hist(load.rda.filters[,4], main="Loadings on RDA4")


cand1_filters <- outliers(load.rda.filters[,1],3) #0, but seems like all SNPs loading on the edges
cand2_filters <- outliers(load.rda.filters[,2],3) #464
cand3_filters <- outliers(load.rda.filters[,3],3) #497
cand4_filters <- outliers(load.rda.filters[,4],3) #23

ncand_filters <- length(cand1_filters) + length(cand2_filters) + length(cand3_filters) + length(cand4_filters)
ncand_filters #984

#Load indfreqs again, but with no minimum MAF filter

indfreq_filters_nomaf <- np$load("pcangsd_nominmaf_bam_list_realigned_mindp132_maxdp4000_minind0.indf.npy")
#400 individuals, still 22,600 SNPs

indfreq_filters_nomaf <- as.data.frame(indfreq_filters_nomaf)
indfreq_filters_nomaf[1:5,1:5] #each individual has a row and each SNP has a column
sum(is.na(indfreq_filters_nomaf)) #No NAs, so no imputation needed

#Combine chromo and position columns
mafs$chromo_position <- paste(mafs$chromo, mafs$position, sep="_")
snp_position <- as.vector(mafs$chromo_position)

colnames(indfreq) <- snp_position
indfreq[1:5, 1:5] #Now each column has the contig and position number as the header

#Check that the entire data frame looks okay by looking at a subset
indfreq_subset <- indfreq[-c(5001:121156)]
indfreq_subset #Looks good, each individual has a row, and each SNP has a column with the right header
