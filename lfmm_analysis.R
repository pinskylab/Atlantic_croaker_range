##GEA Analysis using LFMM
#Following tutorial by Breana Forester

##Individual frequencies using all filters applied dataset (~23K SNPs)

#See gea_analysis.R script for how genetic data was loaded and formatted

#Load in environmental data

env_oct_lfmm <- read.csv("Env_Data_Noaa_Oct.csv", header = T)

#Convert degrees F to degrees C

env_oct_lfmm$noaa_oct_temp <- (env_oct$noaa_oct_temp_f - 32)*(5/9)

#Subset to lat, depth, noaa_oct_temp, and bottom_salinity for analysis
env_oct_lfmm_pred <- subset(env_oct_lfmm, select = c(lat, depth, noaa_oct_temp, bottom_salinity))

#Perform a PCA on the environmental predictors (latitude, depth, temperature, salinity)
#and use the first principal component as a synthetic predictor

pred.pca <- rda(env_oct_lfmm_pred, scale=T)
summary(pred.pca)$cont
#PC1 explains 0.78 of the variance, PC2 0.15

screeplot(pred.pca, main = "Screeplot: Eigenvalues of Croaker Predictor Variables")
round(scores(pred.pca, choices=1:4, display="species", scaling=0), digits=3)

#Store synthetic PC axis predictor as pred.PC1 for use in LFMM
pred.PC1 <- scores(pred.pca, choices=1, display="sites", scaling=0)

#Determine K
#We’ll use a broken stick criterion to determine K. The broken stick stopping rule
#states that principal components should be retained as long as observed eigenvalues
#are higher than corresponding random broken stick components.

screeplot(pred.pca, main = "Screeplot of Croaker Predictor Variables with Broken Stick", bstick=TRUE, type="barplot")
#PC1 has higher eigenvalue than random broken stick components

#Now let’s run a PCA with the genomic data and
#plot the eigenvalues with the broken stick criterion:
#Already ran this in gea_analysis.R script, saved as pca_filters
pca_filters <- rda(indfreq_filters, scale=T)
screeplot(pca_filters, main = "Screeplot of Genetic Data with Broken Stick", bstick=TRUE, type="barplot")
#PC1 has higher eigenvalue than random broken stick components
#So K=2 (this follows the way a pRDA was run conditioning using PC1 + admixture results
#of PCs + 1 = K)

K <- 2

##Run LFMM
library(lfmm)
library(qvalue)

lfmm_oct_filters <- lfmm_ridge(Y=indfreq_filters, X=pred.PC1, K=K)

#Identify LFMM candidates using False Discovery Rate

lfmm_oct_filters.pv <- lfmm_test(Y=indfreq_filters, X=pred.PC1, lfmm=lfmm_oct_filters, calibrate="gif")

names(lfmm_oct_filters.pv) # this object includes raw z-scores and p-values,
#as well as GIF-calibrated scores and p-values

#Look at genomic inflation factor (GIF)
lfmm_oct_filters.pv$gif
#0.08, GIF less than one indicates the test may be too conservative

hist(lfmm_oct_filters.pv$pvalue[,1], main="Unadjusted p-values")
hist(lfmm_oct_filters.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

#Manually adjust the GIF correction factor:
zscore <- lfmm_oct_filters.pv$score[,1]   # zscores for first predictor, we only have one in our case...
(gif <- lfmm_oct_filters.pv$gif[1])       ## d.fault GIF for this predictor

new_gif1 <- 0.5

#Manual adjustment of the p-values:
adj.pv1 <- pchisq(zscore^2/new_gif1, df=1, lower = FALSE)

hist(lfmm_oct_filters.pv$pvalue[,1], main="Unadjusted p-values")
hist(lfmm_oct_filters.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")
hist(adj.pv1, main="REadjusted p-values")
#Doesn't seem to help much

lfmm_oct_filters.qv <- qvalue(lfmm_oct_filters.pv$calibrated.pvalue)$qvalues

length(which(lfmm_oct_filters.pv$calibrated.pvalue < 0.01)) #481

lfmm_oct_filters_p0.01 <- colnames(indfreq_filters)[which(lfmm_oct_filters.pv$calibrated.pvalue < 0.01)]
lfmm_oct_filters_p0.01 <- as.data.frame(lfmm_oct_filters_p0.01)
colnames(lfmm_oct_filters_p0.01) <- "snp"
lfmm_oct_filters_p0.01$calibrated_pvalue <- which(lfmm_oct_filters.pv$calibrated.pvalue < 0.01)

length(intersect(lfmm_oct_filters_p0.01$snp, cand_oct_filters$snp)) #137 candidates in both RDA and LFMM



##Running it with K=1
lfmm_oct_filters_k1 <- lfmm_ridge(Y=indfreq_filters, X=pred.PC1, K=1)

#Identify LFMM candidates using False Discovery Rate

lfmm_oct_filters_k1.pv <- lfmm_test(Y=indfreq_filters, X=pred.PC1, lfmm=lfmm_oct_filters_k1, calibrate="gif")

names(lfmm_oct_filters_k1.pv) # this object includes raw z-scores and p-values,
#as well as GIF-calibrated scores and p-values

#Look at genomic inflation factor (GIF)
lfmm_oct_filters_k1.pv$gif
#0.002, GIF less than one indicates the test may be too conservative

hist(lfmm_oct_filters_k1.pv$pvalue[,1], main="Unadjusted p-values")
hist(lfmm_oct_filters_k1.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

lfmm_oct_filters_k1.qv <- qvalue(lfmm_oct_filters_k1.pv$calibrated.pvalue)$qvalues

length(which(wolf.qv < 0.1))

##Run LFMM with each environmental predictor

#Temperature
lfmm_oct_filters_temp <- lfmm_ridge(Y=indfreq_filters, X=env_oct_lfmm$noaa_oct_temp, K=K)

#Identify LFMM candidates

lfmm_oct_filters_temp.pv <- lfmm_test(Y=indfreq_filters, X=env_oct_lfmm$noaa_oct_temp,
                                      lfmm=lfmm_oct_filters_temp, calibrate="gif")

#Look at genomic inflation factor (GIF)
lfmm_oct_filters_temp.pv$gif
#0.018, GIF less than one indicates the test may be too conservative

hist(lfmm_oct_filters_temp.pv$pvalue[,1], main="Unadjusted p-values")
hist(lfmm_oct_filters_temp.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

length(which(lfmm_oct_filters_temp.pv$calibrated.pvalue < 0.01)) #683

lfmm_oct_filters_p0.01_temp <- colnames(indfreq_filters)[which(lfmm_oct_filters_temp.pv$calibrated.pvalue < 0.01)]

length(intersect(lfmm_oct_filters_p0.01_temp, cand_oct_filters$snp)) #250 candidates in both RDA and temp LFMM

#Salinity
lfmm_oct_filters_salinity <- lfmm_ridge(Y=indfreq_filters, X=env_oct_lfmm$bottom_salinity, K=K)

#Identify LFMM candidates

lfmm_oct_filters_salinity.pv <- lfmm_test(Y=indfreq_filters, X=env_oct_lfmm$bottom_salinity,
                                      lfmm=lfmm_oct_filters_salinity, calibrate="gif")

#Look at genomic inflation factor (GIF)
lfmm_oct_filters_salinity.pv$gif
#2.96

hist(lfmm_oct_filters_salinity.pv$pvalue[,1], main="Unadjusted p-values")
hist(lfmm_oct_filters_salinity.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

length(which(lfmm_oct_filters_salinity.pv$calibrated.pvalue < 0.01)) #0

#Depth
lfmm_oct_filters_depth <- lfmm_ridge(Y=indfreq_filters, X=env_oct_lfmm$depth, K=K)

#Identify LFMM candidates

lfmm_oct_filters_depth.pv <- lfmm_test(Y=indfreq_filters, X=env_oct_lfmm$depth,
                                      lfmm=lfmm_oct_filters_depth, calibrate="gif")

#Look at genomic inflation factor (GIF)
lfmm_oct_filters_depth.pv$gif
#5.45

hist(lfmm_oct_filters_depth.pv$pvalue[,1], main="Unadjusted p-values")
hist(lfmm_oct_filters_depth.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

length(which(lfmm_oct_filters_depth.pv$calibrated.pvalue < 0.01)) #0

#Latitude
lfmm_oct_filters_lat <- lfmm_ridge(Y=indfreq_filters, X=env_oct_lfmm$lat, K=K)

#Identify LFMM candidates

lfmm_oct_filters_lat.pv <- lfmm_test(Y=indfreq_filters, X=env_oct_lfmm$lat,
                                       lfmm=lfmm_oct_filters_lat, calibrate="gif")

#Look at genomic inflation factor (GIF)
lfmm_oct_filters_lat.pv$gif
#0.14

hist(lfmm_oct_filters_lat.pv$pvalue[,1], main="Unadjusted p-values")
hist(lfmm_oct_filters_lat.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

length(which(lfmm_oct_filters_lat.pv$calibrated.pvalue < 0.01)) #592
lfmm_oct_filters_p0.01_lat <- colnames(indfreq_filters)[which(lfmm_oct_filters_lat.pv$calibrated.pvalue < 0.01)]

length(intersect(lfmm_oct_filters_p0.01_lat, cand_oct_filters$snp)) #190 candidates in both RDA and temp LFMM

##Using Copernicus data
#Env data is lat, depth, oxygen, temp, and salinity
env_ocean_filters_cand

#Perform a PCA on the environmental predictors (latitude, depth, oxygen, temperature, salinity)

pred.ocean.pca <- rda(env_ocean_filters_cand, scale=T)
summary(pred.ocean.pca)$cont
#PC1 explains 0.82 of the variance, PC2 0.11

round(scores(pred.ocean.pca, choices=1:5, display="species", scaling=0), digits=3)

#Store synthetic PC axis predictor as pred.PC1 for use in LFMM
pred.ocean.PC1 <- scores(pred.ocean.pca, choices=1, display="sites", scaling=0)

#Run LFMM with PC1 as a synthetic predictor
lfmm_ocean_filters <- lfmm_ridge(Y=indfreq_filters, X=pred.ocean.PC1, K=2)

#Identify LFMM candidates

lfmm_ocean_filters.pv <- lfmm_test(Y=indfreq_filters, X=pred.ocean.PC1,
                                        lfmm=lfmm_ocean_filters, calibrate="gif")

#Look at genomic inflation factor (GIF)
lfmm_ocean_filters.pv$gif
#0.12, GIF less than one indicates the test may be too conservative

hist(lfmm_ocean_filters.pv$pvalue[,1], main="Unadjusted p-values")
hist(lfmm_ocean_filters.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

lfmm_ocean_filters.qv <- qvalue(lfmm_ocean_filters.pv$calibrated.pvalue)$qvalues
length(which(lfmm_ocean_filters.qv < 0.1)) #381

lfmm_ocean_filters_q0.1 <- colnames(indfreq_filters)[which(lfmm_ocean_filters.qv < 0.1)]

length(intersect(lfmm_ocean_filters_q0.1, cand_ocean_filters$snp)) #47 in both RDA and LFMM with PC1


screeplot(pred.ocean.pca, main = "Screeplot of Croaker Predictor Variables with Broken Stick", bstick=TRUE, type="barplot")
#PC1 has higher eigenvalue than random broken stick components

##Run LFMM with each environmental predictor

##Temperature
lfmm_ocean_filters_temp <- lfmm_ridge(Y=indfreq_filters, X=env_ocean_filters_cand$temp, K=2)

#Identify LFMM candidates

lfmm_ocean_filters_temp.pv <- lfmm_test(Y=indfreq_filters, X=env_ocean_filters_cand$temp,
                                      lfmm=lfmm_ocean_filters_temp, calibrate="gif")

#Look at genomic inflation factor (GIF)
lfmm_ocean_filters_temp.pv$gif
#0.33, GIF less than one indicates the test may be too conservative

hist(lfmm_ocean_filters_temp.pv$pvalue[,1], main="Unadjusted p-values")
hist(lfmm_ocean_filters_temp.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

length(which(lfmm_ocean_filters_temp.pv$calibrated.pvalue < 0.01)) #359

lfmm_ocean_filters_p0.01_temp <- colnames(indfreq_filters)[which(lfmm_ocean_filters_temp.pv$calibrated.pvalue < 0.01)]

length(intersect(lfmm_ocean_filters_p0.01_temp, cand_ocean_filters$snp)) #16 candidates in both RDA and temp LFMM

lfmm_ocean_filters_temp.qv <- qvalue(lfmm_ocean_filters_temp.pv$calibrated.pvalue)$qvalues
length(which(lfmm_ocean_filters_temp.qv < 0.1)) #289

lfmm_ocean_filters_q0.1_temp <- colnames(indfreq_filters)[which(lfmm_ocean_filters_temp.qv < 0.1)]

length(intersect(lfmm_ocean_filters_q0.1_temp, cand_ocean_filters$snp)) #9 candidates in both RDA and temp LFMM

##Oxygen
lfmm_ocean_filters_oxygen <- lfmm_ridge(Y=indfreq_filters, X=env_ocean_filters_cand$oxygen, K=2)

#Identify LFMM candidates

lfmm_ocean_filters_oxygen.pv <- lfmm_test(Y=indfreq_filters, X=env_ocean_filters_cand$oxygen,
                                        lfmm=lfmm_ocean_filters_oxygen, calibrate="gif")

#Look at genomic inflation factor (GIF)
lfmm_ocean_filters_oxygen.pv$gif
#0.31, GIF less than one indicates the test may be too conservative

hist(lfmm_ocean_filters_oxygen.pv$pvalue[,1], main="Unadjusted p-values")
hist(lfmm_ocean_filters_oxygen.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

lfmm_ocean_filters_oxygen.qv <- qvalue(lfmm_ocean_filters_oxygen.pv$calibrated.pvalue)$qvalues
length(which(lfmm_ocean_filters_oxygen.qv < 0.1)) #122

lfmm_ocean_filters_q0.1_oxygen <- colnames(indfreq_filters)[which(lfmm_ocean_filters_oxygen.qv < 0.1)]

length(intersect(lfmm_ocean_filters_q0.1_oxygen, cand_ocean_filters$snp)) #9 candidates in both RDA and oxygen LFMM

##Salinity
lfmm_ocean_filters_salinity <- lfmm_ridge(Y=indfreq_filters, X=env_ocean_filters_cand$salinity, K=2)

#Identify LFMM candidates

lfmm_ocean_filters_salinity.pv <- lfmm_test(Y=indfreq_filters, X=env_ocean_filters_cand$salinity,
                                          lfmm=lfmm_ocean_filters_salinity, calibrate="gif")

#Look at genomic inflation factor (GIF)
lfmm_ocean_filters_salinity.pv$gif
#1.03

hist(lfmm_ocean_filters_salinity.pv$pvalue[,1], main="Unadjusted p-values")
hist(lfmm_ocean_filters_salinity.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

lfmm_ocean_filters_salinity.qv <- qvalue(lfmm_ocean_filters_salinity.pv$calibrated.pvalue)$qvalues
length(which(lfmm_ocean_filters_salinity.qv < 0.1)) #0

##Depth
lfmm_ocean_filters_depth <- lfmm_ridge(Y=indfreq_filters, X=env_ocean_filters_cand$depth, K=2)

#Identify LFMM candidates

lfmm_ocean_filters_depth.pv <- lfmm_test(Y=indfreq_filters, X=env_ocean_filters_cand$depth,
                                            lfmm=lfmm_ocean_filters_depth, calibrate="gif")

#Look at genomic inflation factor (GIF)
lfmm_ocean_filters_depth.pv$gif
#5.45

hist(lfmm_ocean_filters_depth.pv$pvalue[,1], main="Unadjusted p-values")
hist(lfmm_ocean_filters_depth.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

lfmm_ocean_filters_depth.qv <- qvalue(lfmm_ocean_filters_depth.pv$calibrated.pvalue)$qvalues
length(which(lfmm_ocean_filters_depth.qv < 0.1)) #0

##Latitude
lfmm_ocean_filters_lat <- lfmm_ridge(Y=indfreq_filters, X=env_ocean_filters_cand$lat, K=2)

#Identify LFMM candidates

lfmm_ocean_filters_lat.pv <- lfmm_test(Y=indfreq_filters, X=env_ocean_filters_cand$lat,
                                         lfmm=lfmm_ocean_filters_lat, calibrate="gif")

#Look at genomic inflation factor (GIF)
lfmm_ocean_filters_lat.pv$gif
#0.14

hist(lfmm_ocean_filters_lat.pv$pvalue[,1], main="Unadjusted p-values")
hist(lfmm_ocean_filters_lat.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

lfmm_ocean_filters_lat.qv <- qvalue(lfmm_ocean_filters_lat.pv$calibrated.pvalue)$qvalues
length(which(lfmm_ocean_filters_lat.qv < 0.1)) #460
lfmm_ocean_filters_q0.1_lat <- colnames(indfreq_filters)[which(lfmm_ocean_filters_lat.qv < 0.1)]

length(intersect(lfmm_ocean_filters_q0.1_lat, cand_ocean_filters$snp)) #90 SNPs in both RDA and lat LFMM
