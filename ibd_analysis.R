##Isolation-by-distance Analysis##

library(geodist)
library(tidyverse)

##Load in average lat and lon for each population

lat_lon_pop <- read.csv("Lat_Lon_byPop.csv")

colnames(lat_lon_pop) <- c('population', 'latitude', 'longitude')

#Make a vector for geodist

lat <- lat_lon_pop$latitude
lon <- lat_lon_pop$longitude

#Calculate geographic distance between populations using geodist package (outputs a matrix)

geo_distance <- geodist_vec(lon, lat, measure = "geodesic")

colnames(geo_distance) <- c("NJ", "DE", "MD", "VA", "NC", "OB", "LB", "SC", "GA", "FL")
rownames(geo_distance) <- c("NJ", "DE", "MD", "VA", "NC", "OB", "LB", "SC", "GA", "FL")

##Load in pairwise Fst matrix

fst <- read.csv("Pairwise_Fst.csv") #Not reading in row names correctly, will enter them manually
fst <- fst[,-1]
rownames(fst) <- c("NJ", "DE", "MD", "VA", "NC", "OB", "LB", "SC", "GA", "FL")

library(gdata)
#Make matrix symmetrical
upperTriangle(fst) <- lowerTriangle(fst, byrow=TRUE)

fst <- as.matrix(fst)

#Take ln of fst matrix

ln_fst <- log(fst)

library(vegan) #version 2.5-7

#Mantel test

mantel(geo_distance, ln_fst)
##r=-0.04563, significance=0.621, 999 permutations

#Get 95% CI
cor.test(geo_distance[lower.tri(geo_distance)],
         ln_fst[lower.tri(ln_fst)])
#-0.3346854 - 0.2512650

##Linear regression

y <- as.numeric(ln_fst)
x <- as.numeric(geo_distance)
mod_fst_geo <- lm(y~x)
summary(mod_fst_geo)
#Multiple R-squared:  0.002083,	Adjusted R-squared:  -0.009257
#F-statistic: 0.1836 on 1 and 88 DF,  p-value: 0.6693

plot(geo_distance/1000, ln_fst,
     xlab="Pairwise Geographic Distance (km)", ylab="ln(Fst)", pch=19)
abline(mod_fst_geo, col="red", lwd=2)


##Combining NJ and DE as one population (only 5 NJ individuals)

##Load in average lat and lon for each population

lat_lon_9pop <- read.csv("Lat_Lon_9Pop.csv")


#Make a vector for geodist

lat_9pop <- lat_lon_9pop$avg_lat
lon_9pop <- lat_lon_9pop$avg_lon

#Calculate geographic distance between populations using geodist package (outputs a matrix)

geo_distance_9pop <- geodist_vec(lon_9pop, lat_9pop, measure = "geodesic")

colnames(geo_distance_9pop) <- c("NJ_DE", "MD", "VA", "NC", "OB", "LB", "SC", "GA", "FL")
rownames(geo_distance_9pop) <- c("NJ_DE", "MD", "VA", "NC", "OB", "LB", "SC", "GA", "FL")

##Load in pairwise Fst matrix

fst_9pop <- read.csv("Pairwise_Fst_9pop.csv") #Not reading in row names correctly, will enter them manually
fst_9pop <- fst_9pop[,-1]
rownames(fst_9pop) <- c("NJ_DE", "MD", "VA", "NC", "OB", "LB", "SC", "GA", "FL")

fst_9pop <- as.matrix(fst_9pop)

library(gdata)
#Make matrix symmetrical
upperTriangle(fst_9pop) <- lowerTriangle(fst_9pop, byrow=TRUE)


#Take ln of fst matrix

ln_fst_9pop <- log(fst_9pop)

library(vegan) #version 2.5-7

#Mantel test

mantel(geo_distance_9pop, ln_fst_9pop)
##r=-0.1818, significance=0.793, 999 permutations

#Get 95% CI
cor.test(geo_distance_9pop[lower.tri(geo_distance_9pop)],
         ln_fst_9pop[lower.tri(ln_fst_9pop)])
#-0.4815939 to 0.1560294

##Linear regression

y9 <- as.numeric(ln_fst_9pop)
x9 <- as.numeric(geo_distance_9pop)
mod_fst_geo_9 <- lm(y9~x9)
summary(mod_fst_geo_9)
#Multiple R-squared:  0.03306,	Adjusted R-squared:  0.01925
#F-statistic: 2.393 on 1 and 70 DF,  p-value: 0.1264

plot(geo_distance_9pop/1000, ln_fst_9pop,
     xlab="Pairwise Geographic Distance (km)", ylab="ln(Fst)", pch=19)
abline(mod_fst_geo_9, col="red", lwd=2)

##Calculating in-water distances using the marmap package (following Benestan et al. 2021)

library(marmap)
library(rgdal) #Had to install rgdal and raster again for the getNOAA.bathy function to work
library(raster)
vignette("marmap-DataAnalysis")

#Need a matrix in format bathy first

colnames(lat_lon_9pop) <- c('population', 'lat', 'lon')
lat_lon_9pop_dist <- lat_lon_9pop[,-1]
lat_lon_9pop_dist <- data.frame(lat_lon_9pop$lon, lat_lon_9pop$lat)
colnames(lat_lon_9pop_dist) <- c('x', 'y')
lat_lon_9pop_dist$z <- c(-30, -30, -30, -30, -30, -30, -30, -30, -30)


eastcoast <- getNOAA.bathy(lon1=-85, lon2=-70, lat1=42, lat2=27, resolution=4)


#Compute transition matrix to be used by lc.dist() to compute least cost distances
#Limit possible paths to the continental shelf by setting max.depth=-200, min-depth=-10 to limit issues with path crossing land

trans1 <- trans.mat(eastcoast)
trans2 <- trans.mat(eastcoast, min.depth=10)
#min.depth=-10, max.depth=-200

#Calculate least-cost distances

lat_lon_9pop_dist <- lat_lon_9pop[,-1]
rownames(lat_lon_9pop_dist) <- c("NJ_DE", "MD", "VA", "NC", "OB", "LB", "SC", "GA", "FL")
colnames(lat_lon_9pop_dist) <- c("y", "x")

lat_lon_9pop_order <- data.frame(lat_lon_9pop_dist$x, lat_lon_9pop_dist$y)
colnames(lat_lon_9pop_order) <- c("x", "y")

out1 <- lc.dist(trans1, lat_lon_9pop_order, res="dist")
out1

out2 <- lc.dist(trans2, lat_lon_9pop_order, res="dist") #min.depth of 10 to avoid distances crossing land
out2

#Compare to geodist output
geo_distance_9pop_km <- geo_distance_9pop/1000
geo_distance_9pop_km #very similar, slightly longer on the in-water distance paths which makes sense

##Mantel test

mantel(out2, ln_fst_9pop)
##r=-0.1854, significance=0.798, 999 permutations

#Get 95% CI
cor.test(out2[lower.tri(out2)],
         ln_fst_9pop[lower.tri(ln_fst_9pop)])
#-0.4755269  0.2914712

##Linear regression

out2_matrix <- as.matrix(out2)
x9w <- as.numeric(out2_matrix)
y9w <- as.numeric(ln_fst_9pop)
mod_fst_geo_9w <- lm(y9w~x9w)
summary(mod_fst_geo_9w)
#Multiple R-squared:  0.03437,	Adjusted R-squared:  0.02057
#F-statistic: 2.491 on 1 and 70 DF,  p-value: 0.119

plot(out2_matrix, ln_fst_9pop,
     xlab="Pairwise Geographic Distance over Water (km)", ylab="ln(Fst)", pch=19)
abline(mod_fst_geo_9w, col="red", lwd=2)

##With linearized Fst (Fst/1-Fst)

fstlin <- as.matrix(fst_9pop/(1-fst_9pop))

##Mantel test

mantel(out2, fstlin)
##r=-0.111, significance=0.699, 999 permutations

##Linear regression

y9lin <- as.numeric(fstlin)
mod_fst_geo_9lin <- lm(y9lin~x9w)
summary(mod_fst_geo_9lin)
#Multiple R-squared:  0.01233,	Adjusted R-squared:  -0.001783
#F-statistic: 0.8736 on 1 and 70 DF,  p-value: 0.3532

plot(out2_matrix, fstlin,
     xlab="Pairwise Geographic Distance over Water (km)", ylab="Fst/(1 - Fst)", pch=19)
abline(mod_fst_geo_9lin, col="red", lwd=2)

##With Fst

##Mantel test

mantel(out2, fst_9pop)
##r=-0.113, significance=0.71, 999 permutations

##Linear regression

y9fst <- as.numeric(fst_9pop)
mod_fst_geo_9fst <- lm(y9fst~x9w)
summary(mod_fst_geo_9fst)
#Multiple R-squared:  0.01278,	Adjusted R-squared:  -0.001328
#F-statistic: 0.9058 on 1 and 70 DF,  p-value: 0.3445

plot(out2_matrix, fst_9pop,
     xlab="Pairwise Geographic Distance over Water (km)", ylab="Fst", pch=19)
abline(mod_fst_geo_9fst, col="red", lwd=2)


###Least-cost over-water distances between each pair of individuals

#Calculate least-cost distances

lat_lon_ind <- data.frame(lat_lon$V4, lat_lon$V3)
colnames(lat_lon_ind) <- c("x", "y")

dist_ind <- lc.dist(trans2, lat_lon_ind, res="dist") #min.depth of 10 to avoid distances crossing land
dist_ind_matrix <- as.matrix(dist_ind)

##Read in PCAngsd kinship matrix

library(reticulate)
np <- import("numpy")

kinship_pcangsd <- np$load("pcangsd_bam_list_realigned_mindp132_maxdp4000_minind0.kinship.npy")

#Mantel for individuals

mantel(dist_ind_matrix, kinship_pcangsd)
#r:-0.06626, p=1, 999 permutations

plot(dist_ind_matrix, kinship_pcangsd,
     xlab="Pairwise Geographic Distance over Water (km)", ylab="Kinship Coefficient", pch=19)

##Fine-scale IBD

require(ggplot2)
require(grid)
require(gridExtra)
require(reshape2)

#'-------------------------------------------------------------------------   @code

##Format matrices
rownames(kinship_pcangsd) <- seq(1:400)
colnames(kinship_pcangsd) <- seq(1:400)

rownames(dist_ind_matrix) <- seq(1:400)
colnames(dist_ind_matrix) <- seq(1:400)

dist_ind_table <- as.data.frame(as.table(dist_ind_matrix))
colnames(dist_ind_table) <- c("INDV1", "INDV2", "DIST")

kinship_pcangsd_table <- as.data.frame(as.table(kinship_pcangsd))
colnames(kinship_pcangsd_table) <- c("INDV1", "INDV2", "KINSHIP")

### Merge both datasets
geo_gen <- merge(dist_ind_table, kinship_pcangsd_table, by=c("INDV1","INDV2"))

### Check the number of pairwise observations in the minimum interval of 5km
min <- nrow(subset(geo_gen, geo_gen$DIST <5)) #5392

### Run the function
Slope_SES <- Slope.boot.SES(dataGen = geo_gen, window = seq(5,max(dataGen$DIST),by=5),min = min,repnull = 2,repboot = 2,core = 1)

### Save the results in a dataframe
Slope_SES_results <- data.frame(apply(lapply(Slope_SES, as.matrix), c(2, 3), mean),apply(lapply(Slope_SES, as.matrix), c(2, 3), sd)[,c(2,4)])
colnames(Slope_SES_results)[c(5,6)]<- c("slope_sd","SES_sd")

### Save the results in a file
write.table(Slope_pinksy_final_results, "results_ses_no_family_serranus.txt", quote=FALSE, row.names=FALSE)

### Save the results in a Rdata
save.image("SES_no_family_serranus.Rdata")
