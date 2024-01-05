#################################################################################################################################
#' Slope.boot.SES function  ---
#'
#' Function to compute relation between geographical distance and genetic
#'
#' For 1 iteration:  compute observed "lm(genetic ~ geographical distance)" and n "lm(genetic ~ geographical distance)" under null hypothesis
#' where genetic value are sample randomly
#'
#' Parameters
#' ---------------
#'  repboot : number of bootstrap
#'  repnull : number of iteration for null model
#'  window : value in km of the different sliding windows
#'  core: number of core for parallel computation
#'  min: minimum number of sample in a window
#'
#' Inputs
#' ---------------
#'  dataGen : data frame with pair, geo distance and genetics
#'
#'           INDV1       INDV2 DIST         Ajk
#'1: dip_bal_001 dip_bal_002  324 -0.00469689
#'2: dip_bal_001 dip_bal_003  324 -0.00465318
#'
#'
#'Outputs
#'---------------
#' list of dataframes (n = repboot) with value of the windows, slope, number of sample and Standardized effect Size
###################################################################################################################################

#'  -------------------------------------------------------------------------   @library
require(vegan)
require(dplyr)
require(data.table)
require(reshape)
require(parallel)
require(plyr)
library(pbapply)

Slope.boot.SES <- function (dataGen,window,min,repnull,repboot,core) {

  #core <- 1
  #min <- 5392
  #dataGen <- geo_gen
  #window <- seq(5,max(dataGen$DIST),by=5)
#repnull <- 2
  #repboot <- 2

  pblapply(1:repboot, function(j){

    corrGeoGen <- do.call(rbind, lapply(window,function(i) {

      results <- data.frame(dist = NA, slope = NA, n_comparisons = NA, SES = NA)[numeric(0), ]
      dataGen_boot <- subset(dataGen, DIST < i)
      dataGen_boot <- dataGen_boot[sample(nrow(dataGen_boot), min), ]

      #Null model on the selected iteration
      corr_GeoGen_random <- do.call(rbind, mclapply(1:repnull,function(k) {
        results_random <- data.frame(dist=NA, slope=NA, n_comparisons=NA)[numeric(0), ]
        dataGen_boot_random <- dataGen_boot
        dataGen_boot_random$KINSHIP <- sample(dataGen$KINSHIP,min)
        lm_random <- lm(dataGen_boot_random$KINSHIP~dataGen_boot_random$DIST)
        results_random <- rbind(results_random, data.frame(slope = lm_random$coefficients[2]))
        return(results_random) },mc.cores=core))

      meanNull <- apply(corr_GeoGen_random,2,mean)
      sdNull <- apply(corr_GeoGen_random,2,sd)

      # make a linear regression
      lm_boot <- lm(dataGen_boot$KINSHIP~dataGen_boot$DIST) #

      results <- rbind(results, data.frame(windows = i, slope = lm_boot$coefficients[2], n_comparisons = nrow(dataGen_boot), SES=(lm_boot$coefficients[2]-meanNull)/sdNull)) #add results to empty data frame that was created at beginning

      return(results)
    }))
  })
}
