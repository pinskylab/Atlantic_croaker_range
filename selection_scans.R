##Selection Scans##

##PCAngsd Selection Scan

library(reticulate)
np <- import("numpy")

selection_pcangsd <- np$load("pcangsd_bam_list_realigned_mindp132_maxdp4000_minind0.selection.npy")

#Obtain p-values from selection scan

p <- pchisq(selection_pcangsd, 1, lower.tail=FALSE)

#Subset p-values <0.05

p <- as.data.frame(p)
p_sig <- subset(p, V1 < 0.05) #subsets to a data frame where SNPs with p <0.05 remain

#1,491 SNPs have p-values < 0.05, so ~6.6% of SNPs show evidence of selection (1,491/22,600)

