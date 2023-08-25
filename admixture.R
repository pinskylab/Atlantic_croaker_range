##Admixture##

##Using outputs from NGSadmix, following code in tutorials from https://baylab.github.io/MarineGenomics/week-9--population-structure-using-ngsadmix.html#using-ngsadmix
## and http://www.popgen.dk/software/index.php/NgsAdmixTutorial

##Finding the best K using log-likelihood values from NGSadmix log outputs

library(stringr)

#Read in the NGSadmix logs

ngsadmix_logs <- list.files("ngsadmix/", pattern = ".log", full.names = T)
ngsadmix_logs

#Use lapply to read in all our log files at once

all_logs <- lapply(1:10, FUN = function(i) readLines(ngsadmix_logs[i]))

#this will pull out the line that starts with "b" from each file and return it as a list

foundset<-sapply(1:10, FUN= function(x) all_logs[[x]][which(str_sub(all_logs[[x]], 1, 1) == 'b')])
foundset

#now we need to pull out the first number in the string, we'll do this with the function sub

as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) )

#now lets store it in a dataframe
#make a dataframe with an index 1:7, this corresponds to our K values
logs <- data.frame(K=rep(1:10, each=1))

#add it to our likelihood values

logs$like<-as.vector(as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) ))

#and now we can calculate our delta K and probability

deltak <- tapply(logs$like, logs$K, FUN= function(x) (mean(x)/sd(x))) #not working as expected

#Based on log-likelihoods, K=1, then K=2 have the lowest values. K=10 is the highest.


##Admixture Plots##

#For K=2
#Load the covariance matrix
admix = read_table(paste0("ngsadmix/ngsadmix_bam_list_realigned_mindp132_maxdp4000_minind0_k2.qopt"), col_names = F)

##Add in pop labels

admix.id = as.data.frame(cbind(pop_label, admix))
names(admix.id) = c("pop","q1","q2")

plot = barplot(t(as.matrix(subset(admix.id, select=q1:q2))), col=1:2, border=NA)

#For K=10

admix_10 = read_table(paste0("ngsadmix/ngsadmix_bam_list_realigned_mindp132_maxdp4000_minind0_k10.qopt"), col_names = F)

##Add in pop labels

admix.id.10 = as.data.frame(cbind(pop_label, admix_10))
names(admix.id.10) = c("pop","q1","q2", "q3", "q4", "q5", "q6", "q7", "q8", "q9", "q10")

plot_10 = barplot(t(as.matrix(subset(admix.id.10, select=q1:q10))), col=1:10, border=NA)


#for K=3

admix_3 = read_table(paste0("ngsadmix/ngsadmix_bam_list_realigned_mindp132_maxdp4000_minind0_k3.qopt"), col_names = F)

##Add in pop labels

admix.id.3 = as.data.frame(cbind(pop_label, admix_3))
names(admix.id.3) = c("pop","q1","q2", "q3")

plot_3 = barplot(t(as.matrix(subset(admix.id.3, select=q1:q3))), col=1:3, border=NA)

#For K=4

admix_4 = read_table(paste0("ngsadmix/ngsadmix_bam_list_realigned_mindp132_maxdp4000_minind0_k4.qopt"), col_names = F)

##Add in pop labels

admix.id.4 = as.data.frame(cbind(pop_label, admix_4))
names(admix.id.4) = c("pop","q1","q2", "q3", "q4")

plot_4 = barplot(t(as.matrix(subset(admix.id.4, select=q1:q4))), col=1:4, border=NA)


##PCAngsd Admixture Resulsts##


library(reticulate)
np <- import("numpy")

admix_pcangsd <- np$load("pcangsd_bam_list_realigned_mindp132_maxdp4000_minind0_e.admix.Q.npy")


##Add in pop labels

pcadmix.id = as.data.frame(cbind(pop_label, admix_pcangsd))
names(pcadmix.id) = c("pop","q1","q2")

plot_pcadmix = barplot(t(as.matrix(subset(pcadmix.id, select=q1:q2))), names=pcadmix.id$pop, col=c("orange", "blue"), border=NA)


