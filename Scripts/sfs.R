##Plotting Site Frequency Spectrum for each population##


sfs_nj <- scan("NJ.sfs")
barplot(sfs_nj[-1])

sfs_de <- scan("DE1.sfs")
barplot(sfs_de[-1])

sfs_md <- scan("MD.sfs")
barplot(sfs_md[-1])

sfs_va <- scan("VA.sfs")
barplot(sfs_va[-1])

sfs_nc <- scan("NC.sfs")
barplot(sfs_nc[-1])

sfs_ob <- scan("OB.sfs")
barplot(sfs_ob[-1])

sfs_lb <- scan("LB.sfs")
barplot(sfs_lb[-1])

sfs_sc <- scan("SC.sfs")
barplot(sfs_sc[-1])

sfs_ga <- scan("GA.sfs")
barplot(sfs_ga[-1])

sfs_fl <- scan("FL.sfs")
barplot(sfs_fl[-1])

##DE_NJ combined population
sfs_nj_de <- scan("DE.sfs")
barplot(sfs_nj_de[-1])
