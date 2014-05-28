require(plyr)
require(reshape2)
require(ggplot2)

gwas2 = read.csv("~/Dropbox/HTA/Results/ProcessedData/GWAS2_complete.csv")
gwas3 = read.csv("~/Dropbox/HTA/Results/ProcessedData/GWAS3_complete.csv")

#################Control vs. control

gwas3 = gwas3[gwas3$drug == "control-lysate" & gwas3$assay == "a",-c(1,2,3,4,5,6,7,94:length(gwas3))]
gwas3[,3:length(gwas3)] = sapply(gwas3[,3:length(gwas3)], as.numeric)
gwas3 = gwas3[gwas3$col %% 2 == 1,]

mean2 = function(x){mean(x, na.rm=TRUE)}

gwas3mean = ddply(gwas3, .(row, col), numcolwise(mean2))
gwas3mean = gwas3mean[,-(1:2)]

df1 = melt(cor(gwas3mean, use = "complete.obs"))

controlsplot = ggplot(df1, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "red", mid = "green") + theme(text = element_text(size=4), axis.text.x=element_text(angle = 90, hjust = 0))
controlsplot


################Dact vs. dact

dact = gwas2[gwas2$drug == "dactinomycin" & gwas2$assay == "a",-c(1,2,3,4,5,6,7,94:length(gwas2))]
dact[,3:length(dact)] = sapply(dact[,3:length(dact)], as.numeric)
dact = dact[dact$col %% 2 == 1,]

mean2 = function(x){mean(x, na.rm=TRUE)}

dactmean = ddply(dact, .(row, col), numcolwise(mean2))
dactmean = dactmean[,-(1:2)]

df2 = melt(cor(dactmean, use = "complete.obs"))

dactplot = ggplot(df2, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "red", mid = "green") + theme(text = element_text(size=4), axis.text.x=element_text(angle = 90, hjust = 0))
dactplot

################Dact vs control

control = gwas2[gwas2$drug == "control-DMSO1" & gwas2$assay == "a",-c(1,2,3,4,5,6,7,94:length(gwas2))]
control[,3:length(control)] = sapply(control[,3:length(control)], as.numeric)
control = control[control$col %% 2 == 1,]
control = control[,-(1:2)]

df3 = melt(cor(dactmean, control, use = "complete.obs"))

dactcontrolplot = ggplot(df3, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient2(low = "purple", high = "red", mid = "green") + theme(text = element_text(size=4), axis.text.x=element_text(angle = 90, hjust = 0))
dactcontrolplot

