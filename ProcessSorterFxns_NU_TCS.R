
#Function to read in txt file from sorter
readSorter <- function(file, tofmin=20, tofmax=2000, extmin=20, extmax=5000)  {
  data<-read.delim(file=file, header=T, na.strings=c("n/a"), as.is=T, stringsAsFactors=F)
  data<-data[!is.na(data$TOF),]
  data<-data[,!is.na(data[1,])]
  data<-data[data$TOF>=tofmin,]
  data<-data[data$EXT>=extmin,]
  data<-data[data$TOF<=tofmax,]
  data<-data[data$EXT<=extmax,]
  data$Column<-as.factor(data$Column)
  data$Row<-as.factor(data$Row)
  return(data)
}

#Function to extract sorter files from the L4 sort setup
procSetup <- function(file, tofmin=20, tofmax=2000, extmin=20, extmax=5000) {
  
  require(plyr)
  
  plate <- readSorter(file, tofmin, tofmax, extmin)
  modplate <- with(plate, data.frame(row=Row, col=as.factor(Column), 
                                     sort = Status.sort, TOF=TOF, EXT=EXT, 
                                     time=Time.Stamp, green=Green, yellow=Yellow, 
                                     red=Red))
  
  proc <- ddply(.data=modplate[modplate$row %in% c("A","B","C","D","E","F","G", "H"),], .var=c("row", "col"), .drop=F, .fun=function(x){
    c(pop = length(x$EXT), sorted = sum(x$sort==6), TOF = mean(x$TOF), EXT = mean(x$EXT))
  })
  
  
  return(proc)
  
}

#Function to make time per well go from 0 to X as opposed to running for the entire plate
extractTime <- function(x) {x$time <- x$time - min(x$time); return(x) }

#Function to take sorter raw dataframe and process to determine worm or bubble using SVM
sortertoDF <- function(file, tofmin=20, tofmax=2000, extmin=20, extmax=5000) {
  require(plyr)
  require(kernlab)
  plate <- readSorter(file, tofmin, tofmax, extmin, extmax)
  modplate <- with(plate, data.frame(row=Row, col=as.factor(Column), TOF=TOF, EXT=EXT, time=Time.Stamp, green=Green, yellow=Yellow, red=Red))
  modplate <- ddply(.data=modplate, .variables=c("row", "col"), .fun=extractTime)
  load("~/Dropbox/Biosort/Scripts and functions/bubbleSVMmodel_noProfiler.RData")
  plateprediction <- predict(bubbleSVMmodel_noProfiler, modplate[,3:8], type="probabilities")
  modplate$worm <- plateprediction[,"1"]
  modplate$call50 <- factor(as.numeric(modplate$worm>0.5), levels=c(1,0), labels=c("worm", "bubble"))
  modplate$norm.red <- modplate$red/modplate$TOF
  modplate$stage <- ifelse(modplate$TOF<=100, "L1", 
                           ifelse(modplate$TOF>100 & modplate$TOF<=250, "L2", 
                                  ifelse(modplate$TOF>250 & modplate$TOF<=400, "L3",
                                         ifelse(modplate$TOF>400 & modplate$TOF<=700, "L4", "adult"))))
  return(modplate)
}

#We might need a function to enter data without an SVM
sortertoDFnoSVM <- function(file, tofmin=20, tofmax=2000, extmin=20, extmax=5000) {
  require(plyr)
  #require(kernlab)
  plate <- readSorter(file, tofmin, tofmax, extmin, extmax)
  modplate <- with(plate, data.frame(row=Row, col=as.factor(Column), TOF=TOF, EXT=EXT, time=Time.Stamp, green=Green, yellow=Yellow, red=Red))
  modplate <- ddply(.data=modplate, .variables=c("row", "col"), .fun=extractTime)
  #load("~/Dropbox/Biosort/Scripts and functions/bubbleSVMmodel_noProfiler.RData")
  #plateprediction <- predict(bubbleSVMmodel_noProfiler, modplate[,3:8], type="probabilities")
  #modplate$worm <- plateprediction[,"1"]
  #modplate$call50 <- factor(as.numeric(modplate$worm>0.5), levels=c(1,0), labels=c("worm", "bubble"))
  modplate$norm.red <- modplate$red/modplate$TOF
  modplate$stage <- ifelse(modplate$TOF<=100, "L1", 
                           ifelse(modplate$TOF>=100 & modplate$TOF<=250, "L2", 
                                  ifelse(modplate$TOF>=250 & modplate$TOF<=400, "L3",
                                         ifelse(modplate$TOF>=400 & modplate$TOF<=700, "L4", "adult"))))
  return(modplate)
}

#Function to process worm data frame to phenotype data frame by row and column
processPheno <- function(modplate, strains) {
  require(plyr)
  processed <- ddply(.data=modplate[modplate$call50=="worm",], .variables=c("row", "col"), 
                     .fun=function(x){c(n=length(x$TOF), 
                                        meanTOF=mean(x$TOF, na.rm=T), 
                                        medianTOF=median(x$TOF, na.rm=T),
                                        minTOF=as.numeric(quantile(x$TOF)[1]),
                                        q05.t=as.numeric(quantile(x$TOF, probs=0.05)),
                                        q10.t=as.numeric(quantile(x$TOF, probs=0.1)[1]),
                                        q15.t=as.numeric(quantile(x$TOF, probs=0.15)[1]), 
                                        q20.t=as.numeric(quantile(x$TOF, probs=0.2)[1]),
                                        q25.t=as.numeric(quantile(x$TOF, probs=0.25)[1]),
                                        q30.t=as.numeric(quantile(x$TOF, probs=0.3)[1]),
                                        q35.t=as.numeric(quantile(x$TOF, probs=0.35)[1]),
                                        q40.t=as.numeric(quantile(x$TOF, probs=0.4)[1]),
                                        q45.t=as.numeric(quantile(x$TOF, probs=0.45)[1]),
                                        q55.t=as.numeric(quantile(x$TOF, probs=0.55)[1]),
                                        q60.t=as.numeric(quantile(x$TOF, probs=0.6)[1]),
                                        q65.t=as.numeric(quantile(x$TOF, probs=0.65)[1]),
                                        q70.t=as.numeric(quantile(x$TOF, probs=0.70)[1]),
                                        q75.t=as.numeric(quantile(x$TOF, probs=0.75)[1]),
                                        q80.t=as.numeric(quantile(x$TOF, probs=0.8)[1]),
                                        q85.t=as.numeric(quantile(x$TOF, probs=0.85)[1]),
                                        q90.t=as.numeric(quantile(x$TOF, probs=0.90)[1]),
                                        q95.t=as.numeric(quantile(x$TOF, probs=0.95)[1]),
                                        maxTOF=as.numeric(quantile(x$TOF)[5]),
                                        meanEXT=mean(x$EXT, na.rm=T), 
                                        medianEXT=median(x$EXT, na.rm=T),
                                        minEXT=as.numeric(quantile(x$EXT)[1]),
                                        q05.e=as.numeric(quantile(x$EXT, probs=0.05)[1]),
                                        q10.e=as.numeric(quantile(x$EXT, probs=0.1)[1]),
                                        q15.e=as.numeric(quantile(x$EXT, probs=0.15)[1]), 
                                        q20.e=as.numeric(quantile(x$EXT, probs=0.2)[1]),
                                        q25.e=as.numeric(quantile(x$EXT, probs=0.25)[1]),
                                        q30.e=as.numeric(quantile(x$EXT, probs=0.3)[1]),
                                        q35.e=as.numeric(quantile(x$EXT, probs=0.35)[1]),
                                        q40.e=as.numeric(quantile(x$EXT, probs=0.4)[1]),
                                        q45.e=as.numeric(quantile(x$EXT, probs=0.45)[1]),
                                        q55.e=as.numeric(quantile(x$EXT, probs=0.55)[1]),
                                        q60.e=as.numeric(quantile(x$EXT, probs=0.6)[1]),
                                        q65.e=as.numeric(quantile(x$EXT, probs=0.65)[1]),
                                        q70.e=as.numeric(quantile(x$EXT, probs=0.70)[1]),
                                        q75.e=as.numeric(quantile(x$EXT, probs=0.75)[1]),
                                        q80.e=as.numeric(quantile(x$EXT, probs=0.8)[1]),
                                        q85.e=as.numeric(quantile(x$EXT, probs=0.85)[1]),
                                        q90.e=as.numeric(quantile(x$EXT, probs=0.90)[1]),
                                        q95.e=as.numeric(quantile(x$EXT, probs=0.95)[1]),
                                        maxEXT=as.numeric(quantile(x$EXT)[5]),
                                        q25.logEXT=as.numeric(quantile(log(x$EXT), probs=0.25)[1]),
                                        median.logEXT=median(log(x$EXT), na.rm=T),
                                        mean.logEXT=mean(log(x$EXT), na.rm=T),
                                        q75.logEXT=as.numeric(quantile(log(x$EXT), probs=0.75)[1]),
                                        mean.red = mean(x$red, na.rm=T),
                                        q25.r = as.numeric(quantile(x$red, probs=0.25)[1]),
                                        median.red = median(x$red, na.rm=T),
                                        q75.r = as.numeric(quantile(x$red, probs=0.75)[1]),
                                        mean.normred = mean(x$norm.red, na.rm=T),
                                        median.normred = median(x$norm.red, na.rm=T),
                                        mean.gr = mean(x$green, na.rm=T),
                                        median.gr = median(x$green, na.rm=T),
                                        mean.y = mean(x$yellow, na.rm=T),
                                        median.y = median(x$yellow, na.rm=T),
                                        f.L1 = length(which(x$stage == "L1"))/length(x$stage),
                                        f.L2 = length(which(x$stage == "L2"))/length(x$stage),
                                        f.L3 = length(which(x$stage == "L3"))/length(x$stage),
                                        f.L4 = length(which(x$stage == "L4"))/length(x$stage),
                                        f.ad = length(which(x$stage == "adult"))/length(x$stage)
                                        )}, .drop=F)
  analysis <- data.frame(strain = as.character(strains), processed)
  analysis <- analysis[order(analysis$strain),]
  analysis <- analysis[order(analysis$row, analysis$col),]
  analysis$resid.meanTOF <- residuals(lm(analysis$meanTOF ~ analysis$n, na.action=na.exclude))
  analysis$resid.medianTOF <- residuals(lm(analysis$medianTOF ~ analysis$n, na.action=na.exclude))
  analysis$resid.q25TOF <- residuals(lm(analysis$q25.t ~ analysis$n, na.action=na.exclude))
  analysis$resid.q75TOF <- residuals(lm(analysis$q75.t ~ analysis$n, na.action=na.exclude))
  analysis$resid.meanEXT <- residuals(lm(analysis$meanEXT ~ analysis$n, na.action=na.exclude))
  analysis$resid.medianEXT <- residuals(lm(analysis$medianEXT ~ analysis$n, na.action=na.exclude))
  analysis$resid.q25EXT <- residuals(lm(analysis$q25.e ~ analysis$n, na.action=na.exclude))
  analysis$resid.q75EXT <- residuals(lm(analysis$q75.e ~ analysis$n, na.action=na.exclude))
  analysis$norm.red.q25 <- analysis$q25.r/analysis$q25.t
  analysis$norm.red.median <- analysis$median.red/analysis$medianTOF
  analysis$norm.red.q75 <- analysis$q75.r/analysis$q75.t
  return(analysis)
}

#Function to process worm data frame to phenotype data frame by row and column also add concentration
processPhenoConc <- function(modplate, strains, conc) {
  require(plyr)
  processed <- ddply(.data=modplate[modplate$call50=="worm",], .variables=c("row", "col"), 
                     .fun=function(x){c(n=length(x$TOF), 
                                        meanTOF=mean(x$TOF), 
                                        medianTOF=median(x$TOF),
                                        minTOF=as.numeric(quantile(x$TOF)[1]),
                                        q05=as.numeric(quantile(x$TOF, probs=0.05)),
                                        q10=as.numeric(quantile(x$TOF, probs=0.1)[1]),
                                        q15=as.numeric(quantile(x$TOF, probs=0.15)[1]), 
                                        q20=as.numeric(quantile(x$TOF, probs=0.2)[1]),
                                        q25=as.numeric(quantile(x$TOF, probs=0.25)[1]),
                                        q30=as.numeric(quantile(x$TOF, probs=0.3)[1]),
                                        q35=as.numeric(quantile(x$TOF, probs=0.35)[1]),
                                        q40=as.numeric(quantile(x$TOF, probs=0.4)[1]),
                                        q45=as.numeric(quantile(x$TOF, probs=0.45)[1]),
                                        q55=as.numeric(quantile(x$TOF, probs=0.55)[1]),
                                        q60=as.numeric(quantile(x$TOF, probs=0.6)[1]),
                                        q65=as.numeric(quantile(x$TOF, probs=0.65)[1]),
                                        q70=as.numeric(quantile(x$TOF, probs=0.70)[1]),
                                        q75=as.numeric(quantile(x$TOF, probs=0.75)[1]),
                                        q80=as.numeric(quantile(x$TOF, probs=0.8)[1]),
                                        q85=as.numeric(quantile(x$TOF, probs=0.85)[1]),
                                        q90=as.numeric(quantile(x$TOF, probs=0.90)[1]),
                                        q95=as.numeric(quantile(x$TOF, probs=0.95)[1]),
                                        maxTOF=as.numeric(quantile(x$TOF)[5]),
                                        meanEXT=mean(x$EXT), 
                                        medianEXT=median(x$EXT)
                     )}, .drop=F)
  analysis <- data.frame(strain = as.character(strains), processed)
  analysis <- analysis[order(analysis$strain),]
  analysis <- analysis[order(analysis$row, analysis$col),]
  analysis <- droplevels(na.omit(analysis))
  analysis$conc <- ifelse(analysis$col %in% c(1, 7), conc[1], 
                          ifelse(analysis$col %in% c(2,8), conc[2], 
                                 ifelse(analysis$col %in% c(3,9), conc[3], 
                                        ifelse(analysis$col %in% c(4,10), conc[4],
                                               ifelse(analysis$col %in% c(5,11), conc[5], NA)))))
  return(analysis)
}


#Process by strain
processFullbyStrain <- function(results){
  ddply(.data=results, .var=c("strain"), .fun=function(x) {c(broodmedian=median(x$n),
                                                             broodmean=mean(x$n),
                                                             TOFmean=median(x$meanTOF), TOFmedian=median(x$medianTOF), 
                                                             EXTmean=median(x$meanEXT), EXTmedian=median(x$medianEXT))})
}

#Process by strain and concentration for dose response summary
processFullbyStrainConc <- function(results){
  ddply(.data=results, .var=c("strain", "conc"), .fun=function(x) {c(broodmedian=median(x$n), TOFmean=median(x$meanTOF), TOFmedian=median(x$medianTOF),
                                                                     TOFmin=median(x$minTOF), TOFmax=median(x$maxTOF), TOF25=median(x$q25), TOF75=median(x$q75),
                                                                     broodmean = mean(x$n))})
}
#Strain lists need to be expanded depending on the setup orientation on plate

#Expand strain list to seven replicates with blank by row H
expandStrainlist_7repsbyCol <- function(setup){
  
  new.setup <- matrix(ncol=12, nrow=8)
  
  for (i in c(1:12)) {
    
    new.setup[,i]<-setup[i]
    new.setup[8,]<- NA
    
  }
  
  return(c(t(new.setup)))
  
}

#Expand strain list to five replicates with blanks in columns 6 and 12
expandStrainlist_5repsperStrain <- function(setup, num.strains){
  
  new.setup <- matrix(ncol=num.strains, nrow=6)
  
  for (i in c(1:num.strains)) {
    
    new.setup[,i]<-setup[i]
    new.setup[6,]<- NA
    
  }
  
  return(c(new.setup))
  
}

#Expand strain list to three replicates with blanks in columns 4, 8, and 12
expandStrainlist_3repsperStrain <- function(setup, num.strains){
  
  new.setup <- matrix(ncol=num.strains, nrow=4)
  
  for (i in c(1:num.strains)) {
    
    new.setup[,i]<-setup[i]
    new.setup[4,]<- NA
    
  }
  
  return(c(new.setup))
  
}

#Expand strain list to 48 strains per plate
expandStrainlist_48strains <- function(setup){
  
  new.setup <- matrix(ncol=length(setup), nrow=2)
  
  for (i in c(1:length(setup))) {
    
    new.setup[,i]<-setup[i]
    new.setup[2,]<- NA
    
  }
  
  return(c(new.setup))
  
}

expandStrainlist_dose <- function(setup){
  
  new.setup <- matrix(ncol=16, nrow=6)
  
  for (i in 1:length(setup)) {
    
    new.setup[,i]<-setup[i]
    
    
  }
  
  return(c(new.setup))
  
}

#Function to remove contaminated wells. #81 columns is the size of the dataframe

removeWells <- function(proc, badwells) {
  sp.bw <- str_split(badwells, "", 3)
  for (i in seq(1, length(sp.bw))) {
    row <- sp.bw[[i]][2]
    col <- sp.bw[[i]][3]
    proc[which(proc$row == row & proc$col == col),-(1:7)] <- NA
  }
  return(proc)
}

#Function to make sure all wells were recorded by the sorter software

checkData <- function(sorterDF){
  complete <- paste(rep(c("A", "B", "C", "D", "E", "F", "G", "H"), each =12), rep(seq(1,12), 8), sep="")
  row.col <- paste(sorterDF$row, sorterDF$col, sep="")
  num <- length(unique(row.col))
  missing <- complete[is.na(match(complete, row.col))]
  return(missing)
}