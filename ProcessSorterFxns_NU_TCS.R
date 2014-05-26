#Extract the plate number
plateno <- function(string)
{
    split<-strsplit(string,"_")
    split<-split[[1]]
    split<-split[1]
    
    split<-strsplit(split,"p")
    split<-split[[1]]
    plate<-split[2]
}

experimentName <- function(filePath){
    splitfp <- strsplit(filePath,"/")
    dirName <- splitfp[[1]][(length(splitfp[[1]]))]
    details <- strsplit(dirName,"_")[[1]][2]
    return(details)
}

#Extract the metadata info
info <- function(filePath){
    splitfp <- strsplit(filePath,"/")
    dirName <- splitfp[[1]][(length(splitfp[[1]])-1)]
    
    date <- strsplit(dirName,"_")[[1]][1]
    
    details <- strsplit(dirName,"_")[[1]][2]
    
    experiment <- strsplit(details,"[0-9]+")[[1]][1]
    round <- strsplit(details,"(?i)[a-z]+")[[1]][2]
    assay <- strsplit(details,"[0-9]+")[[1]][2]
    
    split <- strsplit(splitfp[[1]][(length(splitfp[[1]]))],"_")[[1]]
    drug <- strsplit(split[2],"\\.")[[1]][1]
    plate <- strsplit(split[1],"p")[[1]][2]
    
    frame <- data.frame(date,experiment,round,assay,plate,drug)
    
    return(frame)
}

#Function to process the setup data
procSetup <- function(file, tofmin=60, tofmax=2000, extmin=20, extmax=5000) {
    
    #Read in the sorter data from the file
    plate <- readSorter(file, tofmin, tofmax, extmin)
    modplate <- with(plate, data.frame(row=Row, col=as.factor(Column),
                                       sort = Status.sort, TOF=TOF, EXT=EXT,
                                       time=Time.Stamp, green=Green, yellow=Yellow,
                                       red=Red))
    
    #Process the data
    proc <- ddply(.data=modplate[modplate$row %in% c("A","B","C","D","E","F","G", "H"),], .var=c("row", "col"), .drop=F, .fun=function(x){
        c(pop = length(x$EXT), sorted = sum(x$sort==6), TOF = mean(x$TOF), EXT = mean(x$EXT),TOFmed=median(x$TOF),EXTmed=median(x$EXT)
        )})
    
    return(proc)
}

#Convert the sorter data to a data frame
sortertoDF <- function(file, tofmin=60, tofmax=2000, extmin=20, extmax=5000) {
    plate <- readSorter(file, tofmin, tofmax, extmin, extmax)
    modplate <- with(plate, data.frame(row=Row, col=as.factor(Column), TOF=TOF, EXT=EXT, time=Time.Stamp, green=Green, yellow=Yellow, red=Red))
    modplate <- ddply(.data=modplate, .variables=c("row", "col"), .fun=extractTime)
    load(file.path(dir.root,"Dropbox/Biosort/Scripts and functions/bubbleSVMmodel_noProfiler.RData"))
    plateprediction <- predict(bubbleSVMmodel_noProfiler, modplate[,3:8], type="probabilities")
    modplate$worm <- plateprediction[,"1"]
    modplate$call50 <- factor(as.numeric(modplate$worm>0.5), levels=c(1,0), labels=c("worm", "bubble"))
    modplate$norm.red <- modplate$red/modplate$TOF
    modplate$stage <- ifelse(modplate$TOF<=100, "L1", 
                             ifelse(modplate$TOF>=100 & modplate$TOF<=250, "L2", 
                                    ifelse(modplate$TOF>=250 & modplate$TOF<=400, "L3",
                                           ifelse(modplate$TOF>=400 & modplate$TOF<=700, "L4", "adult"))))
    return(modplate)
}

#added quantiles of minEXT, q05_EXT through q95_EXT, maxEXT
processPheno <- function(modplate, strains) {
    processed <- ddply(.data=modplate[modplate$call50=="worm" | modplate$TOF == -1,], .variables=c("row", "col"),
                       .fun=function(x){c(n=length(x$TOF),
                                          meanTOF=mean(x$TOF),
                                          medianTOF=median(x$TOF),
                                          minTOF=as.numeric(quantile(x$TOF)[1]),
                                          q05_TOF=as.numeric(quantile(x$TOF, probs=0.05)),
                                          q10_TOF=as.numeric(quantile(x$TOF, probs=0.1)[1]),
                                          q15_TOF=as.numeric(quantile(x$TOF, probs=0.15)[1]),
                                          q20_TOF=as.numeric(quantile(x$TOF, probs=0.2)[1]),
                                          q25_TOF=as.numeric(quantile(x$TOF, probs=0.25)[1]),
                                          q30_TOF=as.numeric(quantile(x$TOF, probs=0.3)[1]),
                                          q35_TOF=as.numeric(quantile(x$TOF, probs=0.35)[1]),
                                          q40_TOF=as.numeric(quantile(x$TOF, probs=0.4)[1]),
                                          q45_TOF=as.numeric(quantile(x$TOF, probs=0.45)[1]),
                                          q55_TOF=as.numeric(quantile(x$TOF, probs=0.55)[1]),
                                          q60_TOF=as.numeric(quantile(x$TOF, probs=0.6)[1]),
                                          q65_TOF=as.numeric(quantile(x$TOF, probs=0.65)[1]),
                                          q70_TOF=as.numeric(quantile(x$TOF, probs=0.70)[1]),
                                          q75_TOF=as.numeric(quantile(x$TOF, probs=0.75)[1]),
                                          q80_TOF=as.numeric(quantile(x$TOF, probs=0.8)[1]),
                                          q85_TOF=as.numeric(quantile(x$TOF, probs=0.85)[1]),
                                          q90_TOF=as.numeric(quantile(x$TOF, probs=0.90)[1]),
                                          q95_TOF=as.numeric(quantile(x$TOF, probs=0.95)[1]),
                                          maxTOF=as.numeric(quantile(x$TOF)[5]),
                                          meanEXT=mean(x$EXT),
                                          medianEXT=median(x$EXT),
                                          
                                          minEXT=as.numeric(quantile(x$EXT)[1]),
                                          q05_EXT=as.numeric(quantile(x$EXT, probs=0.05)),
                                          q10_EXT=as.numeric(quantile(x$EXT, probs=0.1)[1]),
                                          q15_EXT=as.numeric(quantile(x$EXT, probs=0.15)[1]),
                                          q20_EXT=as.numeric(quantile(x$EXT, probs=0.2)[1]),
                                          q25_EXT=as.numeric(quantile(x$EXT, probs=0.25)[1]),
                                          q30_EXT=as.numeric(quantile(x$EXT, probs=0.3)[1]),
                                          q35_EXT=as.numeric(quantile(x$EXT, probs=0.35)[1]),
                                          q40_EXT=as.numeric(quantile(x$EXT, probs=0.4)[1]),
                                          q45_EXT=as.numeric(quantile(x$EXT, probs=0.45)[1]),
                                          q55_EXT=as.numeric(quantile(x$EXT, probs=0.55)[1]),
                                          q60_EXT=as.numeric(quantile(x$EXT, probs=0.6)[1]),
                                          q65_EXT=as.numeric(quantile(x$EXT, probs=0.65)[1]),
                                          q70_EXT=as.numeric(quantile(x$EXT, probs=0.70)[1]),
                                          q75_EXT=as.numeric(quantile(x$EXT, probs=0.75)[1]),
                                          q80_EXT=as.numeric(quantile(x$EXT, probs=0.8)[1]),
                                          q85_EXT=as.numeric(quantile(x$EXT, probs=0.85)[1]),
                                          q90_EXT=as.numeric(quantile(x$EXT, probs=0.90)[1]),
                                          q95_EXT=as.numeric(quantile(x$EXT, probs=0.95)[1]),
                                          maxEXT=as.numeric(quantile(x$EXT)[5]),
                                          mean.red = mean(x$red, na.rm=T),
                                          median.red = median(x$red, na.rm=T),
                                          mean.gr = mean(x$green, na.rm=T),
                                          median.gr = median(x$green, na.rm=T),
                                          mean.y = mean(x$yellow, na.rm=T),
                                          median.y = median(x$yellow, na.rm=T),
                                          mean.normred = mean(x$norm.red, na.rm=T),
                                          median.normred = mean(x$norm.red, na.rm=T),
                                          f.L1 = length(which(x$stage == "L1"))/length(x$stage),
                                          f.L2 = length(which(x$stage == "L2"))/length(x$stage),
                                          f.L3 = length(which(x$stage == "L3"))/length(x$stage),
                                          f.L4 = length(which(x$stage == "L4"))/length(x$stage),
                                          f.ad = length(which(x$stage == "adult"))/length(x$stage),
                                          
                                          log.meanEXT=log(mean(x$EXT)),
                                          log.medianEXT=log(median(x$EXT)),
                                          log.minEXT=as.numeric(log(quantile(x$EXT)[1])),
                                          log.q05_EXT=as.numeric(log(quantile(x$EXT, probs=0.05))),
                                          log.q10_EXT=as.numeric(log(quantile(x$EXT, probs=0.1)[1])),
                                          log.q15_EXT=as.numeric(log(quantile(x$EXT, probs=0.15)[1])),
                                          log.q20_EXT=as.numeric(log(quantile(x$EXT, probs=0.2)[1])),
                                          log.q25_EXT=as.numeric(log(quantile(x$EXT, probs=0.25)[1])),
                                          log.q30_EXT=as.numeric(log(quantile(x$EXT, probs=0.3)[1])),
                                          log.q35_EXT=as.numeric(log(quantile(x$EXT, probs=0.35)[1])),
                                          log.q40_EXT=as.numeric(log(quantile(x$EXT, probs=0.4)[1])),
                                          log.q45_EXT=as.numeric(log(quantile(x$EXT, probs=0.45)[1])),
                                          log.q55_EXT=as.numeric(log(quantile(x$EXT, probs=0.55)[1])),
                                          log.q60_EXT=as.numeric(log(quantile(x$EXT, probs=0.6)[1])),
                                          log.q65_EXT=as.numeric(log(quantile(x$EXT, probs=0.65)[1])),
                                          log.q70_EXT=as.numeric(log(quantile(x$EXT, probs=0.70)[1])),
                                          log.q75_EXT=as.numeric(log(quantile(x$EXT, probs=0.75)[1])),
                                          log.q80_EXT=as.numeric(log(quantile(x$EXT, probs=0.8)[1])),
                                          log.q85_EXT=as.numeric(log(quantile(x$EXT, probs=0.85)[1])),
                                          log.q90_EXT=as.numeric(log(quantile(x$EXT, probs=0.90)[1])),
                                          log.q95_EXT=as.numeric(log(quantile(x$EXT, probs=0.95)[1])),
                                          log.maxEXT=as.numeric(log(quantile(x$EXT)[5])),
                                          log.mean.red = log(mean(x$red, na.rm=T)),
                                          log.median.red = log(median(x$red, na.rm=T)),
                                          log.mean.normred = log(mean(x$norm.red, na.rm=T)),
                                          log.median.normred = log(mean(x$norm.red, na.rm=T))
                       )}, .drop=F)
    
    analysis <- data.frame(strain = as.character(strains), processed)
    analysis <- analysis[order(analysis$strain),]
    analysis <- analysis[order(analysis$row, analysis$col),]
    analysis[as.numeric(analysis$meanTOF)==-1 | is.na(analysis$meanTOF),4:ncol(analysis)] = NA
    #analysis <- droplevels(na.omit(analysis))
    return(analysis)
}

meltdf <- function(score){
    newscore<-data.frame(row=rep(score$row,each=1),col=rep(score$col,each=1),n=rep(score$n,each=1),f.L1=rep(score$f.L1,each=1),f.L2=rep(score$f.L2,each=1),f.L3=rep(score$f.L3,each=1),f.L4=rep(score$f.L4,each=1),f.ad=rep(score$f.ad,each=1))
    newscore<-melt(newscore,id.var=c("row","col","n"))
    return(newscore)
}

possContam <- function(procDataFrame){
    strainMean <- mean(procDataFrame$n[!is.na(procDataFrame$strain)], na.rm = TRUE)
    strainSD <- sd(procDataFrame$n[!is.na(procDataFrame$strain)], na.rm = TRUE)
    washMean <- mean(procDataFrame$n[is.na(procDataFrame$strain)], na.rm = TRUE)
    washSD <- sd(procDataFrame$n[is.na(procDataFrame$strain)], na.rm = TRUE)
    possibleContam <- c()
    for(j in seq(1,nrow(procDataFrame),)){
        if(!is.na(procDataFrame[j,"strain"] & !is.na(procDataFrame[j,"n"]))){
            if(procDataFrame[j,"n"] > strainMean + (2*strainSD)){
                row <- as.character(procDataFrame[j, "row"])
                col <- as.numeric(procDataFrame[j, "col"])
                adjacentWash <- procDataFrame[procDataFrame$row==row & procDataFrame$col==(col+1),"n"]
                if(adjacentWash > washMean + (2*washSD)){
                    possibleContam <- append(possibleContam, paste0(row, col))
                }
            }
        }
    } 
}












#Function to read in txt file from sorter
readSorter <- function(file, tofmin=60, tofmax=2000, extmin=20, extmax=5000)  {
  data <- read.delim(file=file, header=T, na.strings=c("n/a"), as.is=T, stringsAsFactors=F)
  data <- data[!is.na(data$TOF),]
  data <- data[,!is.na(data[1,])]
  data <- data[(data$TOF>=tofmin & data$TOF<=tofmax) | data$TOF == -1,]
  data <- data[(data$EXT>=extmin & data$EXT<=extmax) | data$EXT == -1,]
  data$Column <- as.factor(data$Column)
  data$Row <- as.factor(data$Row)
  return(data)
}


#Function to make time per well go from 0 to X as opposed to running for the entire plate
extractTime <- function(x) {x$time <- x$time - min(x$time); return(x) }

#Function to take sorter raw dataframe and process to determine worm or bubble using SVM
sortertoDF <- function(file, tofmin=60, tofmax=2000, extmin=20, extmax=5000) {
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
sortertoDFnoSVM <- function(file, tofmin=60, tofmax=2000, extmin=20, extmax=5000) {
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

#Function to remove contaminated wells

removeWells <- function(proc, badwells) {
  sp.bw <- str_split(badwells, "", 3)
  for (i in seq(1, length(sp.bw))) {
    row <- sp.bw[[i]][2]
    col <- sp.bw[[i]][3]
    proc[which(proc$row == row & proc$col == col),-(1:3)] <- NA
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