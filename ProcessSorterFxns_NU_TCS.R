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
info <- function(filePath, levels = 1){
    splitfp <- strsplit(filePath,"/")
    dirName <- splitfp[[1]][(length(splitfp[[1]])-levels)]
    
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

meltdf <- function(score){
    newscore<-data.frame(row=rep(score$row,each=1),col=rep(score$col,each=1),n=rep(score$n,each=1),f.L1=rep(score$f.L1,each=1),f.L2L3=rep(score$f.L2L3,each=1),f.L4=rep(score$f.L4,each=1),f.ad=rep(score$f.ad,each=1))
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



summarizePlate_worms <- function(plate, strains=NULL, quantiles=FALSE, log=FALSE, ends=FALSE) {
    plate <- plate[plate$call50=="object" | plate$TOF == -1 | is.na(plate$call50),]
    plate <- fillWells(plate)
    processed <- plate %>% group_by(row, col) %>% summarise(n=length(TOF),
                                                            
                                                            mean.TOF=mean(TOF, na.rm=TRUE),
                                                            min.TOF=as.numeric(quantile(TOF, na.rm=TRUE)[1]),
                                                            q10.TOF=as.numeric(quantile(TOF, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.TOF=as.numeric(quantile(TOF, probs=0.25, na.rm=TRUE)[1]),
                                                            median.TOF=median(TOF, na.rm=TRUE),
                                                            q75.TOF=as.numeric(quantile(TOF, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.TOF=as.numeric(quantile(TOF, probs=0.90, na.rm=TRUE)[1]),
                                                            max.TOF=as.numeric(quantile(TOF, na.rm=TRUE)[5]),
                                                            
                                                            mean.EXT=mean(EXT, na.rm=TRUE),
                                                            min.EXT=as.numeric(quantile(EXT, na.rm=TRUE)[1]),
                                                            q10.EXT=as.numeric(quantile(EXT, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.EXT=as.numeric(quantile(EXT, probs=0.25, na.rm=TRUE)[1]),
                                                            median.EXT=median(EXT, na.rm=TRUE),
                                                            q75.EXT=as.numeric(quantile(EXT, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.EXT=as.numeric(quantile(EXT, probs=0.90, na.rm=TRUE)[1]),
                                                            max.EXT=as.numeric(quantile(EXT, na.rm=TRUE)[5]),
                                                            
                                                            mean.red=mean(red, na.rm=TRUE),
                                                            min.red=as.numeric(quantile(red, na.rm=TRUE)[1]),
                                                            q10.red=as.numeric(quantile(red, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.red=as.numeric(quantile(red, probs=0.25, na.rm=TRUE)[1]),
                                                            median.red=median(red, na.rm=TRUE),
                                                            q75.red=as.numeric(quantile(red, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.red=as.numeric(quantile(red, probs=0.9, na.rm=TRUE)[1]),
                                                            max.red=as.numeric(quantile(red, na.rm=TRUE)[5]),
                                                            
                                                            mean.green=mean(green, na.rm=TRUE),
                                                            median.green=median(green, na.rm=TRUE),
                                                            
                                                            mean.yellow=mean(yellow, na.rm=TRUE),
                                                            median.yellow=median(yellow, na.rm=TRUE),
                                                            
                                                            mean.norm.EXT=mean(norm.EXT, na.rm=TRUE),
                                                            min.norm.EXT=as.numeric(quantile(norm.EXT, na.rm=TRUE)[1]),
                                                            q10.norm.EXT=as.numeric(quantile(norm.EXT, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.norm.EXT=as.numeric(quantile(norm.EXT, probs=0.25, na.rm=TRUE)[1]),
                                                            median.norm.EXT=median(norm.EXT, na.rm=TRUE),
                                                            q75.norm.EXT=as.numeric(quantile(norm.EXT, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.norm.EXT=as.numeric(quantile(norm.EXT, probs=0.9, na.rm=TRUE)[1]),
                                                            max.norm.EXT=as.numeric(quantile(norm.EXT, na.rm=TRUE)[5]),
                                                            
                                                            mean.norm.red=mean(norm.red, na.rm=TRUE),
                                                            min.norm.red=as.numeric(quantile(norm.red, na.rm=TRUE)[1]),
                                                            q10.norm.red=as.numeric(quantile(norm.red, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.norm.red=as.numeric(quantile(norm.red, probs=0.25, na.rm=TRUE)[1]),
                                                            median.norm.red=median(norm.red, na.rm=TRUE),
                                                            q75.norm.red=as.numeric(quantile(norm.red, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.norm.red=as.numeric(quantile(norm.red, probs=0.9, na.rm=TRUE)[1]),
                                                            max.norm.red=as.numeric(quantile(norm.red, na.rm=TRUE)[5]),
                                                            
                                                            mean.norm.green=mean(norm.green, na.rm=TRUE),
                                                            median.norm.green=median(norm.green, na.rm=TRUE),
                                                            
                                                            mean.norm.yellow=mean(norm.yellow, na.rm=TRUE),
                                                            median.norm.yellow=median(norm.yellow, na.rm=TRUE),
                                                            
                                                            mean.log.EXT=mean(log(EXT)[!is.infinite(log(EXT))], na.rm=TRUE),

                                                            mean.log.red=mean(log(red)[!is.infinite(log(red))], na.rm=TRUE),
                                                            
                                                            mean.log.green=mean(log(green)[!is.infinite(log(green))], na.rm=TRUE),
                                                            
                                                            mean.log.yellow=mean(log(yellow)[!is.infinite(log(yellow))], na.rm=TRUE),
                                                            
                                                            mean.log.norm.EXT=mean(log(norm.EXT)[!is.infinite(log(norm.EXT))], na.rm=TRUE),
                                                            
                                                            mean.log.norm.red=mean(log(norm.red)[!is.infinite(log(norm.red))], na.rm=TRUE),
                                                            
                                                            mean.log.norm.green=mean(log(norm.green)[!is.infinite(log(norm.green))], na.rm=TRUE),
                                                            
                                                            mean.log.norm.yellow=mean(log(norm.yellow)[!is.infinite(log(norm.yellow))], na.rm=TRUE),
                                                            
                                                            f.L1 = length(which(stage == "L1"))/length(stage),
                                                            f.L2L3 = length(which(stage == "L2/L3"))/length(stage),
                                                            f.L4 = length(which(stage == "L4"))/length(stage),
                                                            f.ad = length(which(stage == "adult"))/length(stage))
    
    if(!ends){
        processed <- processed[,-(grep("min", colnames(processed)))]
        processed <- processed[,-(grep("max", colnames(processed)))]
    }
    if(!quantiles){
        processed <- processed[,-(grep("q", colnames(processed)))]
    }
    if(!log){
        processed <- processed[,-(grep("log", colnames(processed)))]
    }
    if(is.null(strains)){
        analysis <- processed
        analysis <- analysis[order(analysis$row, analysis$col),]
    } else {
        analysis <- data.frame(strain = as.character(strains), processed)
        analysis <- analysis[order(analysis$strain),]
        analysis <- analysis[order(analysis$row, analysis$col),]
    }
    analysis[analysis$mean.TOF==-1 | is.na(analysis$mean.TOF),which(colnames(analysis)=="n"):ncol(analysis)] <- NA
    return(analysis)
}


readPlate_worms <- function(file, tofmin=60, tofmax=2000, extmin=0, extmax=10000, SVM=TRUE) {
    plate <- readSorter(file, tofmin, tofmax, extmin, extmax)
    modplate <- with(plate, data.frame(row=Row, col=as.factor(Column), sort=Status.sort, TOF=TOF, EXT=EXT, time=Time.Stamp, green=Green, yellow=Yellow, red=Red))
    modplate <- modplate %>% group_by(row, col) %>% do(extractTime(.))
    modplate[,10:13] <- apply(modplate[,c(5, 7:9)], 2, function(x){x/modplate$TOF})
    colnames(modplate)[10:13] <- c("norm.EXT", "norm.green", "norm.yellow", "norm.red")
    if(SVM){
        plateprediction <- predict(bubbleSVMmodel_noProfiler, modplate[,3:length(modplate)], type="probabilities")
        modplate$object <- plateprediction[,"1"]
        modplate$call50 <- factor(as.numeric(modplate$object>0.5), levels=c(1,0), labels=c("object", "bubble"))
    }
    modplate$stage <- ifelse(modplate$TOF>=60 & modplate$TOF<90, "L1", 
                             ifelse(modplate$TOF>=90 & modplate$TOF<200, "L2/L3",
                                    ifelse(modplate$TOF>=200 & modplate$TOF<300, "L4",
                                           ifelse(modplate$TOF>=300, "adult", NA))))
    return(modplate)
}

removeWells <- function(plate, badWells, drop=FALSE) {
    sp.bw <- str_split(badWells, "", 3)
    if(!drop){
        for (i in seq(1, length(sp.bw))) {
            row <- sp.bw[[i]][2]
            col <- sp.bw[[i]][3]
            plate[which(plate$row == row & plate$col == col),10:ncol(plate)] <- NA
        }
    } else {
        for (i in seq(1, length(sp.bw))) {
            row <- sp.bw[[i]][2]
            col <- sp.bw[[i]][3]
            plate = plate[plate$row != row & plate$col != col,]
        }
    }
    return(plate)
}
