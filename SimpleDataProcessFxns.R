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

#Process the setup plates
procSetup <- function(file, tofmin=0, tofmax=2000, extmin=20, extmax=5000) {
    
    #Read in the sorter data from the file
    plate <- readSorter(file, tofmin, tofmax, extmin, extmax)
    modplate <- with(plate, data.frame(row=Row, col=as.factor(Column),
                                       sort = Status.sort, TOF=TOF, EXT=EXT,
                                       time=Time.Stamp, green=Green, yellow=Yellow,
                                       red=Red))
    
    #Process the data
    proc <- modplate %>% group_by(row, col) %>% summarise(n = length(EXT), n.sorted = sum(sort==6), mean.TOF = mean(TOF), mean.EXT = mean(EXT), median.TOF=median(TOF), median.EXT=median(EXT))
    
    return(proc)
}

regress <- function(data, completeData, controls){
    plates <- data[!duplicated(data[,c("assay", "plate", "drug")]), c("assay", "plate", "drug")]
    controlValues <- plates %>%
        group_by(assay, plate, drug) %>%
        do(tryCatch({data.frame(filter(completeData,
                                       assay==as.character(.$assay[1]),
                                       as.numeric(plate) %in% as.numeric(unlist(controls[sapply(controls$plates,
                                                                                                function(x){as.numeric(.$plate[1]) %in% as.numeric(x)}) & controls$assay==.$assay[1], "control"])))) %>%
                         group_by(row, col) %>%
                         summarise_each(funs(mean(., na.rm=TRUE)), -date, -experiment, -round, -assay, -plate, -drug)},
                    error = function(err){return(data.frame(matrix(nrow=96)))}))
    
    regressedValues <- data.frame(do.call(cbind, lapply(which(colnames(data)=="n"):ncol(data),
                                                        function(x){
                                                            tryCatch({residuals(lm(data[,x] ~ data$assay + controlValues[,which(colnames(controlValues)==colnames(data)[x])], na.action=na.exclude))},
                                                                     error = function(err){return(NA)})
                                                        })))
    
    
    reactValues <- data.frame(do.call(cbind, lapply(which(colnames(data)=="n"):ncol(data),
                                                    function(x){
                                                        reactNorms <- data[,x] - controlValues[,which(colnames(controlValues)==colnames(data)[x])]
                                                        if(length(reactNorms)==0){
                                                            reactNorms <- NA
                                                        }
                                                        return(reactNorms)
                                                    })))
    
    finalDF <- data.frame(data, regressedValues)
    colnames(finalDF)[(which(colnames(finalDF)=="norm.n")+1):ncol(finalDF)] <- paste0("resid.", colnames(finalDF)[which(colnames(finalDF)=="n"):which(colnames(finalDF)=="norm.n")])
    finalDF <- data.frame(finalDF, reactValues)
    colnames(finalDF)[(which(colnames(finalDF)=="resid.norm.n")+1):ncol(finalDF)] <- paste0("react.", colnames(finalDF)[which(colnames(finalDF)=="n"):which(colnames(finalDF)=="norm.n")])
    return(finalDF)
}


possContam <- function(procDataFrame){
    strainMean <- mean(procDataFrame$n[!is.na(procDataFrame$strain)], na.rm = TRUE)
    strainSD <- sd(procDataFrame$n[!is.na(procDataFrame$strain)], na.rm = TRUE)
    possibleContam <- c()
    for(j in seq(1,nrow(procDataFrame),)){
        if(!is.na(as.character(procDataFrame[j,"strain"])) & !is.na(procDataFrame[j,"n"])){
            if(procDataFrame[j,"n"] > strainMean + (3*strainSD)){
                row <- as.character(procDataFrame[j, "row"])
                col <- as.numeric(procDataFrame[j, "col"])
                possibleContam <- append(possibleContam, paste0(row, col))
            }
        }
    }
    return(possibleContam)
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

summarizePlate_worms <- function(plate, strains=NULL, quantiles=FALSE, log=FALSE, ends=FALSE) {
    plate <- plate[plate$call50=="object" | plate$TOF == -1 | is.na(plate$call50),]
    plate <- fillWells(plate)
    processed <- plate %>% group_by(row, col) %>% summarise(n=ifelse(length(TOF[!is.na(TOF)])==0, NA, length(TOF[!is.na(TOF)])),
                                                            n.sorted=sum(sort==6),
                                                            
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
                                                            min.green=as.numeric(quantile(green, na.rm=TRUE)[1]),
                                                            q10.green=as.numeric(quantile(green, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.green=as.numeric(quantile(green, probs=0.25, na.rm=TRUE)[1]),
                                                            median.green=median(green, na.rm=TRUE),
                                                            q75.green=as.numeric(quantile(green, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.green=as.numeric(quantile(green, probs=0.9, na.rm=TRUE)[1]),
                                                            max.green=as.numeric(quantile(green, na.rm=TRUE)[5]),
                                                            
                                                            mean.yellow=mean(yellow, na.rm=TRUE),
                                                            min.yellow=as.numeric(quantile(yellow, na.rm=TRUE)[1]),
                                                            q10.yellow=as.numeric(quantile(yellow, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.yellow=as.numeric(quantile(yellow, probs=0.25, na.rm=TRUE)[1]),
                                                            median.yellow=median(yellow, na.rm=TRUE),
                                                            q75.yellow=as.numeric(quantile(yellow, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.yellow=as.numeric(quantile(yellow, probs=0.9, na.rm=TRUE)[1]),
                                                            max.yellow=as.numeric(quantile(yellow, na.rm=TRUE)[5]),
                                                            
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
                                                            min.norm.green=as.numeric(quantile(norm.green, na.rm=TRUE)[1]),
                                                            q10.norm.green=as.numeric(quantile(norm.green, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.norm.green=as.numeric(quantile(norm.green, probs=0.25, na.rm=TRUE)[1]),
                                                            median.norm.green=median(norm.green, na.rm=TRUE),
                                                            q75.norm.green=as.numeric(quantile(norm.green, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.norm.green=as.numeric(quantile(norm.green, probs=0.9, na.rm=TRUE)[1]),
                                                            max.norm.green=as.numeric(quantile(norm.green, na.rm=TRUE)[5]),
                                                            
                                                            mean.norm.yellow=mean(norm.yellow, na.rm=TRUE),
                                                            min.norm.yellow=as.numeric(quantile(norm.yellow, na.rm=TRUE)[1]),
                                                            q10.norm.yellow=as.numeric(quantile(norm.yellow, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.norm.yellow=as.numeric(quantile(norm.yellow, probs=0.25, na.rm=TRUE)[1]),
                                                            median.norm.yellow=median(norm.yellow, na.rm=TRUE),
                                                            q75.norm.yellow=as.numeric(quantile(norm.yellow, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.norm.yellow=as.numeric(quantile(norm.yellow, probs=0.9, na.rm=TRUE)[1]),
                                                            max.norm.yellow=as.numeric(quantile(norm.yellow, na.rm=TRUE)[5]),
                                                            
                                                            mean.log.EXT=mean(log(EXT), na.rm=TRUE),
                                                            min.log.EXT=as.numeric(quantile(log(EXT), na.rm=TRUE)[1]),
                                                            q10.log.EXT=as.numeric(quantile(log(EXT), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.EXT=as.numeric(quantile(log(EXT), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.EXT=median(log(EXT), na.rm=TRUE),
                                                            q75.log.EXT=as.numeric(quantile(log(EXT), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.EXT=as.numeric(quantile(log(EXT), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.EXT=as.numeric(quantile(log(EXT), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.red=mean(log(red), na.rm=TRUE),
                                                            min.log.red=as.numeric(quantile(log(red), na.rm=TRUE)[1]),
                                                            q10.log.red=as.numeric(quantile(log(red), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.red=as.numeric(quantile(log(red), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.red=median(log(red), na.rm=TRUE),
                                                            q75.log.red=as.numeric(quantile(log(red), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.red=as.numeric(quantile(log(red), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.red=as.numeric(quantile(log(red), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.green=mean(log(green), na.rm=TRUE),
                                                            min.log.green=as.numeric(quantile(log(green), na.rm=TRUE)[1]),
                                                            q10.log.green=as.numeric(quantile(log(green), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.green=as.numeric(quantile(log(green), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.green=median(log(red), na.rm=TRUE),
                                                            q75.log.green=as.numeric(quantile(log(green), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.green=as.numeric(quantile(log(green), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.green=as.numeric(quantile(log(green), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.yellow=mean(log(yellow), na.rm=TRUE),
                                                            min.log.yellow=as.numeric(quantile(log(yellow), na.rm=TRUE)[1]),
                                                            q10.log.yellow=as.numeric(quantile(log(yellow), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.yellow=as.numeric(quantile(log(yellow), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.yellow=median(log(yellow), na.rm=TRUE),
                                                            q75.log.yellow=as.numeric(quantile(log(yellow), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.yellow=as.numeric(quantile(log(yellow), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.yellow=as.numeric(quantile(log(yellow), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.norm.EXT=mean(log(norm.EXT), na.rm=TRUE),
                                                            min.log.norm.EXT=as.numeric(quantile(log(norm.EXT), na.rm=TRUE)[1]),
                                                            q10.log.norm.EXT=as.numeric(quantile(log(norm.EXT), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.norm.EXT=as.numeric(quantile(log(norm.EXT), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.norm.EXT=median(log(norm.EXT), na.rm=TRUE),
                                                            q75.log.norm.EXT=as.numeric(quantile(log(norm.EXT), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.norm.EXT=as.numeric(quantile(log(norm.EXT), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.norm.EXT=as.numeric(quantile(log(norm.EXT), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.norm.red=mean(log(norm.red), na.rm=TRUE),
                                                            min.log.norm.red=as.numeric(quantile(log(norm.red), na.rm=TRUE)[1]),
                                                            q10.log.norm.red=as.numeric(quantile(log(norm.red), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.norm.red=as.numeric(quantile(log(norm.red), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.norm.red=median(log(norm.red), na.rm=TRUE),
                                                            q75.log.norm.red=as.numeric(quantile(log(norm.red), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.norm.red=as.numeric(quantile(log(norm.red), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.norm.red=as.numeric(quantile(log(norm.red), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.norm.green=mean(log(norm.green), na.rm=TRUE),
                                                            min.log.norm.green=as.numeric(quantile(log(norm.green), na.rm=TRUE)[1]),
                                                            q10.log.norm.green=as.numeric(quantile(log(norm.green), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.norm.green=as.numeric(quantile(log(norm.green), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.norm.green=median(log(norm.red), na.rm=TRUE),
                                                            q75.log.norm.green=as.numeric(quantile(log(norm.green), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.norm.green=as.numeric(quantile(log(norm.green), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.norm.green=as.numeric(quantile(log(norm.green), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.norm.yellow=mean(log(norm.yellow), na.rm=TRUE),
                                                            min.log.norm.yellow=as.numeric(quantile(log(norm.yellow), na.rm=TRUE)[1]),
                                                            q10.log.norm.yellow=as.numeric(quantile(log(norm.yellow), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.norm.yellow=as.numeric(quantile(log(norm.yellow), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.norm.yellow=median(log(norm.yellow), na.rm=TRUE),
                                                            q75.log.norm.yellow=as.numeric(quantile(log(norm.yellow), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.norm.yellow=as.numeric(quantile(log(norm.yellow), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.norm.yellow=as.numeric(quantile(log(norm.yellow), na.rm=TRUE)[5]),
                                                            
                                                            var.TOF=var(TOF),
                                                            iqr.TOF=quantile(TOF, na.rm=TRUE, probs=.75)-quantile(TOF, na.rm=TRUE, probs=.25),
                                                            var.EXT=var(EXT),
                                                            iqr.EXT=quantile(EXT, na.rm=TRUE, probs=.75)-quantile(EXT, na.rm=TRUE, probs=.25),
                                                            var.red=var(red),
                                                            iqr.red=quantile(red, na.rm=TRUE, probs=.75)-quantile(red, na.rm=TRUE, probs=.25),
                                                            var.green=var(green),
                                                            iqr.green=quantile(green, na.rm=TRUE, probs=.75)-quantile(green, na.rm=TRUE, probs=.25),
                                                            var.yellow=var(yellow),
                                                            iqr.yellow=quantile(yellow, na.rm=TRUE, probs=.75)-quantile(yellow, na.rm=TRUE, probs=.25),
                                                            
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

setupReport <- function(df){
    for(row in 1:nrow(df)){
        if(df$n[row] < 15) {
            df$bin[row] = "low"
        } else if(df$n[row] < 25) {
            df$bin[row] = "low-mid"
        } else if(df$n[row] < 50) {
            df$bin[row] = "high-mid"
        } else {
            df$bin[row] = "high"
        }
    }
    
    filename <- paste0("p", df$plate[1], "_", df$drug[1], "_setup.html")
    
    knit2html("~/SorterDataAssembly/MasterSetupReport2.Rmd", filename)
    
    file.rename(filename, file.path("~/Dropbox/HTA/Results", paste0(df$date[1], "_", df$experiment[1], df$round[1], df$assay[1]), "reports", filename))
    return(data.frame(done = 1))
}

meltdf <- function(score){
    newscore<-data.frame(row=rep(score$row,each=1),col=rep(score$col,each=1),n=rep(score$n,each=1),f.L1=rep(score$f.L1,each=1),f.L2L3=rep(score$f.L2L3,each=1),f.L4=rep(score$f.L4,each=1),f.ad=rep(score$f.ad,each=1))
    newscore<-melt(newscore,id.var=c("row","col","n"))
    return(newscore)
}

scoreReport <- function(df, contamination){
    contamination <- filter(contamination, as.character(assay) == as.character(df$assay[1]), as.numeric(as.character(plate)) == as.numeric(as.character(df$plate[1])))
    melted.proc <- meltdf(df)
    filename <- paste0("p", df$plate[1], "_", df$drug[1], "_score.html")
    knit2html("~/SorterDataAssembly/MasterScoreReport2.Rmd", filename)
    file.rename(filename, file.path("~/Dropbox/HTA/Results", paste0(df$date[1], "_", df$experiment[1], df$round[1], df$assay[1]), "reports", filename))
    return(data.frame(done = 1))
}