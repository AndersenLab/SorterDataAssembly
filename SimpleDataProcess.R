require(COPASutils)

# Extract the metadata info
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
    plate <- readSorter(file, tofmin, tofmax, extmin)
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










# Set all experiment directories
directories <- c("~/Dropbox/HTA/Results/20110811_RIAILs0a", "~/Dropbox/HTA/Results/20110812_RIAILs0b", "~/Dropbox/HTA/Results/20110818_RIAILs0c", "~/Dropbox/HTA/Results/20110819_RIAILs0d")

directories <- c("~/Dropbox/HTA/Results/20140616_GWAS8a", "~/Dropbox/HTA/Results/20140617_GWAS8b")

directories <- c("~/Dropbox/HTA/Results/20140609_GWAS7a", "~/Dropbox/HTA/Results/20140610_GWAS7b")

directories <- c("~/Dropbox/HTA/Results/20140421_GWAS6a", "~/Dropbox/HTA/Results/20140422_GWAS6b")

directories <- c("~/Dropbox/HTA/Results/20140414_GWAS5a", "~/Dropbox/HTA/Results/20140415_GWAS5b")

directories <- c("~/Dropbox/HTA/Results/20140407_GWAS4a", "~/Dropbox/HTA/Results/20140408_GWAS4b")

directories <- c("~/Dropbox/HTA/Results/20140331_GWAS3a", "~/Dropbox/HTA/Results/20140401_GWAS3b")

directories <- c("~/Dropbox/HTA/Results/20140324_GWAS2a", "~/Dropbox/HTA/Results/20140325_GWAS2b")

directories <- c("~/Dropbox/HTA/Results/20140317_GWAS1a", "~/Dropbox/HTA/Results/20140318_GWAS1b")


# Get the setup list vector
setupList <- sapply(unlist(sapply(directories,
                    function(x){lapply(list.files(x, pattern="\\.txt", recursive=TRUE, full.names=TRUE),
                                       function(y){if(!grepl("IncompleteData",y) & !grepl("UnstitchedData",y) & !grepl("Archive",y) & grepl("setup",y)){return(y)}})})
                          ), as.character)

# Get the score list vector
scoreList <- sapply(unlist(sapply(directories,
                   function(x){lapply(list.files(x, pattern="\\.txt", recursive=TRUE, full.names=TRUE),
                                      function(y){if(!grepl("IncompleteData",y) & !grepl("UnstitchedData",y) & !grepl("Archive",y) & grepl("score",y)){return(y)}})})
                          ), as.character)

# Read in all setup plates
sortData <- do.call(rbind, lapply(setupList, function(x){data.frame(cbind(info(x,2), procSetup(x, 60, 2000), step="setup"))}))

# Read in all score plates
rawScoreData <- do.call(rbind, lapply(scoreList, function(x){data.frame(cbind(info(x,2), readPlate(x, 60, 2000), step="score"))}))

# Get the strain data
strainFiles <- sapply(unlist(sapply(directories,
                                    function(x){
                                        lapply(list.files(x,
                                                          pattern="\\.R",
                                                          recursive=TRUE,
                                                          full.names=TRUE),
                                               function(y){
                                                   if(!grepl("IncompleteData",y) & !grepl("UnstitchedData",y) & !grepl("Archive",y) & grepl("strains",y)){
                                                       return(y)
                                                   }
                                               })
                                    })), as.character)

strainsData <- do.call(rbind, lapply(strainFiles, function(x){
                                                    source(x)
                                                    data.frame(cbind(assay=as.character(info(x)$assay),
                                                        I(list(strains))))
                                                  }))
colnames(strainsData) <- c("assay", "strains")

# Summarize the score data
summarizedScoreData <- rawScoreData %>% group_by(date, experiment, round, assay, plate, drug) %>% do(summarizePlate(., strains=eval(strainsData[strainsData$assay==.$assay[1],2])[[1]], quantiles=TRUE))

# Join the number sorted and calculate norm.n, then remove the number sorted 
completeData <- left_join(summarizedScoreData, select(sortData, assay, plate, row, col, n.sorted.setup = n.sorted)) %>% mutate(norm.n=n/n.sorted.setup) %>% select(-contains("sort"))

#Read in and handle contamination
contamFiles <- sapply(unlist(sapply(directories,
                                    function(x){
                                        lapply(list.files(x,
                                                          pattern="\\.R",
                                                          recursive=TRUE,
                                                          full.names=TRUE),
                                               function(y){
                                                   if(!grepl("IncompleteData",y) & !grepl("UnstitchedData",y) & !grepl("Archive",y) & grepl("contamination",y)){
                                                       return(y)
                                                    }
                                                })
                                        })), as.character)

contamination <- do.call(rbind, lapply(contamFiles, function(x){data.frame(cbind(assay=info(x)$assay,
                                                      read.delim(x, sep="<", quote=NULL, header=FALSE)))}))
colnames(contamination) <- c("assay", "plate", "contam")

contamination$plate <- sapply(contamination$plate, function(x){as.numeric(strsplit(as.character(x), "p")[[1]][2])})
contamination$contam <- I(sapply(contamination$contam, function(x){eval(parse(text=gsub("- ", "", as.character(x))))}))

#Remove contamination
completeData <- completeData %>% group_by(assay, plate) %>% do(removeWells(., unlist(contamination[contamination$assay==.$assay[1] & contamination$plate==as.numeric(.$plate[1]),3])))

#NA out wash wells, change if patterning is different than for mapping
completeData[is.na(completeData$strain), which(colnames(completeData)=="n"):ncol(completeData)] <- NA

#NA out any wells where norm.n is infinite (i.e. nothing sorted to well)
completeData[is.infinite(completeData$norm.n), which(colnames(completeData)=="n"):ncol(completeData)] <- NA

controls <- do.call(rbind, lapply(lapply(directories, function(x){list.files(x, pattern="controls", recursive=TRUE, full.names=TRUE)}),
                                  function(x){
                                      source(x)
                                      assay <- info(x)$assay
                                      data.frame(assay, control=I(controlPlates), plates=I(testPlates))
                                  }))

finalData <- completeData %>% group_by(drug) %>% do(regress(., completeData, controls))# %>% arrange(assay)





for(i in unique(completeData$drug)){
    print(i)
    data = completeData[completeData$drug==i,]
    print(nrow(regress(data, completeData, controls)))
}







write.csv(finalData, "~/Dropbox/HTA/Results/ProcessedData/GWAS7_complete_simple.csv")

write.csv(finalData, "~/Dropbox/HTA/Results/ProcessedData/GWAS6_complete_simple.csv")

write.csv(finalData, "~/Dropbox/HTA/Results/ProcessedData/GWAS5_complete_simple.csv")

write.csv(finalData, "~/Dropbox/HTA/Results/ProcessedData/GWAS4_complete_simple.csv")

write.csv(finalData, "~/Dropbox/HTA/Results/ProcessedData/GWAS3_complete_simple.csv")

write.csv(finalData, "~/Dropbox/HTA/Results/ProcessedData/GWAS2_complete_simple.csv")

write.csv(finalData, "~/Dropbox/HTA/Results/ProcessedData/GWAS1_complete_simple.csv")
