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

procSetup <- function(file, tofmin=0, tofmax=2000, extmin=20, extmax=5000) {
    
    #Read in the sorter data from the file
    plate <- readSorter(file, tofmin, tofmax, extmin)
    modplate <- with(plate, data.frame(row=Row, col=as.factor(Column),
                                       sort = Status.sort, TOF=TOF, EXT=EXT,
                                       time=Time.Stamp, green=Green, yellow=Yellow,
                                       red=Red))
    
    #Process the data
    proc <- modplate %>% group_by(row, col) %>% summarise(pop = length(EXT), sort = sum(sort==6), TOF = mean(TOF), EXT = mean(EXT),TOFmed=median(TOF),EXTmed=median(EXT))
    
    return(proc)
}


# Set all experiment directories
directories <- c("~/Dropbox/HTA/Results/20110811_RIAILs0a", "~/Dropbox/HTA/Results/20110812_RIAILs0b", "~/Dropbox/HTA/Results/20110818_RIAILs0c", "~/Dropbox/HTA/Results/20110819_RIAILs0d")

# Get the setup list vector
setupList <- sapply(directories,
                    function(x){lapply(list.files(x, pattern="\\.txt", recursive=TRUE, full.names=TRUE),
                                       function(y){if(!grepl("IncompleteData",y) & !grepl("UnstitchedData",y) & !grepl("Archive",y) & grepl("setup",y)){return(y)}})})
setupList <- unlist(setupList)
setupList <- sapply(setupList, as.character)

# Get the score list vector
scoreList <- sapply(directories,
                   function(x){lapply(list.files(x, pattern="\\.txt", recursive=TRUE, full.names=TRUE),
                                      function(y){if(!grepl("IncompleteData",y) & !grepl("UnstitchedData",y) & !grepl("Archive",y) & grepl("score",y)){return(y)}})})
scoreList <- unlist(scoreList)
scoreList <- sapply(scoreList, as.character)

# Read in all setup plates
sortData <- do.call(rbind, lapply(setupList, function(x){data.frame(cbind(info(x,2), procSetup(x, 60, 2000), step="setup"))}))

# Read in all score plates
rawScoreData <- do.call(rbind, lapply(scoreList, function(x){data.frame(cbind(info(x,2), readPlate_worms(x, 60, 2000), step="score"))}))
summarizedScoreData <- rawData %>% group_by(assay, plate) %>% do(summarizePlate(.))





















# raw3 <- do.call(rbind, lapply(fileList, function(x){data.frame(cbind(info(x,2), readPlate_worms(x, 60, 2000)))}))
# raw2 <- do.call(rbind, lapply(fileList, function(x){data.frame(cbind(info(x,2), readPlate(x, 60, 2000)))}))
# 
# sum3 <- raw3 %>% group_by(assay, plate) %>% do(summarizePlate(.))
# sum2 <- raw2 %>% group_by(assay, plate) %>% do(summarizePlate(.))
# 
# 
# write.csv(summarizedData, "RIAILs0_raw2.csv")


