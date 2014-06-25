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
    proc <- modplate %>% group_by(row, col) %>% summarise(n = length(EXT), n.sorted = sum(sort==6), mean.TOF = mean(TOF), mean.EXT = mean(EXT), median.TOF=median(TOF), median.EXT=median(EXT))
    
    return(proc)
}

# Set all experiment directories
directories <- c("~/Dropbox/HTA/Results/20110811_RIAILs0a", "~/Dropbox/HTA/Results/20110812_RIAILs0b", "~/Dropbox/HTA/Results/20110818_RIAILs0c", "~/Dropbox/HTA/Results/20110819_RIAILs0d")

directories <- c("~/Dropbox/HTA/Results/20140616_GWAS8a", "~/Dropbox/HTA/Results/20140617_GWAS8b")

directories <- c("~/Dropbox/HTA/Results/20140609_GWAS7a", "~/Dropbox/HTA/Results/20140610_GWAS7b")



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
rawScoreData <- do.call(rbind, lapply(scoreList, function(x){data.frame(cbind(info(x,2), readPlate(x, 60, 2000), step="score"))}))
summarizedScoreData <- rawScoreData %>% group_by(date, experiment, round, assay, plate, drug) %>% do(summarizePlate(., quantiles=TRUE))

completeData <- left_join(summarizedScoreData, select(sortData, assay, plate, row, col, n.sorted.setup = n.sorted)) %>% mutate(norm.n=n/n.sorted.setup) %>% select(-contains("sort"))

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
                                                      read.delim(x, sep="<", quote="\'", header=FALSE)))}))
colnames(contamination) <- c("assay", "plate", "contam")

contamination$plate <- sapply(contamination$plate, function(x){as.numeric(strsplit(as.character(x), "p")[[1]][2])})
contamination$contam <- I(sapply(contamination$contam, function(x){eval(parse(text=gsub("- ", "", as.character(x))))}))

completeData <- completeData %>% group_by(assay, plate) %>% do(removeWells(., unlist(contamination[contamination$assay==.$assay[1] & contamination$plate==as.numeric(.$plate[1]),3])))

#NA out wash wells, change if patterning is different than for mapping
completeData[as.numeric(completeData$col) %% 2 == 0, which(colnames(completeData)=="n"):ncol(completeData)] <- NA

#NA out any wells where norm.n is infinite (i.e. nothing sorted to well)
completeData[is.infinite(completeData$norm.n), which(colnames(completeData)=="n"):ncol(completeData)] <- NA














regress <- function(data, completeData){
    plates <- data[!duplicated(data[,c("assay", "plate")]), c("assay", "plate")]
    do.call(rbind, apply(plates, 1, function(x){completeData %>% filter(assay==x[["assay"]], plate %in% controls[controls$assay==assay & x[["plate"]] %in% controls$plate, "control"])}))
}


controls <- do.call(rbind, lapply(lapply(directories, function(x){list.files(x, pattern="controls", recursive=TRUE, full.names=TRUE)}),
                   function(x){
                       source(x)
                       assay <- info(x)$assay
                       data.frame(assay, control=I(controlPlates), plates=I(testPlates))
                   }))

completeData %>% group_by(drug)

write.csv(completeData, "~/Dropbox/HTA/Results/ProcessedData/GWAS8_complete_simple.csv")













# raw3 <- do.call(rbind, lapply(fileList, function(x){data.frame(cbind(info(x,2), readPlate_worms(x, 60, 2000)))}))
# raw2 <- do.call(rbind, lapply(fileList, function(x){data.frame(cbind(info(x,2), readPlate(x, 60, 2000)))}))
# 
# sum3 <- raw3 %>% group_by(assay, plate) %>% do(summarizePlate(.))
# sum2 <- raw2 %>% group_by(assay, plate) %>% do(summarizePlate(.))
# 
# 
# write.csv(summarizedData, "RIAILs0_raw2.csv")


