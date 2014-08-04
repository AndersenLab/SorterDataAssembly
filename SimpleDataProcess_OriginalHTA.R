library(COPASutils)
library(plyr)
library(dplyr)
library(kernlab)

source("~/SorterDataAssembly/SimpleDataProcessFxns.R")

# data=completeData[completeData$drug=="paraquat",]

regress <- function(data, completeData, controls){
    plates <- data[!duplicated(data[,c("assay", "plate", "drug")]), c("assay", "plate", "drug")]
    controlValues <- plates %>%
        group_by(assay, plate, drug) %>%
        do(try({data.frame(filter(completeData,
                                  assay==as.character(.$assay[1]),
                                  as.numeric(as.character(plate)) %in% eval(eval(controls[controls$assay==.$assay[1],"control"]))))}))
    controlValues <- controlValues %>% group_by(strain) %>%
        summarise_each(funs(mean(., na.rm=TRUE)), -date, -experiment, -round, -assay, -plate, -drug, -row, -col)
    
    finalControl = data.frame(matrix(ncol=ncol(controlValues), nrow=0))
    colnames(finalControl)=colnames(controlValues)
    
    for(i in 1:nrow(data)){
        if(length(which(controlValues$strain==data$strain[i]))==0){
            finalControl = rbind.fill(finalControl, data.frame(NA))
        }else{
            finalControl = rbind.fill(finalControl, controlValues[which(controlValues$strain==data$strain[i]),])
        }
    }
    
    regressedValues <- data.frame(do.call(cbind, lapply(which(colnames(data)=="n"):ncol(data),
                                                        function(x){
                                                            tryCatch({residuals(lm(data[,x] ~ data$assay + finalControl[,which(colnames(finalControl)==colnames(data)[x])], na.action=na.exclude))},
                                                                     error = function(err){print(err);return(NA)})
                                                        })))
    
    
    reactValues <- data.frame(do.call(cbind, lapply(which(colnames(data)=="n"):ncol(data),
                                                    function(x){
                                                        reactNorms <- data[,x] - finalControl[,which(colnames(finalControl)==colnames(data)[x])]
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










generateReports=FALSE

options(echo=TRUE)

generateReports <- FALSE

# Set all experiment directories
directories <- c("~/Dropbox/HTA/Results/20110811_RIAILs0a", "~/Dropbox/HTA/Results/20110812_RIAILs0b", "~/Dropbox/HTA/Results/20110818_RIAILs0c", "~/Dropbox/HTA/Results/20110819_RIAILs0d")

sapply(directories, function(x){dir.create(file.path(x, "reports"), showWarnings = FALSE)})

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

if(generateReports){
    sortData %>% group_by(assay, plate) %>% do(data.frame(setupReport(.)))
}

# Read in all score plates
rawScoreData <- data.frame(do.call(rbind, lapply(scoreList, function(x){data.frame(cbind(info(x,2), readPlate_worms(x, 60, 2000), step="score"))}))) %>%
    filter(drug != "missing")

# Get the strain data
strainFiles <- do.call(rbind, lapply(lapply(directories,
                                            function(x){
                                                lapply(list.files(x,
                                                                  pattern="\\.Rds",
                                                                  recursive=TRUE,
                                                                  full.names=TRUE),
                                                       function(y){
                                                           if(!grepl("IncompleteData",y) & !grepl("UnstitchedData",y) & !grepl("Archive",y) & grepl("strains",y)){
                                                               return(y)
                                                           }
                                                       })
                                            }), function(x){
                                                y=unlist(x)
                                                data.frame(cbind(y[1], y[2]))
                                            }))

strainsData <- do.call(rbind, apply(strainFiles, 1, function(x){
    load(x[1])
    load(x[2])
    ctrls <- data.frame(assay=info(x[1])$assay, plate=1:length(ctrl.strains), strains=I(ctrl.strains))
    pq <- data.frame(assay=info(x[2])$assay, plate=((length(ctrl.strains)+1):(length(ctrl.strains)+length(pq.strains))), strains=I(pq.strains))
    return(data.frame(rbind(ctrls, pq)))
}))
colnames(strainsData) <- c("assay", "plate", "strains")

# Summarize the score data
summarizedScoreData <- rawScoreData %>%
    group_by(date, experiment, round, assay, plate, drug) %>%
    do({summarizePlate_worms(.,
                            strains=eval(strainsData[as.character(strainsData$assay)==as.character(.$assay[1]) & as.numeric(strainsData$plate)==as.numeric(as.character(.$plate[1])),
                                                                                                 "strains"])[[1]], quantiles=TRUE)})

# Join the number sorted and calculate norm.n, then remove the number sorted 
completeData <- left_join(summarizedScoreData, select(sortData, assay, plate, row, col, n.sorted.setup = n.sorted)) %>% mutate(norm.n=n/n.sorted.setup) %>% select(-contains("sort"))

#NA out wash wells
completeData[is.na(completeData$strain) | is.na(completeData$norm.n), which(colnames(completeData)=="n"):ncol(completeData)] <- NA

controls <- data.frame()

for(i in unique(completeData$assay)){
    data <- filter(completeData, assay==i)
    controlPlates <- unique(as.numeric(as.character(filter(data, drug=="control")$plate)))
    pqPlates <- unique(as.numeric(as.character(filter(data, drug=="paraquat")$plate)))
    row <- data.frame(matrix(nrow=1, ncol=3))
    row[[1,1]] <- i
    row[[1,2]] <- parse(text=as.character(call("c", controlPlates))[2])
    row[[1,3]] <- parse(text=as.character(call("c", pqPlates))[2])
    controls <- rbind(controls, row)
}

colnames(controls) <- c("assay", "control", "plates")

finalData <- completeData %>% group_by(drug) %>% do(regress(., completeData, controls)) %>% arrange(assay)

finalData[finalData$drug=="control",which(colnames(finalData)=="resid.n"):ncol(finalData)] <- NA

experiment <- info(directories[1], levels=0)$experiment
round <- info(directories[1], levels=0)$round

write.csv(finalData, paste0("~/Dropbox/HTA/Results/ProcessedData/", experiment, round, "_complete_simple.csv"))
