library(COPASutils)
library(dplyr)
library(kernlab)
library(knitr)

source("~/SorterDataAssembly/SimpleDataProcessFxns.R")
generateReports=FALSE

options(echo=TRUE)

args <- commandArgs(trailingOnly = TRUE)
generateReports <- as.logical(args[1])
directories <- args[2:length(args)]

# Set all experiment directories
#directories <- c("~/Dropbox/HTA/Results/20110811_RIAILs0a", "~/Dropbox/HTA/Results/20110812_RIAILs0b", "~/Dropbox/HTA/Results/20110818_RIAILs0c", "~/Dropbox/HTA/Results/20110819_RIAILs0d")
# 
# directories <- c("~/Dropbox/HTA/Results/20140616_GWAS8a", "~/Dropbox/HTA/Results/20140617_GWAS8b")
# 
# directories <- c("~/Dropbox/HTA/Results/20140609_GWAS7a", "~/Dropbox/HTA/Results/20140610_GWAS7b")
# 
# directories <- c("~/Dropbox/HTA/Results/20140421_GWAS6a", "~/Dropbox/HTA/Results/20140422_GWAS6b")
# 
# directories <- c("~/Dropbox/HTA/Results/20140414_GWAS5a", "~/Dropbox/HTA/Results/20140415_GWAS5b")
# 
# directories <- c("~/Dropbox/HTA/Results/20140407_GWAS4a", "~/Dropbox/HTA/Results/20140408_GWAS4b")
# 
# directories <- c("~/Dropbox/HTA/Results/20140331_GWAS3a", "~/Dropbox/HTA/Results/20140401_GWAS3b")
# 
# directories <- c("~/Dropbox/HTA/Results/20140324_GWAS2a", "~/Dropbox/HTA/Results/20140325_GWAS2b")
# 
# directories <- c("~/Dropbox/HTA/Results/20140317_GWAS1a", "~/Dropbox/HTA/Results/20140318_GWAS1b")

# directories <- c("~/Dropbox/HTA/Results/20140429_RIAILs1b") 

# directories=c("/Users/tyler/Dropbox/HTA/Results/20140630_RIAILs2a", "/Users/tyler/Dropbox/HTA/Results/20140701_RIAILs2b", "/Users/tyler/Dropbox/HTA/Results/20140707_RIAILs2c", "/Users/tyler/Dropbox/HTA/Results/20140708_RIAILs2d")

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
rawScoreData <- do.call(rbind, lapply(scoreList, function(x){data.frame(cbind(info(x,2), readPlate_worms(x, 60, 2000), step="score"))}))

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
summarizedScoreData <- rawScoreData %>% group_by(date, experiment, round, assay, plate, drug) %>% do(summarizePlate_worms(., strains=eval(strainsData[strainsData$assay==.$assay[1],2])[[1]], quantiles=TRUE))

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
completeData <- completeData %>% group_by(assay, plate) %>% do(removeWells(., unlist(contamination[as.character(contamination$assay)==.$assay[1] & as.numeric(as.character(contamination$plate))==as.numeric(.$plate[1]),3])))

#NA out wash wells
completeData[is.na(completeData$strain), which(colnames(completeData)=="n"):ncol(completeData)] <- NA

#NA out any wells where norm.n is infinite (i.e. nothing sorted to well)
completeData[is.infinite(completeData$norm.n) | is.na(completeData$norm.n), which(colnames(completeData)=="n"):ncol(completeData)] <- NA

if(generateReports){
    completeData %>% group_by(assay, plate) %>% do(data.frame(scoreReport(., contamination)))
}

controls <- do.call(rbind, lapply(lapply(directories, function(x){list.files(x, pattern="controls", recursive=TRUE, full.names=TRUE)}),
                                  function(x){
                                      source(x)
                                      assay <- info(x)$assay
                                      data.frame(assay, control=I(controlPlates), plates=I(testPlates))
                                  }))

finalData <- completeData %>% group_by(drug) %>% do(regress(., completeData, controls)) %>% arrange(assay)

experiment <- info(directories[1], levels=0)$experiment
round <- info(directories[1], levels=0)$round

write.csv(finalData, paste0("~/Dropbox/HTA/Results/ProcessedData/", experiment, round, "_complete_simple.csv"), row.names=FALSE)
