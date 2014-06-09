#Load the required packages
require(tools)
require(plyr)
require(kernlab)
require(stringr)
require(markdown)
require(knitr)
require(ggplot2)
require(reshape)
require(dplyr)
require(COPASutils)

options(echo=TRUE)

args <- commandArgs(trailingOnly = TRUE)
makeReports <- as.logical(args[1])
dataDirs <- args[2:length(args)]

#Set directories here:
##########################
#Root directory (usually the home directory, otherwise it should be 2 steps above the data files)
dir.root <- "~"

#Path to sorter processor functions minus root dir
file.fxns <- "SorterDataAssembly/ProcessSorterFxns_NU_TCS.R"

#Path to presentation style data file minus root dir
file.pres <- "SorterDataAssembly/PresentationStyle.Rdata"

#Path to contamination file minus root dir
file.contam <- "contamination.R"

#Path to strains file minus root dir
file.strains <- "strains.R"

#Path to controls file minus root dir
file.controls <- "controls.R"

#Path to setup dir minus root dir
dir.setup <- "setup"

#Path to score dir minus root dir
dir.score <- "score"

#Whole paths to the setup and score report markdown templates
file.report.setup <- "~/SorterDataAssembly/MasterSetupReport.Rmd"
file.report.score <- "~/SorterDataAssembly/MasterScoreReport.Rmd"
##########################

completeDFs <- list()

for(dir in seq(1,length(dataDirs))){
    dir.data <- dataDirs[dir]
    
    
    #Next five lines establish the output directories for the reports, results, and temporary files
    ##If you change these, it will break the code as it is written, you will need to change some of the
    ##directories that are hard coded in later on
    dir.existing <- file.path(dir.root,dir.data)
    dir.new = "reports"
    dir.create(file.path(dir.existing,dir.new))
    dir.report<-file.path(dir.existing,dir.new)
    dir.create(file.path(dir.existing,"results"))
    dir.create(file.path(dir.existing,"temp"))
    

    #Establish the user's working directory
    setwd(dir.root)
    
    
    #Source the data to access the enclosed functions
    source(file.path(dir.root,file.fxns))
    
    #Get the user's plot style guide file 
    load(file.path(dir.root,file.pres))
    
    #Change the user's working directory to that which contains the data from the sorter
    setwd(file.path(dir.root,dir.data))
    
    #Load the contamination file
    source(file.path(dir.root,dir.data,file.contam))
    
    #Load the strain data
    source(file.path(dir.root,dir.data,file.strains))
    
    #Arrange the strain data into plate format
    strains<-matrix(strains,nrow=8,ncol=12,byrow=TRUE)
    
    #Change working directory to the setup data folder
    setwd(file.path(dir.root,dir.data,dir.setup))
    
    #Get the plate numbers from the file names in the setup directory
    setup.filelist<-dir(pattern="*.txt")
    setup.plate<-llply(setup.filelist,function(x){plateno(x)})
    
    
    
    #For every file in the above list of files, process the setup file with procSetup
    setup.df<-llply(setup.filelist,function(x){procSetup(x)})
    
    fileName <- file.path(dir.root,dir.data,"results",paste0(experimentName(dir.data),"_rawSetupData.Rds"))
    saveRDS(setup.df, fileName)
    
    
    #Do the same with the scoring data
    setwd(file.path(dir.root,dir.data,dir.score))
    score.filelist<-dir(pattern="*.txt")
    score.plate<-llply(score.filelist,function(x){plateno(x)})
    
    
    #Create a list of scoring results with bubbles sccounted for with the SVM
    score.modplate<-llply(score.filelist,function(x){readPlate_worms(x)})
    
    fileName <- file.path(dir.root,dir.data,"results",paste0(experimentName(dir.data),"_rawScoringData.Rds"))
    saveRDS(score.modplate, fileName)
    
    setup.filelist<-as.list(setup.filelist)
    names(setup.filelist)<-setup.plate
    names(setup.df)<-setup.plate
    
    
    #CHECKS FOR MISSING SCORE PLATES (removes setup plates which do not correspond)
    if(length(score.plate)<(length(setup.plate)))
    {
        for (plate in score.plate)
        {
            for (i in 1:(length(score.plate)))
            {
                if(plate==score.plate[[i]])
                {
                    setup.plate[[i]]<-plate
                    setup.filelist[[i]]<-setup.filelist[[paste0(plate)]]
                    setup.df[[i]]<-setup.df[[paste0(plate)]]
                }
            }
        }
        length(setup.plate)<-length(score.plate)
        length(setup.filelist)<-length(score.filelist)
        length(setup.df)<-length(score.plate)
    }
    
    
    #Generate all of the reports and data from the setup files
    for(i in 1:(length(setup.plate))) {
        dir.create(file.path(dir.existing,"temp"))
        setwd(file.path(dir.existing,"temp"))
        
        setup.proc<-setup.df[[i]]
        saveRDS(setup.proc,file=file.path(dir.existing,"temp","setup-proc.rds"))
        
        setup.proc[is.na(setup.proc)]<--1
        setup.proc[-1:-4]<-round(setup.proc[,-1:-4],1) 
     
        split <- setup.plate[[i]]
        plate<-paste0("p",split)
        
        if(isTRUE(exists(plate))){
            contam<-get(plate)
        } else {
            contam<-"NA"
        }
        
        if(makeReports){
            file.setup <- setup.filelist[[i]]
            saveRDS(file.setup,file=file.path(dir.existing,"temp","file-setup.rds"))
            
            fileName <- strsplit(file.setup, "\\.")[[1]][1]
            
            date=Sys.Date()
            date=as.character(format(date,format="%Y%m%d"))
            saveRDS(date,file=file.path(dir.existing,"temp","date.rds"))
            
            plot.setup.sorted<- ggplot(setup.proc)+geom_rect(aes(xmin=0,xmax=5,ymin=0,ymax=5,fill=sorted))+
                scale_fill_gradient2(high = "green", low = "red", mid = "yellow", midpoint = ceiling((max(setup.proc$sorted)-min(setup.proc$sorted))/2))+
                facet_grid(row~col)+geom_text(aes(x=2.5,y=2.5,label=sorted))+
                presentation+
                theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+
                xlab("columns")+ylab("rows")+labs(title=paste0("Setup p",setup.plate[[i]]," # Sorted"))
            saveRDS(plot.setup.sorted,file=file.path(dir.existing,"temp","plot-setup-sorted.rds"))
            
            for(row in 1:nrow(setup.proc)){
                if(setup.proc$pop[row] < 15) {
                    setup.proc$bin[row] = "low"
                } else if(setup.proc$pop[row] < 25) {
                    setup.proc$bin[row] = "low-mid"
                } else if(setup.proc$pop[row] < 50) {
                    setup.proc$bin[row] = "high-mid"
                } else {
                    setup.proc$bin[row] = "high"
                }
            }
            
            plot.setup.pop <- ggplot(setup.proc)+geom_rect(aes(xmin=0,xmax=5,ymin=0,ymax=5,fill=bin))+
                scale_fill_manual(values = c("low" = "white", "low-mid" = "green", "high-mid" = "yellow", "high" = "red"))+
                facet_grid(row~col)+geom_text(aes(x=2.5,y=2.5,label=pop))+
                presentation+
                theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+
                xlab("columns")+ylab("rows")+labs(title=paste0("Setup p",setup.plate[[i]]," Population"))
            
            saveRDS(plot.setup.pop,file=file.path(dir.existing,"temp","plot-setup-pop.rds"))
            
            plot.setup.tofext<-ggplot(setup.proc)+geom_rect(fill=NA,aes(xmin=0,xmax=5,ymin=0,ymax=5))+facet_grid(row~col)+geom_text(aes(x=1,y=4,label=TOF))+geom_text(aes(x=4,y=4,label=EXT))+geom_text(aes(x=1,y=1,label=TOFmed))+geom_text(aes(x=4,y=1,label=EXTmed))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title="Setup p",setup.plate[[i]]," TOF/EXT")
            saveRDS(plot.setup.tofext,file=file.path(dir.existing,"temp","plot-setup-tofext.rds"))
            
            saveRDS(strains, file=file.path(dir.existing,"temp","strains.rds"))
            
            saveRDS(split,file=file.path(dir.existing,"temp","split.rds"))
            saveRDS(contam,file=file.path(dir.existing,"temp","contam.rds"))
            
            fileNameString <- paste0(date,'-',split,'.md')
            
            knit(file.report.setup, file.path(dir.existing,"temp",fileNameString)) 
            markdownToHTML(fileNameString, file.path(dir.report,paste0(fileName,'_setup.html')))
        }
        unlink(file.path(dir.existing,"temp"), recursive = TRUE)
        print(paste0("Done with setup #", i))
    }
    
    
    
    
    #score.pheno = list of processed score datasets with phenotype information (TOF/EXT quantiles etc)
    #without t(strains), strain matrix won't match up with plate setup
    score.pheno<-llply(score.modplate,function(y){summarizePlate_worms(y,t(strains), quantiles=TRUE, log=TRUE)})
    
    
    #adds to score.pheno datasets a column of n normalized by number sorted in setup
    for(i in 1:length(score.plate))
    {
        setup<-setup.df[[i]]
        score<-score.pheno[[i]]
        n <- score$n
        sorted <- setup$sorted
        if (length(n) != length(sorted)){
            stop("The lengths of the setup and score data frames do not match.\nRun PlateStitcher.py on both the setup and score data and try again.")
        }
        norm<-score$n/setup$sorted
        norm<-ifelse(is.infinite(norm),NA,norm)
        score$norm.n<-norm
        score.pheno[[i]]<-score
    }
    

    
    source(file.path(dir.root,dir.data,file.contam))
    
    bad<-score.plate
    
    contamWells <- list()
    
    for(i in 1:length(bad))
    {
        bad[[i]]<-paste0("p",bad[[i]])
        plateNumber = bad[[i]]
        if(exists(bad[[i]])){
            bad[[i]] = get(bad[[i]])
        } else {
            bad[[i]] = ""
        }
        contamWells[[i]] = bad[[i]]
        try(rm(list = plateNumber))
        badWellsDF <- setup.df[[i]][setup.df[[i]]$sorted == 0,c(1,2)]
        for(row in seq(1, nrow(badWellsDF))){
            badWell <- paste0(badWellsDF[row,1], badWellsDF[row,2])
            bad[[i]] = append(bad[[i]], badWell)
        }
        i = i+1
    }
    
    #score.proc is list of totally processed score files
    for (i in 1:length(score.pheno)) {
        score.pheno[[i]] <- removeWells(score.pheno[[i]], bad[[i]])
    }
    
    date=Sys.Date()
    date=as.character(format(date,format="%Y%m%d"))
    
    melted.score.pheno<-llply(score.pheno,function(x){meltdf(x)})
    
    #Create complete data frame
    score.info <- llply(file.path(dir.data,score.filelist), function(x){info(x)})
    for(i in 1:length(score.pheno)){
        score.pheno[[i]] = as.data.frame(cbind(score.info[[i]],score.pheno[[i]]))
    }
    
    #Print out final reports and assemble final data frames
    for(i in 1:(length(score.plate)))
    {
        dir.create(file.path(dir.existing,"temp"))
        setwd(file.path(dir.existing,"temp"))
        
        
        file.score<-score.filelist[[i]]
        saveRDS(file.score,file=file.path(dir.existing,"temp","file-score.rds"))
        
        file.name <- strsplit(file.score, "\\.")[[1]][1]
        
        saveRDS(strains,file=file.path(dir.existing,"temp","strains.rds"))
        
        
        date=Sys.Date()
        date=as.character(format(date,format="%Y%m%d"))
        saveRDS(date,file=file.path(dir.existing,"temp","date.rds"))
        
        
        proc<-score.pheno[[i]]
        saveRDS(proc,file=file.path(dir.existing,"temp","proc.rds"))
        
        naCount <- nrow(proc[as.numeric(as.numeric(proc$col)) %% 2 == 1 & is.na(proc$n),])
        naStrains <- as.character(proc[as.numeric(proc$col) %% 2 == 1 & is.na(proc$n),]$strain)
        saveRDS(naCount, file=file.path(dir.existing,"temp","naCount.rds"))
        saveRDS(naStrains, file=file.path(dir.existing,"temp","naStrains.rds"))
        
        possibleContam <- possContam(proc)
        
        saveRDS(possibleContam,file=file.path(dir.existing,"temp","possibleContam.rds"))
        
        
        fileName <- paste0(i,".csv")
        write.csv(proc,file <- file.path(dir.existing,"results",fileName))
        
        
        split<-score.plate[[i]]
        saveRDS(split,file=file.path(dir.existing,"temp","split.rds"))
        
        plate<-paste0("p",split)
        contam<-contamWells[[i]]
        saveRDS(contam,file=file.path(dir.existing,"temp","contam.rds"))
        
        proc[is.na(proc)]<--1
        proc$norm.n<-round(proc$norm.n,0)
        proc$mean.norm.red<-round(proc$mean.norm.red,2)
        proc$median.norm.red<-round(proc$median.norm.red,2)
        proc$mean.TOF<-round(proc$mean.TOF,1)
        proc$mean.EXT<-round(proc$mean.EXT,1)
        
        
        
        if(makeReports){
            
            plot.score.pop<-ggplot(proc)+geom_rect(aes(xmin=0,xmax=5,ymin=0,ymax=5,fill=norm.n))+facet_grid(row~col)+geom_text(aes(x=2.5,y=2.5,label=norm.n,colour="white"))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Score p",score.plate[[i]]," Population"))
            
            saveRDS(plot.score.pop,file=file.path(dir.existing,"temp","plot-score-pop.rds"))
            
            plot.score.red<-ggplot(proc)+geom_rect(aes(xmin=0,xmax=5,ymin=0,ymax=5),fill=NA)+facet_grid(row~col)+geom_text(aes(x=1,y=4,label=mean.norm.red))+geom_text(aes(x=4,y=1,label=median.norm.red))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Score p",score.plate[[i]]," Red Fluorescence"))
            saveRDS(plot.score.red,file=file.path(dir.existing,"temp","plot-score-red.rds"))
            
            plot.score.tofext<-ggplot(proc)+geom_rect(fill=NA,aes(xmin=0,xmax=5,ymin=0,ymax=5))+facet_grid(row~col)+geom_text(aes(x=1,y=4,label=mean.TOF))+geom_text(aes(x=4,y=4,label=mean.EXT))+geom_text(aes(x=1,y=1,label=median.TOF))+geom_text(aes(x=4,y=1,label=median.EXT))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Score p",score.plate[[i]]," TOF/EXT"))
            saveRDS(plot.score.tofext,file=file.path(dir.existing,"temp","plot-score-tofext.rds"))
            
            melted.proc<-melted.score.pheno[[i]]
            plot.score.dist<-ggplot(melted.proc,aes(as.factor(col),value,fill=variable))+geom_bar(aes(x=3),stat="identity",position="stack")+facet_grid(row~col)+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Score p",score.plate[[i]]," Life-stage Distribution"))
            saveRDS(plot.score.dist,file=file.path(dir.existing,"temp","plot-score-dist.rds"))
            
            fileNameString <- paste0(date,'-',split,'-score.md')
            
            knit(file.report.score, file.path(dir.existing,"temp",fileNameString))
            markdownToHTML(fileNameString, file.path(dir.report,paste0(file.name,"_score.html")))
        }
        unlink(file.path(dir.existing,"temp"), recursive = TRUE)
        print(paste0("Done with score #", i))
    } 
    
    
    
    filelist <- dir(file.path(dir.existing,"results"), pattern = "[0-9]+.csv")
    
    df <- read.csv(file.path(dir.existing,"results",filelist[1]))
    
    for(i in 2:length(filelist)){
        df <- rbind.fill(df, read.csv(file.path(dir.existing,"results",filelist[i])))
    }
    
    splitfp <- strsplit(dir.data,"/")
    dirName <- splitfp[[1]][(length(splitfp[[1]]))]
    
    details <- strsplit(dirName,"_")[[1]][2]
    
    name <- paste0(details, "_complete.csv")
    
    write.csv(df, file.path(dir.existing,"results",name), row.names = FALSE)

    completeDFs[[dir]] = df
    
    file.remove(file.path(dir.existing,"results",filelist))
}

completeDF <- ldply(completeDFs)
completeDF <- completeDF[,2:ncol(completeDF)]

# Remove infinite values from taking log of columns with 0 as value
completeDF <- data.frame(lapply(completeDF, function(x) replace(x, is.infinite(x),NA)))

startCol <- 10 

assays <- dlply(completeDF, .variables = "assay")

plateData = list()

controlsData = list()

for(dir in seq(1,length(dataDirs))){
    dir.data <- dataDirs[dir]
    data <- assays[[dir]]
    
    columnNames = colnames(data)
    
    plateData = append(plateData, list(data))
    
    #Creating dataframes of control population and growth (n, TOF/EXT 25/50/75/mean)
    source(file.path(dir.root, dir.data, file.controls))
    
    controlDFs <- list()
    controlData = list()
    length(controlData) = length(unique(data$plate))
    
    assay = unique(data$assay)
    
    if (length(controlPlates) != 0){
        for(i in 1:length(controlPlates)){
            controls <- controlPlates[[i]]
            means <- list()
            length(means) = ncol(data)
            plates <- data[data$plate %in% controls,]
            meanPlate = plates %>% group_by(row, col) %>% do(data.frame(means = colMeans(.[10:ncol(.)], na.rm = TRUE)))
            meanPlate = as.data.frame(matrix(meanPlate$means, ncol=ncol(plates)-9, byrow=TRUE))
            controlPlate = as.data.frame(cbind(data.frame(matrix(nrow=96, ncol=9)), meanPlate))
            
            for(j in controlPlates[[i]]){
                controlData[[j]] = data.frame(matrix(nrow=96, ncol=ncol(data)))
                colnames(controlData[[j]]) = columnNames
            }
            for(j in testPlates[[i]]){
                controlData[[j]] = controlPlate
                colnames(controlData[[j]]) = columnNames
            }
        }
        for(i in 1:length(unique(data$plate))){
            if(is.null(controlData[[i]])){
                controlData[[i]] = data.frame(matrix(nrow=96, ncol=ncol(data)))
                colnames(controlData[[i]]) = columnNames
            }
        }
        controlsData = append(controlsData, list(do.call(rbind,controlData)))
    }
}



plateData = data.frame(do.call(rbind, plateData))
controlsData = data.frame(do.call(rbind, controlsData))

plateData[is.na(plateData$strain), 10:ncol(plateData)] = NA
controlsData[controlsData$col %% 2 == 0 | is.na(controlsData$col), 10:ncol(controlsData)] = NA

#important to make sure that the plate data matches the order of the control data
plateData = plateData[order(plateData$assay, plateData$plate),]

completePlates = list()

for(i in unique(plateData$drug)){
    plate = plateData[plateData$drug == i,]
    control = controlsData[plateData$drug == i,]
    for(j in startCol:ncol(plate)){
        residuals = tryCatch({residuals(lm(plate[,j]~plate$assay+control[,j], na.action = na.exclude))},
                             error = function(err){tryCatch({residuals(lm(plate[,j]~control[,j], na.action = na.exclude))},
                                                            error = function(err){return(rep(NA, nrow(plate)))})})
        plate = as.data.frame(cbind(plate, residuals))
        colnames(plate)[ncol(plate)] = paste0("resid.", colnames(plate)[j])
    }
    completePlates = append(completePlates, list(plate))
}

finalDF = ldply(completePlates)
finalDF = finalDF[order(finalDF$assay, finalDF$plate),]

nameFrame = info(dir.data, 0)
fileName = paste0(nameFrame$experiment[1], nameFrame$round[1], "_complete.csv")

write.csv(finalDF, file.path("~/Dropbox/HTA/Results/ProcessedData", fileName), row.names=FALSE)
