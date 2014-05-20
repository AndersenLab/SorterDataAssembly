#Load the required packages
require(tools)
require(plyr)
require(kernlab)
require(stringr)
require(markdown)
require(knitr)
require(ggplot2)
require(reshape)


###############Change for each experiment###################

#Path to experiment folder minus root dir
dataDirs = c("Dropbox/HTA/Results/20140422_N2timecourse")

#Set to false if you want to skip making the reports (saves a lot of time)
makeReports = TRUE

###############Change for each experiment###################






#Set directories here:
##########################
#Root directory (usually the home directory, otherwise it should be 2 steps above the data files)
dir.root = "~"

#Path to sorter processor functions minus root dir
file.fxns = "SorterDataAssembly/ProcessSorterFxns_NU_TCS.R"

#Path to presentation style data file minus root dir
file.pres = "SorterDataAssembly/PresentationStyle.Rdata"

#Path to contamination file minus root dir
file.contam = "contamination.R"

#Path to strains file minus root dir
file.strains = "strains.R"

#Path to controls file minus root dir
file.controls = "controls.R"

#Path to setup dir minus root dir
dir.setup = "setup"

#Path to score dir minus root dir
dir.score = "score"

#Whole paths to the setup and score report markdown templates
file.report.setup = "~/SorterDataAssembly/MasterSetupReport.Rmd"
file.report.score = "~/SorterDataAssembly/MasterScoreReport.Rmd"
##########################

completeDFs = list()

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

experimentName = function(filePath){
    splitfp = strsplit(filePath,"/")
    dirName = splitfp[[1]][(length(splitfp[[1]]))]
    details = strsplit(dirName,"_")[[1]][2]
    return(details)
}

#Extract the metadata info
info = function(filePath){
    splitfp = strsplit(filePath,"/")
    dirName = splitfp[[1]][(length(splitfp[[1]])-1)]
    
    date = strsplit(dirName,"_")[[1]][1]
    
    details = strsplit(dirName,"_")[[1]][2]
    
    experiment = strsplit(details,"[0-9]+")[[1]][1]
    round = strsplit(details,"(?i)[a-z]+")[[1]][2]
    assay = strsplit(details,"[0-9]+")[[1]][2]
    
    split = strsplit(splitfp[[1]][(length(splitfp[[1]]))],"_")[[1]]
    drug = strsplit(split[2],"\\.")[[1]][1]
    plate = strsplit(split[1],"p")[[1]][2]
    
    frame = data.frame(date,experiment,round,assay,plate,drug)
    
    return(frame)
}

#Function to process the setup data
procSetup <- function(file, tofmin=20, tofmax=2000, extmin=20, extmax=5000) {
    
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

possContam = function(procDataFrame){
    strainMean = mean(procDataFrame$n[!is.na(procDataFrame$strain)], na.rm = TRUE)
    strainSD = sd(procDataFrame$n[!is.na(procDataFrame$strain)], na.rm = TRUE)
    washMean = mean(procDataFrame$n[is.na(procDataFrame$strain)], na.rm = TRUE)
    washSD = sd(procDataFrame$n[is.na(procDataFrame$strain)], na.rm = TRUE)
    possibleContam = c()
    for(j in seq(1,nrow(procDataFrame),)){
        if(!is.na(procDataFrame[j,"strain"] & !is.na(procDataFrame[j,"n"]))){
            if(procDataFrame[j,"n"] > strainMean + (2*strainSD)){
                row = as.character(procDataFrame[j, "row"])
                col = as.numeric(procDataFrame[j, "col"])
                adjacentWash = procDataFrame[procDataFrame$row==row & procDataFrame$col==(col+1),"n"]
                if(adjacentWash > washMean + (2*washSD)){
                    possibleContam = append(possibleContam, paste0(row, col))
                }
            }
        }
    } 
}









for(dir in seq(1,length(dataDirs))){
    dir.data = dataDirs[dir]
    
    
    #Next five lines establish the output directories for the reports, results, and temporary files
    ##If you change these, it will break the code as it is written, you will need to change some of the
    ##directories that are hard coded in later on
    dir.existing = file.path(dir.root,dir.data)
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
    
    fileName = file.path(dir.root,dir.data,"results",paste0(experimentName(dir.data),"_rawSetupData.Rds"))
    saveRDS(setup.df, fileName)
    
    
    #Do the same with the scoring data
    setwd(file.path(dir.root,dir.data,dir.score))
    score.filelist<-dir(pattern="*.txt")
    score.plate<-llply(score.filelist,function(x){plateno(x)})
    score.df<-llply(score.filelist,function(x){readSorter(x)})
    
    
    
    #Create a list of scoring results with bubbles sccounted for with the SVM
    score.modplate<-llply(score.filelist,function(x){sortertoDF(x)})
    
    fileName = file.path(dir.root,dir.data,"results",paste0(experimentName(dir.data),"_rawScoringData.Rds"))
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
        
        
        
        pop1 <- subset(setup.proc, pop>=50)
        pop2 <- subset(setup.proc, pop>=25&pop<50)
        pop3 <- subset(setup.proc,pop>=15&pop<25)
        pop4 <- subset(setup.proc,pop<15)
     
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
            
            plot.setup.pop <- ggplot(pop1)+geom_rect(aes(xmin=0,xmax=5,ymin=0,ymax=5), fill="red")+geom_rect(data=pop2,aes(xmin=0,xmax=5,ymin=0,ymax=5),fill="yellow")+geom_rect(data=pop3,aes(xmin=0,xmax=5,ymin=0,ymax=5),fill="green")+geom_rect(data=pop4,aes(xmin=0,xmax=5,ymin=0,ymax=5),fill="white")+facet_grid(row~col) +geom_text(aes(x=2.5,y=2.5,label=pop))+geom_text(data=pop2,aes(x=2.5,y=2.5,label=pop))+geom_text(data=pop3,aes(x=2.5,y=2.5,label=pop))+geom_text(data=pop4,aes(x=2.5,y=2.5,label=pop))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Setup p",setup.plate[[i]]," Population"))
            saveRDS(plot.setup.pop,file=file.path(dir.existing,"temp","plot-setup-pop.rds"))
            
            plot.setup.tofext<-ggplot(setup.proc)+geom_rect(fill=NA,aes(xmin=0,xmax=5,ymin=0,ymax=5))+facet_grid(row~col)+geom_text(aes(x=1,y=4,label=TOF))+geom_text(aes(x=4,y=4,label=EXT))+geom_text(aes(x=1,y=1,label=TOFmed))+geom_text(aes(x=4,y=1,label=EXTmed))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title="Setup p",setup.plate[[i]]," TOF/EXT")
            saveRDS(plot.setup.tofext,file=file.path(dir.existing,"temp","plot-setup-tofext.rds"))
            
            saveRDS(strains, file=file.path(dir.existing,"temp","strains.rds"))
            
            saveRDS(split,file=file.path(dir.existing,"temp","split.rds"))
            saveRDS(contam,file=file.path(dir.existing,"temp","contam.rds"))
            
            fileNameString = paste0(date,'-',split,'.md')
            
            knit(file.report.setup, file.path(dir.existing,"temp",fileNameString)) 
            markdownToHTML(fileNameString, file.path(dir.report,paste0(split,'_setup.html')))
        }
        unlink(file.path(dir.existing,"temp"), recursive = TRUE)
        print(paste0("Done with setup #", i))
    }
    
    
    
    
    
    
    
    
    
    #score.pheno = list of processed score datasets with phenotype information (TOF/EXT quantiles etc)
    #without t(strains), strain matrix won't match up with plate setup
    score.pheno<-llply(score.modplate,function(x){processPheno(x,t(strains))})
    
    
    
    
    
    #adds to score.pheno datasets a column of n normalized by number sorted in setup
    for(i in 1:length(score.plate))
    {
        setup<-setup.df[[i]]
        score<-score.pheno[[i]]
        n = score$n
        sorted = setup$sorted
        if (length(n) != length(sorted)){
            stop("The lengths of the setup and score data frames do not match.\nRun PlateStitcher.py on both the setup and score data and try again.")
        }
        norm<-score$n/setup$sorted
        norm<-ifelse(is.infinite(norm),NA,norm)
        score$norm.n<-norm
        score.pheno[[i]]<-score
    }
    
    bad<-score.plate
    for(i in 1:length(bad))
    {
        bad[[i]]<-paste0("p",bad[[i]])
        bad[[i]]<-get(bad[[i]])
        badWellsDF = setup.df[[i]][setup.df[[i]]$sorted == 0,c(1,2)]
        for(row in seq(1, nrow(badWellsDF))){
            badWell = paste0(badWellsDF[row,1], badWellsDF[row,2])
            bad[[i]] = append(bad[[i]], badWell)
        }
    }
    
    #score.proc is list of totally processed score files
    for (i in 1:length(score.pheno)) {
        score.pheno[[i]] <- removeWells(score.pheno[[i]], bad[[i]])
    }
    
    
    date=Sys.Date()
    date=as.character(format(date,format="%Y%m%d"))
    
    
    score.pheno[is.na(score.pheno)]<-0
    
    melted.score.pheno<-llply(score.pheno,function(x){meltdf(x)})
    
    
    
    #Create complete data frame
    score.info = llply(file.path(dir.data,score.filelist), function(x){info(x)})
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
        
        
        saveRDS(strains,file=file.path(dir.existing,"temp","strains.rds"))
        
        
        date=Sys.Date()
        date=as.character(format(date,format="%Y%m%d"))
        saveRDS(date,file=file.path(dir.existing,"temp","date.rds"))
        
        
        proc<-score.pheno[[i]]
        saveRDS(proc,file=file.path(dir.existing,"temp","proc.rds"))
        
        possibleContam = possContam(proc)
        
        saveRDS(possibleContam,file=file.path(dir.existing,"temp","possibleContam.rds"))
        
        
        fileName = paste0(i,".csv")
        write.csv(proc,file = file.path(dir.existing,"results",fileName))
        
        
        split<-score.plate[[i]]
        saveRDS(split,file=file.path(dir.existing,"temp","split.rds"))
        
        plate<-paste0("p",split)
        contam<-get(plate)
        saveRDS(contam,file=file.path(dir.existing,"temp","contam.rds"))
        
        proc[is.na(proc)]<--1
        proc$norm.n<-round(proc$norm.n,0)
        proc$mean.normred<-round(proc$mean.normred,2)
        proc$median.normred<-round(proc$median.normred,2)
        proc$meanTOF<-round(proc$meanTOF,1)
        proc$meanEXT<-round(proc$meanEXT,1)
        
        
        
        if(makeReports){
            
            plot.score.pop<-ggplot(proc)+geom_rect(aes(xmin=0,xmax=5,ymin=0,ymax=5,fill=norm.n))+facet_grid(row~col)+geom_text(aes(x=2.5,y=2.5,label=norm.n,colour="white"))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Score p",score.plate[[i]]," Population"))
            
            saveRDS(plot.score.pop,file=file.path(dir.existing,"temp","plot-score-pop.rds"))
            
            plot.score.red<-ggplot(proc)+geom_rect(aes(xmin=0,xmax=5,ymin=0,ymax=5),fill=NA)+facet_grid(row~col)+geom_text(aes(x=1,y=4,label=mean.normred))+geom_text(aes(x=4,y=1,label=median.normred))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Score p",score.plate[[i]]," Red Fluorescence"))
            saveRDS(plot.score.red,file=file.path(dir.existing,"temp","plot-score-red.rds"))
            
            plot.score.tofext<-ggplot(proc)+geom_rect(fill=NA,aes(xmin=0,xmax=5,ymin=0,ymax=5))+facet_grid(row~col)+geom_text(aes(x=1,y=4,label=meanTOF))+geom_text(aes(x=4,y=4,label=meanEXT))+geom_text(aes(x=1,y=1,label=medianTOF))+geom_text(aes(x=4,y=1,label=medianEXT))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Score p",score.plate[[i]]," TOF/EXT"))
            saveRDS(plot.score.tofext,file=file.path(dir.existing,"temp","plot-score-tofext.rds"))
            
            melted.proc<-melted.score.pheno[[i]]
            plot.score.dist<-ggplot(melted.proc,aes(as.factor(col),value,fill=variable))+geom_bar(aes(x=3),stat="identity",position="stack")+facet_grid(row~col)+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Score p",score.plate[[i]]," Life-stage Distribution"))
            saveRDS(plot.score.dist,file=file.path(dir.existing,"temp","plot-score-dist.rds"))
            
            fileNameString = paste0(date,'-',split,'-score.md')
            
            knit(file.report.score, file.path(dir.existing,"temp",fileNameString))
            markdownToHTML(fileNameString, file.path(dir.report,paste0(split,"_score.html")))
        }
        unlink(file.path(dir.existing,"temp"), recursive = TRUE)
        print(paste0("Done with score #", i))
    } 
    
    
    
    filelist = dir(file.path(dir.existing,"results"), pattern = "[0-9]+.csv")
    
    df = read.csv(file.path(dir.existing,"results",filelist[1]))
    
    for(i in 2:length(filelist)){
        df = rbind.fill(df, read.csv(file.path(dir.existing,"results",filelist[i])))
    }
    
    splitfp = strsplit(dir.data,"/")
    dirName = splitfp[[1]][(length(splitfp[[1]]))]
    
    details = strsplit(dirName,"_")[[1]][2]
    
    name = paste0(details, "_complete.csv")
    
    write.csv(df, file.path(dir.existing,"results",name), row.names = FALSE)

    completeDFs[[dir]] = df
    
    file.remove(file.path(dir.existing,"results",filelist))
}

completeDF = ldply(completeDFs)
completeDF = completeDF[,2:ncol(completeDF)]
if(length(levels(completeDF$assay)) > 1){
    startCol = ncol(completeDF) + 1
    for (i in seq(10,ncol(completeDF))){
        colName = colnames(completeDF)[i]
        vectorName = paste0("resid.assay.",colName)
        assign("vector", residuals(lm(completeDF[,i]~completeDF$assay, na.action = na.exclude)))
        completeDF = cbind(completeDF, vector)
        colnames(completeDF)[ncol(completeDF)] = vectorName
    }
} else {
    startCol = 10 
}



assays = dlply(completeDF, .variables = "assay")
for(i in seq(1, length(assays))){
    assays[[i]] = dlply(assays[[i]], .variables = "plate")
}

output = list()

for(dir in seq(1,length(dataDirs))){
    dir.data = dataDirs[dir]
    data = assays[[dir]]
    
    #Creating dataframes of control population and growth (n, TOF/EXT 25/50/75/mean)
    source(file.path(dir.root, dir.data, file.controls))
    
    controlDFs = list()
    
    finalPlates = list()
    length(finalPlates) = length(data)
    if (length(controlPlates != 0)){
        for(i in 1:length(controlPlates)){
            controls = controlPlates[[i]]
            means = list()
            length(means) = ncol(data[[controls[1]]])
            for(j in 1:length(controls)){
                plate = data[[controls[j]]]
                finalPlates[[controls[j]]] = plate
                for(k in 10:ncol(plate)){
                    means[[k]] = cbind(means[[k]], plate[,k])
                }
            }
            controlDF = plate[,1:9]
            for(j in 10:ncol(data[[controls[1]]])){
                controlDF = as.data.frame(cbind(controlDF, rowMeans(means[[j]], na.rm = TRUE)))
            }
            colnames(controlDF) = colnames(data[[controls[1]]])
            controlDF = as.data.frame(controlDF)
            controlDFs = append(controlDFs, list(controlDF))
        }
        
        rsq = data.frame()
        
        for(i in 1:length(testPlates)){
            plates = testPlates[[i]]
            controlPlate = controlDFs[[i]]
            for(j in 1:length(plates)){
                plate = data[[plates[j]]]
                for(k in startCol:ncol(plate)){
                    plate = cbind(plate, residuals(lm(plate[,k]~controlPlate[,k], na.action = na.exclude)),
                                  residuals(lm(plate[,k]~plate$n, na.action = na.exclude)),
                                  residuals(lm(plate[,k]~plate$n+controlPlate[,k], na.action = na.exclude)) #Do you want additive or interaction effects??
                    )
                    colnames(plate)[ncol(plate)-2] = paste0("resid.control.",colnames(plate)[k])
                    colnames(plate)[ncol(plate)-1] = paste0("resid.n.",colnames(plate)[k])
                    colnames(plate)[ncol(plate)] = paste0("resid.control_n.",colnames(plate)[k])
                    fit1 = summary(lm(plate[,k]~controlPlate[,k], na.action = na.exclude))
                    rsq = rbind(rsq, data.frame(Variable = paste0("resid.control.",colnames(plate)[k]), RSquared = summary(lm(plate[,k]~controlPlate[,k], na.action = na.exclude))$r.squared,
                                                pval = fit1$coefficients[2,4]))
                    fit2 = summary(lm(plate[,k]~plate$n, na.action = na.exclude))
                    rsq = rbind(rsq, data.frame(Variable = paste0("resid.n.",colnames(plate)[k]), RSquared = summary(lm(plate[,k]~plate$n, na.action = na.exclude))$r.squared,
                                                pval = fit2$coefficients[2,4]))
                    fit3 = summary(lm(plate[,k]~plate$n+controlPlate[,k], na.action = na.exclude))
                    rsq = rbind(rsq, data.frame(Variable = paste0("resid.control_n.",colnames(plate)[k]), RSquared = summary(lm(plate[,k]~plate$n+controlPlate[,k], na.action = na.exclude))$r.squared,
                                                pval = fit3$coefficients[2,4]))
                    
                }
                plate = cbind(plate, residuals(lm(plate$n~controlPlate$n, na.action = na.exclude)),
                              residuals(lm(plate$n~controlPlate$norm.n, na.action = na.exclude)),
                              residuals(lm(plate$norm.n~controlPlate$n, na.action = na.exclude)),
                              residuals(lm(plate$norm.n~controlPlate$norm.n, na.action = na.exclude)))
                colnames(plate)[ncol(plate)-3] = "resid.control.n.resid.assay.n"
                colnames(plate)[ncol(plate)-2] = "resid.control.norm.n.resid.assay.n"
                colnames(plate)[ncol(plate)-1] = "resid.control.n.resid.assay.norm.n"
                colnames(plate)[ncol(plate)] = "resid.control.norm.n.resid.assay.norm.n"
                finalPlates[[plates[j]]] = plate
            }
        }
        output = append(output, finalPlates)
    }
}

if(exists(output)){
    finalDF = ldply(output)
} else {
    finalDF = completeDF
}


#sumRSQ = ddply(rsq, .variables = "Variable", summarize, meanRSQ = mean(RSquared), sdRSQ = sd(RSquared), medianRSQ = median(RSquared), meanP = mean(pval), sdP = sd(pval), medianP = median(pval))
