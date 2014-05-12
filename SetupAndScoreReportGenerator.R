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
dir.data = "SorterDataAssembly"

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

#Path to setup dir minus root dir
dir.setup = "setup"

#Path to score dir minus root dir
dir.score = "score"

#Whole paths to the setup and score report markdown templates
file.report.setup = "~/SorterDataAssembly/MasterSetupReport.Rmd"
file.report.score = "~/SorterDataAssembly/MasterScoreReport.Rmd"


dir.existing = file.path(dir.root,dir.data)

#Next four lines establish the output directories for the reports, results, and temporary files
##If you change these, it will break the code as it is written, you will need to change some of the
##directories that are hard coded in later on
dir.new = "reports"
dir.create(file.path(dir.existing,dir.new))
dir.report<-file.path(dir.existing,dir.new)
dir.create(file.path(dir.existing,"results"))
dir.create(file.path(dir.existing,"temp"))
##########################


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

#For every file in the above list of files, process the setup file with procSetup
setup.df<-llply(setup.filelist,function(x){procSetup(x)})


#Do the same with the scoring data
setwd(file.path(dir.root,dir.data,dir.score))
score.filelist<-dir(pattern="*.txt")
score.plate<-llply(score.filelist,function(x){plateno(x)})
score.df<-llply(score.filelist,function(x){readSorter(x)})

#Convert the sorter data to a data frame
sortertoDF <- function(file, tofmin=20, tofmax=2000, extmin=20, extmax=5000) {
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

#Create a list of scoring results with bubbles sccounted for with the SVM
score.modplate<-llply(score.filelist,function(x){sortertoDF(x)})


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
for(i in 1:(length(setup.plate)))
{
    dir.create(file.path(dir.existing,"temp"))
    setwd(file.path(dir.existing,"temp"))
    
    file.setup <- setup.filelist[[i]]
    saveRDS(file.setup,file=file.path(dir.existing,"temp","file-setup.rds"))
    
    date=Sys.Date()
    date=as.character(format(date,format="%Y%m%d"))
    saveRDS(date,file=file.path(dir.existing,"temp","date.rds"))
    
    setup.proc<-setup.df[[i]]
    saveRDS(setup.proc,file=file.path(dir.existing,"temp","setup-proc.rds"))
    
    setup.proc[is.na(setup.proc)]<-0
    setup.proc[-1:-4]<-round(setup.proc[,-1:-4],1) 
    
    sorted1 <- subset(setup.proc, sorted==1|sorted==0)
    sorted2 <- subset(setup.proc, sorted==2)
    sorted3 <- subset(setup.proc,sorted==3)
    
    plot.setup.sorted<- ggplot(sorted1)+geom_rect(aes(xmin=0,xmax=5,ymin=0,ymax=5),fill="red")+geom_rect(data=sorted2,aes(xmin=0,xmax=5,ymin=0,ymax=5),fill="yellow")+geom_rect(data=sorted3,aes(xmin=0,xmax=5,ymin=0,ymax=5),fill="green")+facet_grid(row~col) +geom_text(aes(x=2.5,y=2.5,label=sorted))+geom_text(data=sorted2,aes(x=2.5,y=2.5,label=sorted))+geom_text(data=sorted3,aes(x=2.5,y=2.5,label=sorted))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Setup p",setup.plate[[i]]," # Sorted"))
    saveRDS(plot.setup.sorted,file=file.path(dir.existing,"temp","plot-setup-sorted.rds"))
    
    
    pop1 <- subset(setup.proc, pop>=50)
    pop2 <- subset(setup.proc, pop>=25&pop<50)
    pop3 <- subset(setup.proc,pop>=15&pop<25)
    pop4 <- subset(setup.proc,pop<15)
    plot.setup.pop <- ggplot(pop1)+geom_rect(aes(xmin=0,xmax=5,ymin=0,ymax=5), fill="red")+geom_rect(data=pop2,aes(xmin=0,xmax=5,ymin=0,ymax=5),fill="yellow")+geom_rect(data=pop3,aes(xmin=0,xmax=5,ymin=0,ymax=5),fill="green")+geom_rect(data=pop4,aes(xmin=0,xmax=5,ymin=0,ymax=5),fill="white")+facet_grid(row~col) +geom_text(aes(x=2.5,y=2.5,label=pop))+geom_text(data=pop2,aes(x=2.5,y=2.5,label=pop))+geom_text(data=pop3,aes(x=2.5,y=2.5,label=pop))+geom_text(data=pop4,aes(x=2.5,y=2.5,label=pop))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Setup p",setup.plate[[i]]," Population"))
    saveRDS(plot.setup.pop,file=file.path(dir.existing,"temp","plot-setup-pop.rds"))
    
    plot.setup.tofext<-ggplot(setup.proc)+geom_rect(fill=NA,aes(xmin=0,xmax=5,ymin=0,ymax=5))+facet_grid(row~col)+geom_text(aes(x=1,y=4,label=TOF))+geom_text(aes(x=4,y=4,label=EXT))+geom_text(aes(x=1,y=1,label=TOFmed))+geom_text(aes(x=4,y=1,label=EXTmed))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title="Setup p",setup.plate[[i]]," TOF/EXT")
    saveRDS(plot.setup.tofext,file=file.path(dir.existing,"temp","plot-setup-tofext.rds"))
    
    saveRDS(strains, file=file.path(dir.existing,"temp","strains.rds"))
    
    split <- setup.plate[[i]]
    plate<-paste0("p",split)
    
    if(isTRUE(exists(plate))){
        contam<-get(plate)
    } else {
        contam<-"NA"
    }
    
    saveRDS(split,file=file.path(dir.existing,"temp","split.rds"))
    saveRDS(contam,file=file.path(dir.existing,"temp","contam.rds"))
    
    
    knit(file.report.setup, file.path(dir.existing,"temp",paste0(date,'-',split,'.md'))) 
    markdownToHTML(paste0(date,'-',split,'.md'), file.path(dir.report,paste0(split,'_setup.html')))
    
    unlink(file.path(dir.existing,"temp"), recursive = TRUE)
}






#added quantiles of minEXT, q05_EXT through q95_EXT, maxEXT
processPheno <- function(modplate, strains) {
    processed <- ddply(.data=modplate[modplate$call50=="worm",], .variables=c("row", "col"),
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
                                          f.ad = length(which(x$stage == "adult"))/length(x$stage)
                       )}, .drop=F)
    analysis <- data.frame(strain = as.character(strains), processed)
    analysis <- analysis[order(analysis$strain),]
    analysis <- analysis[order(analysis$row, analysis$col),]
    analysis = analysis[analysis$meanTOF!=-1,]
    #analysis <- droplevels(na.omit(analysis))
    return(analysis)
}


#score.pheno = list of processed score datasets with phenotype information (TOF/EXT quantiles etc)
#without t(strains), strain matrix won't match up with plate setup
score.pheno<-llply(score.modplate,function(x){processPheno(x,t(strains))})





#adds to score.pheno datasets a column of n normalized by number sorted in setup
for(i in 1:length(score.plate))
{
    setup<-setup.df[[i]]
    score<-score.pheno[[i]]
    norm<-score$n/setup$sorted
    norm<-ifelse(is.infinite(norm),NA,norm)
    score$norm.n<-norm
    score.pheno[[i]]<-score
}



#creating dataframes of control population and growth (n, TOF/EXT 25/50/75/mean)
ctrl.n <- data.frame(n1 = as.numeric(score.pheno[[1]]$norm.n),
                     n2 = as.numeric(score.pheno[[2]]$norm.n),
                     n3 = as.numeric(score.pheno[[3]]$norm.n))
ctrl.n$m<-rowMeans(ctrl.n,na.rm=T)

ctrl.q25_TOF <- data.frame(n1 = as.numeric(score.pheno[[1]]$q25_TOF),
                           n2 = as.numeric(score.pheno[[2]]$q25_TOF),
                           n3 = as.numeric(score.pheno[[3]]$q25_TOF))
ctrl.q25_TOF$m<-rowMeans(ctrl.q25_TOF,na.rm=T)

ctrl.q50_TOF <- data.frame(n1 = as.numeric(score.pheno[[1]]$medianTOF),
                           n2 = as.numeric(score.pheno[[2]]$medianTOF),
                           n3 = as.numeric(score.pheno[[3]]$medianTOF))
ctrl.q50_TOF$m<-rowMeans(ctrl.q50_TOF,na.rm=T)

ctrl.q75_TOF <- data.frame(n1 = as.numeric(score.pheno[[1]]$q75_TOF),
                           n2 = as.numeric(score.pheno[[2]]$q75_TOF),
                           n3 = as.numeric(score.pheno[[3]]$q75_TOF))
ctrl.q75_TOF$m<-rowMeans(ctrl.q75_TOF,na.rm=T)

ctrl.mean_TOF <- data.frame(n1 = as.numeric(score.pheno[[1]]$meanTOF),
                            n2 = as.numeric(score.pheno[[2]]$meanTOF),
                            n3 = as.numeric(score.pheno[[3]]$meanTOF))
ctrl.mean_TOF$m<-rowMeans(ctrl.mean_TOF,na.rm=T)

ctrl.q25_EXT <- data.frame(n1 = as.numeric(score.pheno[[1]]$q25_EXT),
                           n2 = as.numeric(score.pheno[[2]]$q25_EXT),
                           n3 = as.numeric(score.pheno[[3]]$q25_EXT))
ctrl.q25_EXT$m<-rowMeans(ctrl.q25_EXT,na.rm=T)

ctrl.q50_EXT <- data.frame(n1 = as.numeric(score.pheno[[1]]$medianEXT),
                           n2 = as.numeric(score.pheno[[2]]$medianEXT),
                           n3 = as.numeric(score.pheno[[3]]$medianTOF))
ctrl.q50_EXT$m<-rowMeans(ctrl.q50_EXT,na.rm=T)

ctrl.q75_EXT <- data.frame(n1 = as.numeric(score.pheno[[1]]$q75_EXT),
                           n2 = as.numeric(score.pheno[[2]]$q75_EXT),
                           n3 = as.numeric(score.pheno[[3]]$q75_EXT))
ctrl.q75_EXT$m<-rowMeans(ctrl.q75_EXT,na.rm=T)

ctrl.mean_EXT <- data.frame(n1 = as.numeric(score.pheno[[1]]$meanEXT),
                            n2 = as.numeric(score.pheno[[2]]$meanEXT),
                            n3 = as.numeric(score.pheno[[3]]$meanEXT))
ctrl.mean_EXT$m<-rowMeans(ctrl.mean_EXT,na.rm=T)


#adding control condition columns and regressing out control conditions
addCtrl <- function(score)
{
    score$ctrl.n <- ctrl.n$m
    score[is.na(score)]<-0
    resid.n <- residuals(lm(score$norm.n~score$ctrl.n))
    score$resid.n <- resid.n
    
    score$ctrl.q25_TOF <- ctrl.q25_TOF$m
    resid.q25_TOF <- residuals(lm(score$q25_TOF~score$ctrl.q25_TOF))
    score$resid.q25_TOF <- resid.q25_TOF
    
    score$ctrl.q50_TOF <- ctrl.q50_TOF$m
    resid.q50_TOF <- residuals(lm(score$medianTOF~score$ctrl.q50_TOF))
    score$resid.q50_TOF <- resid.q50_TOF
    
    score$ctrl.q75_TOF <- ctrl.q75_TOF$m
    resid.q75_TOF <- residuals(lm(score$q75_TOF~score$ctrl.q75_TOF))
    score$resid.q75_TOF <- resid.q75_TOF
    
    score$ctrl.mean_TOF <- ctrl.mean_TOF$m
    resid.mean_TOF <- residuals(lm(score$meanTOF~score$ctrl.mean_TOF))
    score$resid.mean_TOF <- resid.mean_TOF
    
    score$ctrl.q25_EXT <- ctrl.q25_EXT$m
    resid.q25_EXT <- residuals(lm(score$q25_EXT~score$ctrl.q25_EXT))
    score$resid.q25_EXT <- resid.q25_EXT
    
    score$ctrl.q50_EXT <- ctrl.q50_EXT$m
    resid.q50_EXT <- residuals(lm(score$medianEXT~score$ctrl.q50_EXT))
    score$resid.q50_EXT <- resid.q50_EXT
    
    score$ctrl.q75_EXT <- ctrl.q75_EXT$m
    resid.q75_EXT <- residuals(lm(score$q75_EXT~score$ctrl.q75_EXT))
    score$resid.q75_EXT <- resid.q75_EXT
    
    score$ctrl.mean_EXT <- ctrl.mean_EXT$m
    resid.mean_EXT <- residuals(lm(score$meanEXT~score$ctrl.mean_EXT))
    score$resid.mean_EXT <- resid.mean_EXT
    
    return(score)
}

#score.proc = dataset of fully processed score datasets
score.proc <-llply(score.pheno,function(x){addCtrl(x)})



bad<-score.plate
for(i in 1:length(bad))
{
    bad[[i]]<-paste0("p",bad[[i]])
    bad[[i]]<-get(bad[[i]])
}

#score.proc is list of totally processed score files
for (i in 1:length(score.proc))
{
    score.proc[[i]] <- removeWells(score.proc[[i]], bad[[i]])
}


date=Sys.Date()
date=as.character(format(date,format="%Y%m%d"))


score.proc[is.na(score.proc)]<-0



meltdf <- function(score)
{
    newscore<-data.frame(row=rep(score$row,each=1),col=rep(score$col,each=1),n=rep(score$n,each=1),f.L1=rep(score$f.L1,each=1),f.L2=rep(score$f.L2,each=1),f.L3=rep(score$f.L3,each=1),f.L4=rep(score$f.L4,each=1),f.ad=rep(score$f.ad,each=1))
    newscore<-melt(newscore,id.var=c("row","col","n"))
    return(newscore)
}

melted.score.proc<-llply(score.proc,function(x){meltdf(x)})



#Create complete data frame
score.info = llply(file.path(dir.data,score.filelist), function(x){info(x)})
for(i in 1:length(score.proc)){
    score.proc[[i]] = as.data.frame(cbind(score.info[[i]],score.proc[[i]]))
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
    
    
    proc<-score.proc[[i]]
    saveRDS(proc,file=file.path(dir.existing,"temp","proc.rds"))
    
    
    strainMean = mean(proc$n[!is.na(proc$strain)])
    strainSD = sd(proc$n[!is.na(proc$strain)])
    
    
    washMean = mean(proc$n[is.na(proc$strain)])
    washSD = sd(proc$n[is.na(proc$strain)])
    
    
    possibleContam = c()
    for(j in seq(1,nrow(proc),)){
        if(!is.na(proc[j,"strain"])){
            if(proc[j,"n"] > strainMean + (2*strainSD)){
                row = as.character(proc[j, "row"])
                col = as.numeric(proc[j, "col"])
                adjacentWash = proc[proc$row==row & proc$col==(col+1),"n"]
                if(adjacentWash > washMean + (2*washSD)){
                    possibleContam = append(possibleContam, paste0(row, col))
                }
            }
        }
    } 
    
    saveRDS(possibleContam,file=file.path(dir.existing,"temp","possibleContam.rds"))
    
    
    fileName = paste0(i,".csv")
    write.csv(proc,file = file.path(dir.existing,"results",fileName))
    
    
    split<-score.plate[[i]]
    saveRDS(split,file=file.path(dir.existing,"temp","split.rds"))
    
    plate<-paste0("p",split)
    contam<-get(plate)
    saveRDS(contam,file=file.path(dir.existing,"temp","contam.rds"))
    
    proc[is.na(proc)]<-0
    proc$norm.n<-round(proc$norm.n,0)
    proc$mean.normred<-round(proc$mean.normred,2)
    proc$median.normred<-round(proc$median.normred,2)
    proc$meanTOF<-round(proc$meanTOF,1)
    proc$meanEXT<-round(proc$meanEXT,1)
    
    plot.score.pop<-ggplot(proc)+geom_rect(aes(xmin=0,xmax=5,ymin=0,ymax=5,fill=norm.n))+facet_grid(row~col)+geom_text(aes(x=2.5,y=2.5,label=norm.n,colour="white"))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Score p",score.plate[[i]]," Population"))
    saveRDS(plot.score.pop,file=file.path(dir.existing,"temp","plot-score-pop.rds"))
    
    plot.score.red<-ggplot(proc)+geom_rect(aes(xmin=0,xmax=5,ymin=0,ymax=5),fill=NA)+facet_grid(row~col)+geom_text(aes(x=1,y=4,label=mean.normred))+geom_text(aes(x=4,y=1,label=median.normred))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Score p",score.plate[[i]]," Red Fluorescence"))
    saveRDS(plot.score.red,file=file.path(dir.existing,"temp","plot-score-red.rds"))
    
    plot.score.tofext<-ggplot(proc)+geom_rect(fill=NA,aes(xmin=0,xmax=5,ymin=0,ymax=5))+facet_grid(row~col)+geom_text(aes(x=1,y=4,label=meanTOF))+geom_text(aes(x=4,y=4,label=meanEXT))+geom_text(aes(x=1,y=1,label=medianTOF))+geom_text(aes(x=4,y=1,label=medianEXT))+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Score p",score.plate[[i]]," TOF/EXT"))
    saveRDS(plot.score.tofext,file=file.path(dir.existing,"temp","plot-score-tofext.rds"))
    
    melted.proc<-melted.score.proc[[i]]
    plot.score.dist<-ggplot(melted.proc,aes(as.factor(col),value,fill=variable))+geom_bar(aes(x=3),stat="identity",position="stack")+facet_grid(row~col)+presentation+theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+xlab("columns")+ylab("rows")+labs(title=paste0("Score p",score.plate[[i]]," Life-stage Distribution"))
    saveRDS(plot.score.dist,file=file.path(dir.existing,"temp","plot-score-dist.rds"))
    
    knit(file.report.score, file.path(dir.existing,"temp",paste0(date,'-',split,'-score.md')))
    markdownToHTML(paste0(date,'-',split,'-score.md'), file.path(dir.report,paste0(split,"_score.html")))
    
    unlink(file.path(dir.existing,"temp"), recursive = TRUE)
} 



filelist = dir(file.path(dir.existing,"results"), pattern = "[0-9]+.csv")

df = read.csv(file.path(dir.existing,"results",filelist[1]))

for(i in 2:length(filelist)){
    df = rbind(df, read.csv(file.path(dir.existing,"results",filelist[i])))
}

splitfp = strsplit(dir.data,"/")
dirName = splitfp[[1]][(length(splitfp[[1]]))]

details = strsplit(dirName,"_")[[1]][2]

name = paste0(details, "_complete.csv")

write.csv(df, file.path(dir.existing,"results",name))

file.remove(file.path(dir.existing,"results",filelist))

