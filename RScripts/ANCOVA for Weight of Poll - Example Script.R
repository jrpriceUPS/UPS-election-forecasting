#ANCOVA for Weight of Poll - Example Script
#06/30/2020

graphics.off() # This closes all of R's graphics windows.
#rm(list=ls())  # Clear all of R's memory!

#load the cleaned data
mydata=read.csv("Data/raw-polls_538_weekprior.csv")



#Order the data

newdata <- mydata[order(mydata$race_id),]

#create small data frame for reference
races = unique(newdata$race_id)
onlyUnique=(data.table::setDT(newdata)[,.SD[which.max(cand1_actual)],keyby=race_id])
actual = onlyUnique$cand1_actual

 refdataframe=data.frame(races,actual) 
 
#create the whichrace list. 
 
 
 whichrace <- vector(mode = "numeric", 0)
  for (i in unique(newdata$race_id)){
    subsetrace = subset (newdata, race_id==i)
    counter=nrow(subsetrace)
    whichrace = rlist::list.append(whichrace, counter)
  }
 
  
  whichrace=c(1,cumsum(whichrace))
 
  #create list to feed info with
  
  predictorsframe = newdata[,c("race_id","cand1_actual", "cand1_pct",
                               "delMode","transparency", "samplesize","LV")]
  
myDataFrame=predictorsframe
   
fileNameRootSim = "Simulations/Weight-Pollster-Bayes-ANCOVA-" 
fileNameRoot = "Markdown/Figures/Weight-Pollster-Bayes-ANCOVA-" 
graphFileType = "png" 

#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
#source("Jags-Ymet-Xnom1met1-MnormalHom.R")
source("RScripts/WeightPollANCOVA-V03.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( refFrame=refdataframe, datFrmPredictor=myDataFrame, pollName="cand1_pct" ,
                    actualName="actual", 
                    delModeName="delMode" , LVName="LV" , transparencyName="transparency", 
                    samplesizeName ="samplesize", raceIDName = "races", whichrace=whichrace,
                    numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in c("scoreSpread",   
                   "delMode[1]", "delMode[717]" ) ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , datFrm=myDataFrame , delModeName="delMode" , LVName="LV" , transparencyName="transparency", 
                        samplesizeName ="samplesize", 
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information: At this point just for delMode and sample size.
plotMCMC( mcmcCoda , datFrm=myDataFrame , 
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 

