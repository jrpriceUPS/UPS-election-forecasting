#example/practice script for "JAGS-PartyErrorPredictor.R"


graphics.off() # This closes all of R's graphics windows.
#rm(list=ls())  # Clear all of R's memory!

#load the cleaned data
mydata=read.csv("Data/raw-polls_538_cleaned.csv")


#create empty data frame, while maintaining all columns from the mydata structure
#Pick min. number of acceptable polls
minPolls=30
myDataMyPollsters=  mydata[0,]
for(myPollster in unique(mydata$pollster))
{
  
  #subset to a dataset with just each pollster
  subpoll=subset(mydata, pollster==myPollster)
  
  #if that subset has more than thirty entries, use it:
  if(nrow(subpoll)>minPolls){
    #combine each of these subsets together
    myDataMyPollsters=rbind(myDataMyPollsters, subpoll)
  }
}
myDataFrame = myDataMyPollsters
Undecided = 100 - myDataFrame$cand1_pct - myDataFrame$cand2_pct
myDataFrame$undecided = Undecided
myDataFrame$pollster = factor( myDataFrame$pollster)


#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("RScripts/JAGS-PartyErrorPredictor.R")

fileNameRoot = "Markdown/Figures/Jags-PartyErrorPredictorPractice-" 
fileNameRootSim= "Simulations/Jags-PartyErrorPredictorPractice-" 
graphFileType = "png"
# Generate the MCMC chain:
mcmcCodaPE = genMCMC( datFrm=myDataFrame , biasName = "demBias" , pollsterName = "pollster" , yearName = "year", undecidedName = "undecided",
                       numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )


parameterNames = varnames(mcmcCodaPE) 
show( parameterNames ) # show all parameter names, for reference

plotDiagnostics()


#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCodaPE , 
                        datFrm=myDataFrame ,
                        saveName=fileNameRootSim  )
show(summaryInfo)
#------------------------------------------------------------------------------- 

