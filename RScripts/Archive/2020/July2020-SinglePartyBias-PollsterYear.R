# Version 3 of JAGS Pollster-Year 2 Factor ANOVA
# These are all kind of bleeding together for me so I hope taking a closer look next week will help remind me
# me of the development path. I think this was from when we ran the basic 2 Factor Pollster-Year analysis
# for democratic bias only.
#
# Runs!
 
 
 
 
 


#Version 3

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
myDataFrame$pollster = factor( myDataFrame$pollster)


#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("RScripts/Archive/2020/July2020-BiasPredictedByPollsterAndYear-NonGeneric.R")

fileNameRoot = "Markdown/Figures/Jags-2FactorPractice-PollsterV03-" 
fileNameRootSim= "Simulations/Jags-2FactorPractice-PollsterV03-"
graphFileType = "png"
# Generate the MCMC chain:
mcmcCodaV03 = genMCMC( datFrm=myDataFrame , biasName = "demBias" , pollsterName = "pollster" , yearName = "year",
                    numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )


parameterNames = varnames(mcmcCodaV03) 
show( parameterNames ) # show all parameter names, for reference

plotDiagnostics()

#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCodaV03 , 
                        datFrm=myDataFrame ,  pollsterName = "pollster", yearName="year",
                        saveName=fileNameRootSim )
show(summaryInfo)
#------------------------------------------------------------------------------- 
plotDiagnostics()
plotPosteriorPredictive( mcmcCodaV03, 
          datFrm=myDataFrame , biasName="demBias" ,
          saveName=fileNameRoot , saveType=graphFileType )
plotPollsterPosterior( mcmcCodaV03, 
                         datFrm=myDataFrame , biasName="demBias",
                         saveName=fileNameRoot , saveType=graphFileType )
plotYearPosterior( mcmcCodaV03, 
                   datFrm=myDataFrame , biasName="demBias" ,
                   saveName=fileNameRoot , saveType=graphFileType )


FoundMean = mean(mydata$demBias)
FoundMean=-FoundMean
show(FoundMean)
