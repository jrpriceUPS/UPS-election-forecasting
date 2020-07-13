# Tester script for posterior plot with Dem. Error, Rep. Error, and Marginal Bias.

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
source("RScripts/3PartPosterior.R")

fileNameRoot = "Markdown/Figures/3PartPosterior-" 
fileNameRootSim= "Simulations/3PartPosterior-"
graphFileType = "png"
# Generate the MCMC chain:
mcmcCodabias = genMCMC( datFrm=myDataFrame , biasName = "bias" , pollsterName = "pollster" , yearName = "year",
                        numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )
mcmcCodademBias = genMCMC( datFrm=myDataFrame , biasName = "demBias" , pollsterName = "pollster" , yearName = "year",
                           numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )
mcmcCodarepBias = genMCMC( datFrm=myDataFrame , biasName = "repBias" , pollsterName = "pollster" , yearName = "year",
                           numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )

#------------------------------------------------------------------------------- 

plotPosteriorPredictive( codaSamplesbias=mcmcCodabias, codaSamplesdemBias = mcmcCodademBias, codaSamplesrepBias = mcmcCodarepBias,
                         datFrm=myDataFrame , biasName="demBias" ,
                         saveName=fileNameRoot , saveType=graphFileType )








############################################################################################################





#Lets Create a posterior plot for the presentation!

graphics.off() # This closes all of R's graphics windows.
#rm(list=ls())  # Clear all of R's memory!

#load the cleaned data
mydata=read.csv("Data/raw-polls_538_cleaned.csv")


#create empty data frame, while maintaining all columns from the mydata structure

myDataMyPollsters=  mydata[0,]
for(myPollster in unique(mydata$pollster))
{
  
  #subset to a dataset with just each pollster
  subpoll=subset(mydata, pollster==myPollster)
  
 
  #Add preferred pollsters
  if(subpoll$pollster=="YouGov"){
    myDataMyPollsters=rbind(myDataMyPollsters, subpoll)
  }
  if(subpoll$pollster=="RRp-POR"){
    myDataMyPollsters=rbind(myDataMyPollsters, subpoll)
  }
  if(subpoll$pollster=="SrvyUSA"){
    myDataMyPollsters=rbind(myDataMyPollsters, subpoll)
  }
  if(subpoll$pollster=="MrstCll"){
    myDataMyPollsters=rbind(myDataMyPollsters, subpoll)
  }
}


#limit to just 2008,2012, and 2016
YearSub = myDataMyPollsters[0,]
for(myYear in unique(myDataMyPollsters$year)){
  subYear=subset(myDataMyPollsters, year==myYear)
  if(subYear$year>2005){
    YearSub=rbind(YearSub, subYear)
  }
}


myDataFrame = YearSub
myDataFrame$pollster = factor( myDataFrame$pollster)


#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("RScripts/3PartPosterior.R")

#fileNameRoot = "Markdown/Figures/PosteriorsForPresentation-" 
#fileNameRootSim= "Simulations/PosteriorsForPresentation-"
fileNameRoot = "Markdown/Figures/PosteriorsForPresentation-4Pollsters" 
fileNameRootSim= "Simulations/PosteriorsForPresentation-4Pollster"
graphFileType = "png"
# Generate the MCMC chain:
mcmcCodabias = genMCMC( datFrm=myDataFrame , biasName = "bias" , pollsterName = "pollster" , yearName = "year",
                        numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )
mcmcCodademBias = genMCMC( datFrm=myDataFrame , biasName = "demBias" , pollsterName = "pollster" , yearName = "year",
                           numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )
mcmcCodarepBias = genMCMC( datFrm=myDataFrame , biasName = "repBias" , pollsterName = "pollster" , yearName = "year",
                           numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )

#------------------------------------------------------------------------------- 

plotPosteriorPredictiveV02( codaSamplesbias=mcmcCodabias, codaSamplesdemBias = mcmcCodademBias, codaSamplesrepBias = mcmcCodarepBias,
                         datFrm=myDataFrame , biasName="demBias" ,pollsterName="pollsterFull",
                         saveName=fileNameRoot , saveType=graphFileType )



source("RScripts/Jags-2Factor-V03.R")
fileNameRoot = "Markdown/Figures/HistogramsForPresentation-4Pollsters" 
fileNameRootSim= "Simulations/HistogramssForPresentation-4Pollster"
plotPollsterPosterior( mcmcCodabias, 
                       datFrm=myDataFrame , biasName="bias",
                       saveName=fileNameRoot , saveType=graphFileType )
plotYearPosterior( mcmcCodabias, 
                   datFrm=myDataFrame , biasName="bias" ,
                   saveName=fileNameRoot , saveType=graphFileType )

