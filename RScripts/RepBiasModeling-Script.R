# Calls July2020-RepBiasModeling.R to run simulation
# that estimates republican bias using dem bias and undecideds
# can instead call RScripts/Archive/2020/July2020-RepBiasModeling-demResponseYear.R
# to allow the response to dem bias to differ by year
#
# Runs!
#
# Notes by Jake 8/10/20

#example/practice script for "JAGS-PartyErrorPredictor.R"



graphics.off() # This closes all of R's graphics windows.
#rm(list=ls())  # Clear all of R's memory!

#load the cleaned data
#mydata=read.csv("Data/raw-polls_538_cleaned.csv")
mydata=read.csv("Data/raw-polls_538_weekprior.csv")


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

myDataFrame = myDataMyPollsters

myDataFrame$pollster = factor( myDataFrame$pollster)





#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
#source("RScripts/JAGS-PartyErrorPredictor.R")
source("RScripts/RepBiasModeling.R")
fileNameRoot = "Markdown/Figures/Jags-PartyErrorPredictorPractice-" 
fileNameRootSim= "Simulations/Jags-PartyErrorPredictorPractice-" 

# fileNameRoot = "Markdown/Figures/Jags-PartyErrorPredictorPractice-" 
# fileNameRootSim= "Simulations/Jags-PartyErrorPredictorPractice-" 
graphFileType = "png"
# Generate the MCMC chain:
mcmcCodaPE = genMCMC( datFrm=myDataFrame , dembiasName = "demBias" , pollsterName = "pollster" , yearName = "year", undecidedName = "undecided",
                       numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )


parameterNames = varnames(mcmcCodaPE) 
show( parameterNames ) # show all parameter names, for reference

#plotDiagnostics()


#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCodaPE , 
                        datFrm=myDataFrame ,
                        saveName=fileNameRootSim  )
show(summaryInfo)
#------------------------------------------------------------------------------- 
#Get some plotting done
plotPosteriorPredictiveDem( mcmcCodaPE, 
                         datFrm=myDataFrame , dembiasName="demBias" ,
                         saveName=fileNameRoot , saveType=graphFileType )
plotMCMCrep( mcmcCodaPE, saveName=fileNameRoot , saveType=graphFileType )
