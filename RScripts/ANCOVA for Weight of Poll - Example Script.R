#ANCOVA for Weight of Poll - Example Script
#06/30/2020

graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Clear all of R's memory!

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

#Add Score column to the dataframe, will use as dependent variable for model script.
score=100-myDataFrame$error
myDataFrame$score = score

fileNameRootSim = "Simulations/Weight-Pollster-Bayes-ANCOVA-" 
fileNameRoot = "Markdown/Figures/Weight-Pollster-Bayes-ANCOVA-" 
graphFileType = "png" 

#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
#source("Jags-Ymet-Xnom1met1-MnormalHom.R")
source("RScripts/WeightPollANCOVA.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm , scoreName="score" , daysuntilName="daysuntil", 
                    delModeName="delMode" , LVName="LV" , transparencyName="transparency", 
                    samplesizeName ="samplesize",
                    numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , datFrm=myDataFrame , xNomName=xNomName , 
                        xMetName=xMetName , contrasts=contrasts , 
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , datFrm=myDataFrame , yName=yName , xNomName=xNomName , 
          xMetName=xMetName , contrasts=contrasts , 
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 

