


graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Clear all of R's memory!

#load the cleaned data
#mydata=read.csv("Data/raw-polls_538_cleaned.csv")
mydata=read.csv("Data/raw-polls_538_weekprior.csv")

myPollsters = c("RRp-POR","SrvyUSA","ZIn-JZA","MrstCll")

#mydata = subset(mydata,pollster %in% pollstersList)





#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("RScripts/pollsterBiasModel.R")

dir.create("Markdown/Figures/Jags-pollsterBias/")
dir.create("Simulations/Jags-pollsterBias/")

fileNameRoot = "Markdown/Figures/Jags-pollsterBias/" 
fileNameRootSim= "Simulations/Jags-pollsterBias/" 
graphFileType = "png"

# Generate the MCMC chain:
demPrimarySim = runSimulation( mydata , party1 = "demBias", party2 = "repBias",
                               numSavedSteps=100000 , thinSteps=1 , saveName=fileNameRootSim) 
plotDiagnostics(demPrimarySim)

repPrimarySim = runSimulation( mydata , party1 = "repBias", party2 = "demBias",
                               numSavedSteps=100000 , thinSteps=1 , saveName=fileNameRootSim) 
plotDiagnostics(repPrimarySim)


#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
# summaryInfo = smryMCMC( mcmcCodaPE , 
#                         datFrm=myDataFrame ,
#                         saveName=fileNameRootSim  )
# show(summaryInfo)
#------------------------------------------------------------------------------- 
#Get some plotting done
plotParty1PosteriorPredictive( demPrimarySim, p1Name = "demBias",
                               datFrm=mydata, whichPollsters = myPollsters, saveName=fileNameRoot , saveType=graphFileType)

plotParty1PosteriorPredictive( repPrimarySim, p1Name = "repBias",
                               datFrm=mydata, whichPollsters = myPollsters, saveName=fileNameRoot , saveType=graphFileType)


plotMarginalDistributions( demPrimarySim, datFrm = mydata, whichPollsters = myPollsters,
                           p1Name = "demBias", saveName=fileNameRoot , saveType=graphFileType )

plotMarginalDistributions( repPrimarySim, datFrm = mydata,  whichPollsters = myPollsters,
                           p1Name = "repBias", saveName=fileNameRoot , saveType=graphFileType )