# Author: Jake Price
# Date: 8/25/20
#
# A script that models the margin bias of a pollster in a given year ("house effect")
#  * Overall year lean pulled from parent distribution of year leans
#  * Pollster bias in a given year comes from parent distribution for that pollster
#     - Pollster parent dist. comes from larger parent distribution of pollster centers



graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Clear all of R's memory!

#load the cleaned data
#mydata=read.csv("Data/raw-polls_538_cleaned.csv")
mydata=read.csv("Data/raw-polls_538_weekprior.csv")

myPollsters = c("RRp-POR","SrvyUSA","ZIn-JZA","MrstCll")

#mydata = subset(mydata,pollster %in% pollstersList)





#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("RScripts/marginBiasModel.R")

fileNameRoot = "Markdown/Figures/Jags-pollsterBias/" 
fileNameRootSim= "Simulations/Jags-pollsterBias/" 
graphFileType = "png"

dir.create(fileNameRoot)
dir.create(fileNameRootSim)

# Generate the MCMC chain:
marginSim = runSimulation( mydata ,
                               numSavedSteps=100000 , thinSteps=1 , saveName=fileNameRootSim) 
plotDiagnostics(marginSim)


#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
# summaryInfo = smryMCMC( mcmcCodaPE , 
#                         datFrm=myDataFrame ,
#                         saveName=fileNameRootSim  )
# show(summaryInfo)
#------------------------------------------------------------------------------- 
#Get some plotting done
plotPosteriorPredictive( marginSim, 
                               datFrm=mydata, whichPollsters = myPollsters, saveName=fileNameRoot , saveType=graphFileType)


plotMarginalDistributions( marginSim, datFrm = mydata, whichPollsters = myPollsters,
                            saveName=fileNameRoot , saveType=graphFileType )
