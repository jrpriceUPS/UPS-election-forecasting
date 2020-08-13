



graphics.off() # This closes all of R's graphics windows.
#rm(list=ls())  # Clear all of R's memory!

#load the cleaned data
#mydata=read.csv("Data/raw-polls_538_cleaned.csv")
mydata=read.csv("Data/raw-polls_538_weekprior.csv")

pollstersList = c("RRp-POR","SrvyUSA","ZIn-JZA","MrstCll")

mydata = subset(mydata,pollster %in% pollstersList)





#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("RScripts/pollsterBiasModel.R")

fileNameRoot = "Markdown/Figures/Jags-pollsterBias--" 
fileNameRootSim= "Simulations/Jags-pollsterBias-" 
graphFileType = "png"

# Generate the MCMC chain:
sim1 = runSimulation( mydata , party1 = "demBias", party2 = "repBias",
                          numSavedSteps=50000 , thinSteps=1 , saveName=fileNameRootSim) 
sim2 = runSimulation( mydata , party1 = "repBias", party2 = "demBias",
                      numSavedSteps=50000 , thinSteps=1 , saveName=fileNameRootSim) 



#plotDiagnostics()


#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
# summaryInfo = smryMCMC( mcmcCodaPE , 
#                         datFrm=myDataFrame ,
#                         saveName=fileNameRootSim  )
# show(summaryInfo)
#------------------------------------------------------------------------------- 
#Get some plotting done
plotParty1PosteriorPredictive( sim1, p1Name = "demBias",
                            datFrm=mydata)
plotParty1PosteriorPredictive( sim2, p1Name = "repBias",
                               datFrm=mydata)
plotMarginalDistributions( sim1, saveName=fileNameRoot , saveType=graphFileType )
plotMarginalDistributions( sim2, saveName=fileNameRoot , saveType=graphFileType )
