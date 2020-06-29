# #Script to test JAGS-2Factor.R
# 
# 
# graphics.off() # This closes all of R's graphics windows.
# rm(list=ls())  # Clear all of R's memory!
# 
# #load the cleaned data
# mydata=read.csv("Data/raw-polls_538_cleaned.csv")
# 
# 
# #create empty data frame, while maintaining all columns from the mydata structure
# #Pick min. number of acceptable polls
# minPolls=30
# myDataMyPollsters=  mydata[0,]
# for(myPollster in unique(mydata$pollster))
# {
#   
#   #subset to a dataset with just each pollster
#   subpoll=subset(mydata, pollster==myPollster)
#   
#   #if that subset has more than thirty entries, use it:
#   if(nrow(subpoll)>minPolls){
#     #combine each of these subsets together
#     myDataMyPollsters=rbind(myDataMyPollsters, subpoll)
#   }
# }
# myDataFrame = myDataMyPollsters
# myDataFrame$pollster = factor( myDataFrame$pollster)
# # Specify the column names in the data file relevant to the analysis:
# yName="bias" 
# x1Name="pollster" 
# x2Name= "year"
# 
# # Specify desired contrasts.
# # Each main-effect contrast is a list of 2 vectors of level names, 
# # a comparison value (typically 0.0), and a ROPE (which could be NULL):
# contrasts = list( 
#   list( c("RRp-POR") , c("SrvyUSA" ) , 
#         compVal=0.0 , ROPE=c(-1,1) ) ,
#   list( c( "YouGov" ) , c("ZIn-JZA" ) , 
#         compVal=0.0 , ROPE=c(-1,1) ) 
# )
# # Specify filename root and graphical format for saving output.
# # Otherwise specify as NULL or leave saveName and saveType arguments 
# # out of function calls.
# fileNameRoot = "Markdown/Figures/Jags-2FactorPractice-Pollster-" 
# fileNameRootSim= "Simulations/Jags-2FactorPractice-Pollster-"
# graphFileType = "png" 
# #------------------------------------------------------------------------------- 
# 
# # Load the relevant model into R's working memory:
# source("RScripts/Jags-2Factor.R")
# 
# 
# # Generate the MCMC chain:
# mcmcCoda = genMCMC( datFrm=myDataFrame , yName=yName , x1Name=x1Name , x2Name=x2Name,
#                     numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )
# #------------------------------------------------------------------------------- 
# # Display diagnostics of chain, for specified parameters:
# parameterNames = varnames(mcmcCoda) 
# show( parameterNames ) # show all parameter names, for reference
# for ( parName in c("a1[1]" ,    "a2[1]",
#                     "ySigma",  "nuY" ,
#                    "a1Sigma" ,"a2Sigma" ) ) {
#   diagMCMC( codaObject=mcmcCoda , parName=parName , 
#             saveName=fileNameRoot , saveType=graphFileType )
# }
# #------------------------------------------------------------------------------- 
# # Get summary statistics of chain:
# summaryInfo = smryMCMC( mcmcCoda , 
#                         datFrm=myDataFrame , x1Name=x1Name , x2Name=x2Name,
#                         contrasts=contrasts , 
#                         saveName=fileNameRootSim )
# show(summaryInfo)
# # Display posterior information:
# plotMCMC( mcmcCoda , 
#           datFrm=myDataFrame , yName=yName , x1Name=x1Name , x2Name=x2Name,
#           #contrasts=contrasts , 
#           saveName=fileNameRoot , saveType=graphFileType )
# 
 
 
 
 
 


#Version 3

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


#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("RScripts/Jags-2Factor-V03.R")

fileNameRoot = "Markdown/Figures/Jags-2FactorPractice-PollsterV03-" 
fileNameRootSim= "Simulations/Jags-2FactorPractice-PollsterV03-"
graphFileType = "png"
# Generate the MCMC chain:
mcmcCodaV03 = genMCMC( datFrm=myDataFrame , biasName = "bias" , pollsterName = "pollster" , yearName = "year",
                    numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )


parameterNames = varnames(mcmcCodaV03) 
show( parameterNames ) # show all parameter names, for reference

for ( parName in c("biasSpread",  "nuY" , "pollsterSpread" , "yearSpread" , "yearLean[1]", "pollsterBias[16,3]"  ) ) {
  diagMCMC( codaObject=mcmcCodaV03 , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCodaV03 , 
                        datFrm=myDataFrame ,  pollsterName = "pollster", yearName="year",
                        #contrasts=contrasts , 
                        saveName=fileNameRootSim )
show(summaryInfo)
plotMCMC( mcmcCodaV03, 
          datFrm=myDataFrame , biasName="bias" , pollsterName="pollster" , yearName="year",
          #contrasts=contrasts , 
          saveName=fileNameRoot , saveType=graphFileType )
