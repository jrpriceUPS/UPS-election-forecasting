# Example for Jags-Ymet-Xnom1fac-MnormalHom.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
#------------------------------------------------------------------------------- 
# Load The data file 



#load the data
#in this case from 538

setwd("~/Desktop/research/UPS-election-forecasting")
onlyRecentData <- read.csv("raw-polls_538.csv")
#To check for partisan bias, will subset data down to just the races where candidate 1 and candidate 2 are running on two different party tickets.
twoPartyTicket=subset(onlyRecentData,cand1_party!=cand2_party)
#defensive programming against NA values for bias (appears when one of the candidates is an independent) 
twoPartyTicket=subset(twoPartyTicket, !is.na(bias))
#plot bias towards republicans in polling as a dependent variable of the date of the election.

#provide summary statistics for the bias.
summary(twoPartyTicket$bias)
myDataFrame = twoPartyTicket

# Specify the column names in the data file relevant to the analysis:
yName="bias" 
xName="pollster" 


# Specify desired contrasts.
# Each main-effect contrast is a list of 2 vectors of level names, 
# a comparison value (typically 0.0), and a ROPE (which could be NULL):
contrasts = list( 
  list( c("Rasmussen Reports/Pulse Opinion Research" ) , c("Marist College" ) , compVal=0.0 , ROPE=c(-1.5,1.5) ) ,
  list( c(  "ABC News/The Washington Post" ) , c("SurveyUSA" ) , 
        compVal=0.0 , ROPE=c(-1.5,1.5) ) ,
  list( c( "Quinnipiac University" ) , c("Marist College" ) , 
        compVal=0.0 , ROPE=c(-1.5,1.5) ) ,
  list( c("Quinnipiac University" ) , c("Rasmussen Reports/Pulse Opinion Research" ) , compVal=0.0 , ROPE=c(-1.5,1.5) ) 
)
# Specify filename root and graphical format for saving output.
# Otherwise specify as NULL or leave saveName and saveType arguments 
# out of function calls.
fileNameRoot = "PollsterData-NormalHom-" 
graphFileType = "eps" 

# myDataFrame = read.csv( file="NonhomogVarData.csv" )
# yName="Y" 
# xName="Group" 
# contrasts = list( 
#   list( c("D") , c("A") , compVal=0.0 , ROPE=c(-1,1) ) ,
#   list( c("C") , c("B") , compVal=0.0 , ROPE=c(-1,1) ) 
# )
# fileNameRoot = "NonhomogVarData-NormalHom-" 
# graphFileType = "eps" 

# myDataFrame = read.csv( file="AnovaShrinkageData.csv" )
# yName="Y" 
# xName="Group" 
# contrasts = list( 
#   list( c("U") , c("A") , compVal=0.0 , ROPE=c(-1,1) ) ,
#   list( c("M") , c("A") , compVal=0.0 , ROPE=c(-1,1) ) ,
#   list( c("G") , c("A") , compVal=0.0 , ROPE=c(-1,1) ) 
# )
# fileNameRoot = "AnovaShrinkageData-FixSigmaB-" 
# graphFileType = "eps" 

#------------------------------------------------------------------------------- 
setwd("~/Downloads/DBDA2Eprograms 2")
# Load the relevant model into R's working memory:
source("Jags-Ymet-Xnom1fac-MnormalHom.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , yName=yName , xName=xName ,
                    numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in c("ySigma","b0","b[1]","aSigma") ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , 
                        datFrm=myDataFrame , xName=xName ,
                        contrasts=contrasts , 
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , 
          datFrm=myDataFrame , yName=yName , xName=xName ,
          contrasts=contrasts , 
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
