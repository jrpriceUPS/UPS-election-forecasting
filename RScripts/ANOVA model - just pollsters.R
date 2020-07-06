# Example for Jags-Ymet-Xnom1fac-MnormalHom.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
#------------------------------------------------------------------------------- 
# Load The data file 
#load the data in this case from 538
#load the cleaned data
mydata=read.csv("Data/raw-polls_538_cleaned.csv")


library(magrittr)
#limit to just pollsters with a certain number of polls on record.
myDataMyPollsters=  mydata[0,]
minPolls=10
for(myPollster in unique(mydata$pollster))
{
  #subset to a dataset with just each pollster
  subpoll=subset(mydata, pollster==myPollster)
  #if that subset has more than thirty entries, use it:
  if(nrow(subpoll)>minPolls){
    myDataMyPollsters=rbind(myDataMyPollsters, subpoll)
  }
}

mydata=myDataMyPollsters
#provide summary statistics for the bias.
summary(mydata$bias)



#set up for Bayesian analysis. 
myDataFrame = mydata

# Specify the column names in the data file relevant to the analysis:
yName="repBias" 
xName="pollster" 


# Specify desired contrasts.
# Each main-effect contrast is a list of 2 vectors of level names, 
# a comparison value (typically 0.0), and a ROPE (which could be NULL):
contrasts = list( 
  list( c("Public Policy Polling") , c("TCJ Research" ) , 
        compVal=0.0 , ROPE=c(-1,1) ) ,
  list( c( "YouGov" ) , c("Grove Insight" ) , 
        compVal=0.0 , ROPE=c(-1,1) ) 
)
# Specify filename root and graphical format for saving output.
# Otherwise specify as NULL or leave saveName and saveType arguments 
# out of function calls.
fileNameRoot = "Markdown/Figures/PollsterData1-NormalHom-" 
fileNameRootSim= "Simulations/PollsterData1-NormalHom-"
graphFileType = "png" 



#------------------------------------------------------------------------------- 

# Load the relevant model into R's working memory:
source("RScripts/Jags-Ymet-Xnom1fac-MnormalHom.R")

#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , yName=yName , xName=xName ,
                    numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )
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
                        #contrasts=contrasts , 
                        saveName=fileNameRootSim )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , 
          datFrm=myDataFrame , yName=yName , xName=xName ,
          contrasts=contrasts , 
          saveName=fileNameRoot , saveType=graphFileType )

FindMean=as.data.frame(summaryInfo)
FoundMean=FindMean[1, 1]

#------------------------------------------------------------------------------- 


#Helper for "Margin versus Party Error.Rmd"

demBiasPollsterSummary=as.data.frame(demBiasPollsterSummary)
demBiasPollsterSummary=demBiasPollsterSummary[2:83,3]
demBiasPollsterSummary=demBiasPollsterSummary-2.387878233
write.csv(demBiasPollsterSummary, "Simulations/demBiasPollsterSummary.csv")


repBiasPollsterSummary=as.data.frame(repBiasPollsterSummary)
repBiasPollsterSummary=repBiasPollsterSummary[2:83,3]
repBiasPollsterSummary=repBiasPollsterSummary-2.6987124
write.csv(repBiasPollsterSummary, "Simulations/repBiasPollsterSummary.csv")
