# 2 way ANOVA (Pollster & Year)
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
#Load The data file 

figuresRoot="Markdown/Figures/YearPollsterNoramalHom-"
fileNameRoot = "Simulations/YearPollsterNormalHom-" 
graphFileType = "png" 

mydata=read.csv("Data/raw-polls_538_cleaned.csv")

earlyYear=2008
lateYear=2016
#Subset for year constraints given above.
mydata=(subset(mydata, mydata$year>=earlyYear))
mydata=(subset(mydata, mydata$year<=lateYear))



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

# Re-label and re-order the pollster factor:
myDataFrame$pollster = factor( myDataFrame$pollster)

# Specify the column names in the data file relevant to the analysis:
yName="bias" 


# x1 should be factor with fewer levels, to plot in single pane:
x1Name="pollster" 
x2Name="year" 

# Specify desired contrasts.
# Each main-effect contrast is a list of 2 vectors of level names, 
# a comparison value (typically 0.0), and a ROPE (which could be NULL):
x1contrasts = list( 
  list( c("SrvyUSA") , c("GrvsMrk") , compVal=0.0 , ROPE=c(-1,1) ) ,
  list( c("YouGov") , c("PblcPlP") , compVal=0.0 , ROPE=c(-1,1) ) 
)
x2contrasts = list( 
  list( c("2012") , c("2008") , compVal=0.0 , ROPE=c(-1,1) ) ,
  list( c("2016") , c("2012") , compVal=0.0 , ROPE=c(-1,1) ) 
 
)
# Each interaction contrast is a list of 2 lists of 2 vectors of level names, 
# a comparison value (typically 0.0), and a ROPE (which could be NULL)::
x1x2contrasts = list( 
  list( list( c("PblcPlP") , c("RRp-POR") ) ,
        list( c("2012") , c("2008") ) ,
        compVal=0.0 , ROPE=c(-1,1) ) ,
  list( list( c("YouGov") , c("SrvyUSA") ) ,
        list( c("2016") , c("2012") ) ,
        compVal=0.0 , ROPE=c(-1,1) ) 
) 

#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("RScripts/Jags-Ymet-Xnom2fac-MnormalHom.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , 
                    yName=yName , x1Name=x1Name , x2Name=x2Name ,
                    numSavedSteps=15000 , thinSteps=5 , 
                    saveName=fileNameRoot )
#x = na.omit(x)
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in c("b0","b1[1]","b2[1]","b1b2[1,1]","ySigma") ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=figuresRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , 
                        datFrm=myDataFrame , x1Name=x1Name , x2Name=x2Name ,
                        x1contrasts=x1contrasts , 
                        x2contrasts=x2contrasts , 
                        x1x2contrasts=x1x2contrasts ,
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , 
          datFrm=myDataFrame , yName=yName , x1Name=x1Name , x2Name=x2Name ,
          x1contrasts=x1contrasts , 
          x2contrasts=x2contrasts , 
          x1x2contrasts=x1x2contrasts ,
          saveName=figuresRoot , saveType=graphFileType )


#create small table to provide a pollster name key.
PollsterNames =  myDataFrame [, c("X","pollster","pollsterFull")]
#Get rid of repeats by choosing max X identifier within each pollster group
keyPollsterNames=(data.table::setDT(PollsterNames)[,.SD[which.max(X)],keyby=pollster])
#Get rid of X identifier
keyPollsterNames$X      <- NULL
#rename columns
keyPollsterNames=plyr::rename(keyPollsterNames, c("pollster"="Abbrevation", "pollsterFull"="Full Pollster Name"))

show(keyPollsterNames)
