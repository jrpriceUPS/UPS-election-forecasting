# 2 way ANOVA (Pollster & Year)
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
#Load The data file 

fileNameRootFig="Markdown/Figures/YearPollsterNoramalHom-"
fileNameRoot = "Simulations/YearPollsterNormalHom-" 
graphFileType = "png" 

mydata=read.csv("Data/raw-polls_538.csv")
#Limit to just president
preferredType="Pres-G"
mydata=subset(mydata, type_simple==preferredType)
#use lubridate to change dates from character to date data type for functionality. 
mydata$electiondate = lubridate::mdy(mydata$electiondate)
mydata$polldate = lubridate::mdy(mydata$polldate)

earlyYear=2010
lateYear=2020
#Subset for year constraints given above.
mydata=(subset(mydata, mydata$year>=earlyYear))
mydata=(subset(mydata, mydata$year<=lateYear))

#subset to make sure it is a democratic v. republican race
mydata=subset(mydata,cand1_party=="DEM")
mydata=subset(mydata,cand2_party=="REP")


mydata$year = as.factor(mydata$year)

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
  list( c("SurveyUSA") , c("Marist College") , compVal=0.0 , ROPE=c(-1,1) ) ,
  list( c("YouGov") , c("Grove Insight") , compVal=0.0 , ROPE=c(-1,1) ) 
)
x2contrasts = list( 
  list( c("2012") , c("2020") , compVal=0.0 , ROPE=c(-1,1) ) ,
  list( c("2016") , c("2020") , compVal=0.0 , ROPE=c(-1,1) ) 
 
)
# Each interaction contrast is a list of 2 lists of 2 vectors of level names, 
# a comparison value (typically 0.0), and a ROPE (which could be NULL)::
x1x2contrasts = list( 
  list( list( c("Zogby Interactive/JZ Analytics") , c("Quinnipiac University") ) ,
        list( c("2012") , c("2020") ) ,
        compVal=0.0 , ROPE=c(-1,1) ) ,
  list( list( c("YouGov") , c("Grove Insight") ) ,
        list( c("2016") , c("2020") ) ,
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
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in c("b0","b1[1]","b2[1]","b1b2[1,1]","ySigma") ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
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
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
# Other specific comparisons of cells:
if ( fileNameRoot == "YearPollsterNormalHom-" ) {
  # THIS x1level minus THAT x1level at AT x2level:
  THISx1 = "Full"
  THATx1 = "Assis"
  ATx2 = "CHEM"
  THISidx = which(levels(myDataFrame[,x1Name])==THISx1)
  THATidx = which(levels(myDataFrame[,x1Name])==THATx1)
  ATidx   = which(levels(myDataFrame[,x2Name])==ATx2)
  openGraph(height=4,width=4)
  compInfo = plotPost( 
    as.matrix(mcmcCoda)[,paste("m[",THISidx,",",ATidx,"]",sep="")] -
      as.matrix(mcmcCoda)[,paste("m[",THATidx,",",ATidx,"]",sep="")] , 
    main=paste(THISx1,"-",THATx1,"@",ATx2) , 
    xlab=paste("Difference in",yName) , 
    compVal=0 ,ROPE=c(-1000,1000) )
  show(compInfo)
  saveGraph(file=paste(fileNameRoot,THISx1,"-",THATx1,"At",ATx2,sep=""),
            type=graphFileType)
  # THIS x1level minus THAT x1level at AT x2level:
  THISx1 = "Full"
  THATx1 = "Assis"
  ATx2 = "PSY"
  THISidx = which(levels(myDataFrame[,x1Name])==THISx1)
  THATidx = which(levels(myDataFrame[,x1Name])==THATx1)
  ATidx   = which(levels(myDataFrame[,x2Name])==ATx2)
  openGraph(height=4,width=4)
  compInfo = plotPost( 
    as.matrix(mcmcCoda)[,paste("m[",THISidx,",",ATidx,"]",sep="")] -
      as.matrix(mcmcCoda)[,paste("m[",THATidx,",",ATidx,"]",sep="")] , 
    main=paste(THISx1,"-",THATx1,"@",ATx2) , 
    xlab=paste("Difference in",yName) , 
    compVal=0 ,ROPE=c(-1000,1000) )
  show(compInfo)
  saveGraph(file=paste(fileNameRoot,THISx1,"-",THATx1,"At",ATx2,sep=""),
            type=graphFileType)
  # THIS x2level minus THAT x2level at AT x1level:
  THISx2 = "PSY"
  THATx2 = "ENG"
  ATx1 = "Full"
  THISidx = which(levels(myDataFrame[,x2Name])==THISx2)
  THATidx = which(levels(myDataFrame[,x2Name])==THATx2)
  ATidx   = which(levels(myDataFrame[,x1Name])==ATx1)
  openGraph(height=4,width=4)
  compInfo = plotPost( 
    as.matrix(mcmcCoda)[,paste("m[",ATidx,",",THISidx,"]",sep="")] -
      as.matrix(mcmcCoda)[,paste("m[",ATidx,",",THATidx,"]",sep="")] , 
    main=paste(THISx2,"-",THATx2,"@",ATx1) , 
    xlab=paste("Difference in",yName) , 
    compVal=0 ,ROPE=c(-1000,1000) )
  show(compInfo)
  saveGraph(file=paste(fileNameRoot,THISx2,"-",THATx2,"At",ATx1,sep=""),
            type=graphFileType)
}