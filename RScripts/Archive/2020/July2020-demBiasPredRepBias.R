# A lot going on here - this was when we spent some time beginning to look at the relationship
# dem and repub error and undecided. There's some initial Fisherian exploration, then some 
# Bayesian sims that run dem + rep separately and compare found bias
#
# Finally, an ANCOVA model that tries to estimate rep bias using dem bias

#Margin vs. Dem Error
#Haley Reed
#July 2nd 2020


#load the cleaned data
mydata=read.csv("Data/raw-polls_538_cleaned.csv")
#is year rep. lean = -(dem. year lean) ?
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
lm(formula = myDataFrame$repBias ~ myDataFrame$demBias)

# Find out how many undecided voters there are in each poll
Undecided = 100 - myDataFrame$cand1_pct - myDataFrame$cand2_pct
myDataFrame$Undecided = Undecided
# What is the relationship between demBias and undecided Voters?
demUndecidedRelationship = lm( demBias ~ Undecided, data=myDataFrame)
summary(demUndecidedRelationship)
# What is the relationship between repBias and undecided voters?
repUndecidedRelationship = lm( repBias ~ Undecided, data=myDataFrame)
summary(repUndecidedRelationship)

plot(myDataFrame$Undecided,myDataFrame$demBias)
abline(demUndecidedRelationship)
plot(myDataFrame$Undecided,myDataFrame$repBias)
abline(repUndecidedRelationship)
#stronger relationship between Republican bias and Undecided voters
plot(myDataFrame$Undecided,myDataFrame$bias)
marginUndecidedRelationship = lm( bias ~ Undecided, data=myDataFrame)
summary(marginUndecidedRelationship)
abline(marginUndecidedRelationship)
#In general, the more undecided voters, the more polls underestimate the power of republicans.



#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("RScripts/Archive/2020/July2020-BiasPredictedByPollsterAndYear-NonGeneric.R")


#Run through the Democratic Verision first 
fileNameRoot = "Markdown/Figures/Jags-2FactorPractice-PollsterV03-" 
fileNameRootSim= "Simulations/Jags-2FactorPractice-PollsterV03-"
graphFileType = "png"
# Generate the MCMC chain:
mcmcCodaVD = genMCMC( datFrm=myDataFrame , biasName = "demBias" , pollsterName = "pollster" , yearName = "year",
                       numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )


parameterNames = varnames(mcmcCodaVD) 
summaryInfoD = smryMCMC( mcmcCodaVD , 
                        datFrm=myDataFrame ,  pollsterName = "pollster", yearName="year",
                        saveName=fileNameRootSim )
dtD <- summaryInfoD[4:10, 1:7]
plotYearPosterior( mcmcCodaVD, 
                   datFrm=myDataFrame , biasName="demBias" ,
                   saveName=fileNameRoot , saveType=graphFileType )
show(dtD)

fileNameRoot = "Markdown/Figures/Jags-2FactorPractice-PollsterV03-Rep" 
fileNameRootSim= "Simulations/Jags-2FactorPractice-PollsterV03-Rep"

mcmcCodaVR = genMCMC( datFrm=myDataFrame , biasName = "repBias" , pollsterName = "pollster" , yearName = "year",
                       numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )


parameterNames = varnames(mcmcCodaVR) 
summaryInfoR = smryMCMC( mcmcCodaVR , 
                        datFrm=myDataFrame ,  pollsterName = "pollster", yearName="year",
                        saveName=fileNameRootSim )
dtR <- summaryInfoR[4:10, 1:7]
plotYearPosterior( mcmcCodaVR, 
                   datFrm=myDataFrame , biasName="repBias" ,
                   saveName=fileNameRoot , saveType=graphFileType )
show(dtR)
yearLabels=c("2000","2004","2008","2012","2016")
dtR=as.data.frame(dtR[2:6,])
dtR=cbind(dtR,yearLabels)
dtD=as.data.frame(dtD[2:6,])
dtD=cbind(dtD,yearLabels)
show(dtD)
#points are matched by year
plot(dtD$Mode,dtR$Mode, main="Democratic Versus Rebulican Year Leans, 2000-2016", xlab="Democratic Year Lean", 
     ylab="Republican Year Lean")
RelationshipBetween = lm(dtR$Mode~dtD$Mode)
summary(RelationshipBetween)
abline(RelationshipBetween)
bestfit=lm(formula =  myDataFrame$demBias ~ myDataFrame$repBias)
abline(bestfit, col="red")
legend(-2.3, 0, legend=c("Best Fit for Bayesian Year Leans", "Best Fit for Fisherian Polls"), col=c("black", "red"), lty=1:2, cex=0.8)
#25.87% R-squared Value



plot(dtD$yearLabels,dtD$Mode, col="blue", main="Year Leans", xlab="Election Year", 
     ylab="Most Credible Error Value", ylim=c(-5, 1) , xlim=c(1996, 2020), xaxt = "n"
     )
# Define the position of tick marks
axis(side = 1, 
     at = yearLabels, 
     labels = yearLabels,
     tck=-.05)
points(dtR$yearLabels, dtR$Mode, col="red")
# Add a legend
legend(1996, 1, legend=c("Democratic Error", "Republican Error"),
       col=c("blue", "red"), lty=1:2, cex=0.8)



#Add lines of best fit?
dtR$yearLabels = as.numeric(yearLabels)
lR=lm(Mode ~ yearLabels, data=dtR)
abline(lR, col="red")



dtD$yearLabels = as.numeric(yearLabels)
lD=lm(Mode ~ yearLabels, data=dtD)
abline(lD, col="blue")

write.csv(dtD,"Simulations/dtD.csv")
write.csv(dtR,"Simulations/dtR.csv")

summary(lR)
R.Rsquared=summary(lR)$r.squared
show(R.Rsquared)
summary(lD)
D.Rsquared=summary(lD)$r.squared
show(D.Rsquared)
#such weak r-squared value for republicans though that the linear model isn't quite representative

#---------------------------------------------------------------------------------------
#Is democratic polling bias= - Republican polling bias?

mydata=read.csv("Data/raw-polls_538_cleaned.csv")
head(mydata[,33:34])

#so far it does not look like it

plotPollsterPosterior( mcmcCodaVD, 
                       datFrm=myDataFrame , biasName="demBias",
                       saveName=fileNameRoot , saveType=graphFileType )
plotPollsterPosterior( mcmcCodaVR, 
                       datFrm=myDataFrame , biasName="repBias",
                       saveName=fileNameRoot , saveType=graphFileType )

dtRpollster <- as.data.frame(summaryInfoR[10:94, 1:7])

dtDpollster <- as.data.frame(summaryInfoD[10:94, 1:7])

#points are matched by year

plot(dtDpollster$Mode,dtRpollster$Mode, main="Democratic Versus Rebulican Pollster Biases, 2000-2016", xlab="Democratic Pollster Bias", 
     ylab="Republican Pollster Bias")

RelationshipBetween1 = lm(dtRpollster$Mode~dtDpollster$Mode)
summary(RelationshipBetween1)
abline(RelationshipBetween1)
#adjusted r-squared is -.0046, not good...
#little relationship




# Example for Jags-Ymet-Xmet-Mrobust.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.

#------------------------------------------------------------------------------- 
# Load data file and specity column names of x (predictor) and y (predicted):
myData = myDataMyPollsters
xName = "demBias" ; yName = "repBias"
fileNameRoot = "Markdown/Figures/DoubleOneMetric-ErrorandBias-Jags-" 
fileNameRootSim = "Simulations/DoubleOneMetric-ErrorandBias-Jags-" 

xName = "demBias" ; yName = "repbias"
fileNameRoot = "Markdown/Figures/DoubleOneMetric-ErrorandBias-Jags-Margin" 
fileNameRootSim = "Simulations/DoubleOneMetric-ErrorandBias-Jags-Margin"

xName = "Undecided" ; yName = "repBias"
fileNameRoot = "Markdown/Figures/DoubleOneMetric-ErrorandBias-Jags-Undecided" 
fileNameRootSim = "Simulations/DoubleOneMetric-ErrorandBias-Jags-Undecided" 
myData=myDataFrame
#............................................................................

#............................................................................
graphFileType = "png" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("RScripts/Archive/DBDA2E-Code/Jags-Ymet-Xmet-Mrobust.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
#startTime = proc.time()
mcmcCoda = genMCMC( data=myData , xName=xName , yName=yName , 
                    numSavedSteps=20000 , saveName=fileNameRootSim )
#stopTime = proc.time()
#duration = stopTime - startTime
#show(duration)
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , 
                        compValBeta1=0.0 , ropeBeta1=c(-0.5,0.5) ,
                        saveName=fileNameRootSim )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , data=myData , xName=xName , yName=yName , 
         # compValBeta1=0.0 , ropeBeta1=c(-0.5,0.5) ,
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 



graphics.off() # This closes all of R's graphics windows.
myDataFrame= myData
# Specify the column names in the data file relevant to the analysis:
yName="repBias" 
xNomName="year" 
#xNomName="pollster"
xMetName="demBias"             # the covariate
# Specify desired contrasts of slopes.

# Specify filename root and graphical format for saving output.
# Otherwise specify as NULL or leave saveName and saveType arguments 
# out of function calls.

fileNameRoot="Markdown/Figures/yearLean-partyBias-ANCOVA-slopecontrast"
fileNameRootSim="Simulations/yearLean-partyBias-ANCOVA-slopecontrast"
#fileNameRoot = "Markdown/Figures/pollsterError-ANCOVA-" 
#fileNameRootSim = "Simulations/pollsterError-ANCOVA-" 
graphFileType = "png" 

#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
#source("Jags-Ymet-Xnom1met1-MnormalHom.R")
source("RScripts/Archive/2020/July2020-ANCOVA.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , 
                    yName=yName , xNomName=xNomName , xMetName=xMetName ,
                    numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )
#mcmcCodaPollsters=mcmcCoda
mcmcCodaYears=mcmcCoda
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}

#contrasts
contrasts = list( 
  list( c("TCJRsrc") , c("YouGov") , compVal=0.0 , ROPE=c(-1,1) ) ,
  list( c("RRp-POR") , c("QnnpcUn") , 
        compVal=0.0 , ROPE=c(-1,1) ) ,
  list( c("GrvsMrk") , c("SrvyUSA") , 
        compVal=0.0 , ROPE=c(-1,1) ) 
)


#year Contrasts
contrasts = list( 
  list( c("2000") , c("2004") , compVal=0.0 , ROPE=c(-1,1) ) ,
  list( c("2004") , c("2008") , 
        compVal=0.0 , ROPE=c(-.1,.1) ) ,
  list( c("2012") , c("2016") , 
        compVal=0.0 , ROPE=c(-.1,.1) ) 
)

#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , datFrm=myDataFrame , xNomName=xNomName , 
                        xMetName=xMetName , contrasts=contrasts , 
                        saveName=fileNameRootSim )
show(summaryInfo)

PollsterSumm=summaryInfo[1:34,1:7]

write.csv(PollsterSumm,"Simulations/PollsterSumm.csv")

YearLeanSumm=summaryInfo[1:10,1:7]
# 
write.csv(YearLeanSumm,"Simulations/YearLeanSumm.csv")

# Display posterior information:
plotMCMC( mcmcCoda , datFrm=myDataFrame , yName=yName , xNomName=xNomName , 
          xMetName=xMetName , contrasts=contrasts , 
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 

