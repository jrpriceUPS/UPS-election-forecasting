#Margin vs. Dem Error
#Haley Reed
#July 2nd 2020



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


#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("RScripts/Jags-2Factor-V03.R")


#Run through the Democratic Verision first 
fileNameRoot = "Markdown/Figures/Jags-2FactorPractice-PollsterV03-" 
fileNameRootSim= "Simulations/Jags-2FactorPractice-PollsterV03-"
graphFileType = "png"
# Generate the MCMC chain:
mcmcCodaV03 = genMCMC( datFrm=myDataFrame , biasName = "demBias" , pollsterName = "pollster" , yearName = "year",
                       numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )

mcmcCodaVD=mcmcCodaV03
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


summary(lR)
R.Rsquared=summary(lR)$r.squared
show(R.Rsquared)
summary(lD)
D.Rsquared=summary(lD)$r.squared
show(D.Rsquared)
#such weak r-squared value for republicans though that the linear model isn't quite representative

#---------------------------------------------------------------------------------------
#Is democratic polling bias= - Rebuplican polling bias?

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
#little realtionship




plot(dtDpollster$yearLabels,dtDpollster$Mode, col="blue", main="Pollster Biases", xlab="Election Year", 
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


