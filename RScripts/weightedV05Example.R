

#ANCOVA for Weight of Poll - Example Script
#06/30/2020

graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Clear all of R's memory!

#load the cleaned data
#mydata=read.csv("Data/raw-polls_538_cleaned.csv")
mydata=read.csv("Data/raw-polls_538_weekprior.csv")
#remove mail and landline (3 total polls)
mydata=mydata[mydata$delMode!="Landline",]
mydata=mydata[mydata$delMode!="Mail",]

#limit to just US overall, no more consideration of state elections
#mydata=mydata[mydata$location=="US",]


#Order the data

newdata <- mydata[order(mydata$race_id),]
newdata2 = data.frame(newdata)

#create small data frame for reference
races = unique(newdata$race_id)
onlyUnique=(data.table::setDT(newdata2)[,.SD[which.max(margin_actual)],keyby=race_id])
actual = onlyUnique$margin_actual
raceFullName=unique(newdata$race)
refdataframe=data.frame(races,raceFullName,actual) 

for (i in 1:nrow(refdataframe)){
  
  myFullName=refdataframe[i,2]
  
  myFullName= stringr::str_replace_all(myFullName, "_", " ")
  refdataframe[i,2]=myFullName
}

# refdataframe=data.frame(races,raceFullName,actual) 

#create the whichrace list. 


# whichrace <- vector(mode = "numeric", 0)
#  for (i in unique(newdata$race_id)){
#    subsetrace = subset (newdata, race_id==i)
#    counter=nrow(subsetrace)
#    whichrace = rlist::list.append(whichrace, counter)
#  }


#create list to feed info with

predictorsframe = newdata[,c("race_id","margin_actual", "margin_poll",
                             "delMode","transparency", "samplesize","LV", "pollster","year","bias")]




predictorsframe=cbind(predictorsframe, IVR="FALSE")
predictorsframe=cbind(predictorsframe, online="FALSE")
predictorsframe=cbind(predictorsframe, live="FALSE")
predictorsframe=cbind(predictorsframe, text="FALSE")

myDataFrame=predictorsframe
#Less Options for Del Mode is Helpful:
for (i in 1:nrow(myDataFrame)){
  myMode=myDataFrame[i,4]
  if(myMode=="IVR/Online"||myMode=="IVR/Online/Live"||myMode=="IVR/Online/Text"||myMode=="IVR/Online/Live/Text"||myMode=="IVR/Online/Text"||myMode=="IVR/Text"||myMode=="Online/Live"){
    myDataFrame[i,9]="TRUE"
  }
  if(myMode=="Live*"||myMode=="Live/Text"||myMode=="Online/Live"||myMode=="IVR/Live"||myMode=="IVR/Online/Live"||myMode=="IVR/Online/Text"||myMode=="IVR/Online/Live/Text"){
    myDataFrame[i,10]="TRUE"
  }
  if(myMode=="IVR"||myMode=="IVR/Online"||myMode=="IVR/Online/Live"||myMode=="IVR/Online/Text"||myMode=="IVR/Online/Live/Text"||myMode=="	IVR/Online/Text"||myMode=="IVR/Text"){
    myDataFrame[i,8]="TRUE"
  }
  if(myMode=="IVR/Online/Text"||myMode=="IVR/Text"||myMode=="IVR/Online/Live/Text"||myMode=="Live/Text"){
    myDataFrame[i,11]="TRUE"
  }
}



predictorsframe=myDataFrame
#which race indexing:
whichrace=match(unique(predictorsframe$race_id), predictorsframe$race_id)
whichrace=whichrace-1
whichrace=c(whichrace,nrow(newdata))


fileNameRootSim = "Simulations/WeightedandCorrected" 
fileNameRoot = "Markdown/Figures/WeightedModels/WeightedandCorrected/V05" 
graphFileType = "png" 

myDataFrame$samplesize =myDataFrame $ samplesize/1000
#myDataFrame$samplesize = myDataFrame$log(samplesize)

#MarginOfError = sqrt(.25/myDataFrame$samplesize)*100
#myDataFrame$samplesize =MarginOfError
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
#source("Jags-Ymet-Xnom1met1-MnormalHom.R")
source("RScripts/WeightPollANCOVA-V05.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( refFrame=refdataframe, datFrmPredictor=myDataFrame, pollName="margin_poll" ,
                    actualName="actual", 
                    IVRName="IVR", onlineName="online", liveName="live", textName="text",
                    LVName="LV" , transparencyName="transparency", 
                    samplesizeName ="samplesize", raceIDName = "races", whichrace=whichrace,
                    numSavedSteps=21000 , thinSteps=10 , saveName=fileNameRootSim )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in c("actualSpread",   
                   "IVRImpact", "samplesizeImpact" , "mu[1]","LVImpact") ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( codaSamples=mcmcCoda , datFrm=myDataFrame , LVName="LV" , transparencyName="transparency", 
                        samplesizeName ="samplesize",   IVRName="IVR" ,onlineName="online" ,liveName="live" ,textName="text" ,
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information: At this point just for delMode and sample size.
#plotPosteriorPredictive( mcmcCoda , datFrm=myDataFrame , 
#      saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
#plot the posterior predictive distrubtions

plotPosteriorPredictive(codaSample=mcmcCoda, refFrame=refdataframe, datFrmPredictor = myDataFrame, 
                        pollName = "margin_poll", raceIDName="races", raceplots=1:10, whichrace=whichrace, saveName=fileNameRoot , 
                        saveType=graphFileType)

meanPosterior(codaSample=mcmcCoda, refFrame=refdataframe, datFrmPredictor = myDataFrame, 
              pollName = "margin_poll", raceIDName="races", raceplots=1:10, whichrace=whichrace, saveName=fileNameRoot , 
              saveType=graphFileType)
#plot the LV distribution
plotLVPosterior(codaSample=mcmcCoda,  datFrmPredictor = myDataFrame , saveName=fileNameRoot , 
                saveType=graphFileType)

#plot delModeImpact Posteriors
plotIVRPosterior(codaSample=mcmcCoda,  datFrmPredictor = myDataFrame , saveName=fileNameRoot , 
                 saveType=graphFileType)
plotonlinePosterior(codaSample=mcmcCoda,  datFrmPredictor = myDataFrame , saveName=fileNameRoot , 
                    saveType=graphFileType)
plotlivePosterior(codaSample=mcmcCoda,  datFrmPredictor = myDataFrame , saveName=fileNameRoot , 
                  saveType=graphFileType)
plottextPosterior(codaSample=mcmcCoda,  datFrmPredictor = myDataFrame , saveName=fileNameRoot , 
                  saveType=graphFileType)

#plot Transparenct Posterior
plotTransparencyPosterior(codaSample=mcmcCoda,  datFrmPredictor = myDataFrame , saveName=fileNameRoot , 
                          saveType=graphFileType)

#plot the samplesize 
plotSampleSizePosterior(mcmcCoda, datFrm=myDataFrame,  saveName=fileNameRoot , 
                        saveType=graphFileType, title="Sample Size Impact ")

#plot SampleSize Prior 

agammaShRa = unlist( gammaShRaFromModeSD( mode=sd(actual)/2 , sd=2*sd(actual) ) )

#X <- rgamma( shape=1.2832, scale=0.0624, n=2750)
X <- rgamma( shape=1.2832, scale=0.0624, n=11000)
mean(X)
mode(X)

hist(X,prob=T,main='Gamma Prior', breaks=47)
#lines(density(X),col='red',lwd=2)


#plot the different samplesize impacts against each other:
mcmcMat = as.matrix(mcmcCoda,chains=TRUE)
samplesizeImpact= mcmcMat[,"samplesizeImpact"]
LVImpact1 = mcmcMat[,"LVImpact[1]"]
LVImpact2 = mcmcMat[,"LVImpact[2]"]
transparencyImpact1 = mcmcMat[,"transparencyImpact[1]"]
transparencyImpact2 = mcmcMat[,"transparencyImpact[2]"]

IVRImpact1=mcmcMat[,"IVRImpact[1]"]
IVRImpact2=mcmcMat[,"IVRImpact[2]"]
onlineImpact1=mcmcMat[,"onlineImpact[1]"]
onlineImpact2=mcmcMat[,"onlineImpact[2]"]
liveImpact1=mcmcMat[,"liveImpact[1]"]
liveImpact2=mcmcMat[,"liveImpact[2]"]
textImpact1=mcmcMat[,"textImpact[1]"]
textImpact2=mcmcMat[,"textImpact[2]"]


plot(IVRImpact1,IVRImpact2)
plot(IVRImpact1,transparencyImpact1)
plot(IVRImpact1,onlineImpact1)
plot(textImpact1,textImpact2)
cor(textImpact1,textImpact2)
cor(IVRImpact1,IVRImpact2)
cor(IVRImpact1,transparencyImpact1)
cor(IVRImpact2,onlineImpact2)

plot(samplesizeImpact,LVImpact1, 
     ylab="LV Impact [1]", xlab="(Sample Size/100) Impact")
plot(samplesizeImpact,LVImpact2, 
     ylab="LV Impact [2]", xlab="(Sample Size/100) Impact")
plot(samplesizeImpact,transparencyImpact1, 
     ylab="Transparenct Impact [1]", xlab="(Sample Size/100) Impact")
plot(samplesizeImpact,transparencyImpact2, 
     ylab="Transparenct Impact [2]", xlab="(Sample Size/100) Impact")


#Look at the correlations between LV and delMode


plot(LVImpact1,delModeImpact2, xlab="Non LV Model Impact", ylab="IVR Combo Impact") 
plot(LVImpact2,delModeImpact7, xlab="LV Model Impact", ylab="Online Impact") 

plot(LVImpact1, LVImpact2, xlab="Non LV Model Impact", ylab="LV Model Impact" )
abline(a=0, b=1)
# mcmcMat = as.matrix(mcmcCoda,chains=TRUE)
# samplesizeImpactLog= mcmcMat[,"samplesizeImpact"]
# 
# 
# mcmcMat = as.matrix(mcmcCoda,chains=TRUE)
# samplesizeImpactMOE= mcmcMat[,"samplesizeImpact"]



subdata=subset(mydata, year==2000)
lattice::densityplot(~bias, data = subdata, groups = LV, main = "2000 LV v. Non LV Comparsion")
lattice::densityplot(~bias, data = subdata, main = "2000 LV v. Non LV Comparsion")
