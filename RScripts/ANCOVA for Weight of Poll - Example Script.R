#ANCOVA for Weight of Poll - Example Script
#06/30/2020

graphics.off() # This closes all of R's graphics windows.
#rm(list=ls())  # Clear all of R's memory!

#load the cleaned data
mydata=read.csv("Data/raw-polls_538_weekprior.csv")



#Order the data

newdata <- mydata[order(mydata$race_id),]
newdata2 = data.frame(newdata)

#create small data frame for reference
races = unique(newdata$race_id)
onlyUnique=(data.table::setDT(newdata2)[,.SD[which.max(cand1_actual)],keyby=race_id])
actual = onlyUnique$cand1_actual
raceFullName=unique(newdata$race)
refdataframe=data.frame(races,raceFullName,actual) 

for (i in 1:nrow(refdataframe)){
 
myFullName=refdataframe[i,2]

 myFullName= str_replace_all(myFullName, "_", " ")
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
  
  predictorsframe = newdata[,c("race_id","cand1_actual", "cand1_pct",
                               "delMode","transparency", "samplesize","LV")]
 
  
myDataFrame=predictorsframe


#Less Options for Del Mode is Helpful:
for (i in 1:nrow(myDataFrame)){
  myMode=myDataFrame[i,4]
  if(myMode=="IVR/Online"||myMode=="IVR/Online/Live"||myMode=="IVR/Online/Text"||myMode=="IVR/Online/Live/Text"||myMode=="	IVR/Online/Text"||myMode=="IVR/Text"){
    myDataFrame[i,4]="Online"
  }
  if(myMode=="Live*"||myMode=="Live/Text"||myMode=="Online/Live"||myMode=="IVR/Live"){
    myDataFrame[i,4]="Live"
  }
}

#remove mail and landline (3 total polls)
myDataFrame=myDataFrame[myDataFrame$delMode!="Landline",]
myDataFrame=myDataFrame[myDataFrame$delMode!="Mail",]

#which race indexing:
whichrace=match(unique(predictorsframe$race_id), predictorsframe$race_id)
whichrace=whichrace-1
whichrace=c(whichrace,nrow(newdata))


fileNameRootSim = "Simulations/Weight-Pollster-Bayes-ANCOVA-" 
fileNameRoot = "Markdown/Figures/Weight-Pollster-Bayes-ANCOVA-" 
graphFileType = "png" 

myDataFrame$samplesize =myDataFrame $ samplesize/1000
#myDataFrame$samplesize = myDataFrame$log(samplesize)

#MarginOfError = sqrt(.25/myDataFrame$samplesize)*100
#myDataFrame$samplesize =MarginOfError
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
#source("Jags-Ymet-Xnom1met1-MnormalHom.R")
source("RScripts/WeightPollANCOVA-V04.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( refFrame=refdataframe, datFrmPredictor=myDataFrame, pollName="cand1_pct" ,
                    actualName="actual", 
                    delModeName="delMode" , LVName="LV" , transparencyName="transparency", 
                    samplesizeName ="samplesize", raceIDName = "races", whichrace=whichrace,
                    numSavedSteps=11000 , thinSteps=10 , saveName=fileNameRootSim )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in c("actualSpread",   
                   "delModeImpact[1]", "samplesizeImpact" , "mu[1]") ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( codaSamples=mcmcCoda , datFrm=myDataFrame , delModeName="delMode" , LVName="LV" , transparencyName="transparency", 
                        samplesizeName ="samplesize", 
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information: At this point just for delMode and sample size.
#plotPosteriorPredictive( mcmcCoda , datFrm=myDataFrame , 
    #      saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
#plot the posterior predictive distrubtions

plotPosteriorPredictive(codaSample=mcmcCoda, refFrame=refdataframe, datFrmPredictor = myDataFrame, 
                        pollName = "cand1_pct", raceIDName="races", raceplots=1:10, whichrace=whichrace, saveName=fileNameRoot , 
                        saveType=graphFileType)


#plot the LV distribution
plotLVPosterior(codaSample=mcmcCoda,  datFrmPredictor = myDataFrame , saveName=fileNameRoot , 
                saveType=graphFileType)

#plot delModeImpact Posterior
plotdelModePosterior(codaSample=mcmcCoda,  datFrmPredictor = myDataFrame , saveName=fileNameRoot , 
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

plot(samplesizeImpact,LVImpact1, 
     ylab="LV Impact [1]", xlab="(Sample Size/100) Impact")
plot(samplesizeImpact,LVImpact2, 
     ylab="LV Impact [2]", xlab="(Sample Size/100) Impact")
plot(samplesizeImpact,transparencyImpact1, 
     ylab="Transparenct Impact [1]", xlab="(Sample Size/100) Impact")
plot(samplesizeImpact,transparencyImpact2, 
     ylab="Transparenct Impact [2]", xlab="(Sample Size/100) Impact")


#Look at the correlations between LV and delMode
mcmcMat = as.matrix(mcmcCoda,chains=TRUE)
delModeImpact2 = mcmcMat[,"delModeImpact[2]"]
delModeImpact7 = mcmcMat[,"delModeImpact[7]"]

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


plot(samplesizeImpact,samplesizeImpactLog, 
     ylab="Sample Size Impact - Log Transform", xlab="(Sample Size/100) Impact")
plot(samplesizeImpact,samplesizeImpactMOE,
     xlab="(Sample Size/100) Impact", ylab="Margin of Error Impact")
plot(samplesizeImpactLog,samplesizeImpactMOE,
     ylab="Margin of Error Impact",xlab="Sample Size Impact - Log Transform")

subdata=subset(mydata, year==2000)
lattice::densityplot(~bias, data = subdata, groups = LV, main = "2000 LV v. Non LV Comparsion")
lattice::densityplot(~bias, data = subdata, main = "2000 LV v. Non LV Comparsion")
