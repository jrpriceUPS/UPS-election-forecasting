

polls08 = read.csv("https://raw.githubusercontent.com/TheEconomist/us-potus-model/master/data/all_polls_2008.csv")
polls12 = read.csv("https://raw.githubusercontent.com/TheEconomist/us-potus-model/master/data/all_polls_2012.csv")
polls16 = read.csv("https://raw.githubusercontent.com/TheEconomist/us-potus-model/master/data/all_polls.csv")
polls20 = read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQ56fySJKLL18Lipu1_i3ID9JE06voJEz2EXm6JW4Vh11zmndyTwejMavuNntzIWLY0RyhA1UsVEen0/pub?gid=0&single=true&output=csv")

elecResults = read.csv("Data/1976-2016-president.csv")
elecResults = subset(elecResults, year > 2007)

polls = data.frame(year = c(rep(2008, nrow(polls08)),rep(2012,nrow(polls12)),rep(2016,nrow(polls16)),rep(2020,nrow(polls20))),
                   state = c(polls08$state,polls12$state,polls16$state,polls20$state),
                   pollster = c(polls08$pollster,polls12$pollster,polls16$pollster,polls20$pollster),
                   dem = c(polls08$obama,polls12$obama,polls16$clinton,polls20$biden),
                   rep = c(polls08$mccain,polls12$romney,polls16$trump,polls20$trump),
                   mode = c(polls08$mode,polls12$mode,polls16$mode,polls20$mode),
                   n = c(polls08$number.of.observations,polls12$number.of.observations,polls16$number.of.observations,polls20$number.of.observations),
                   start.date = c(lubridate::mdy(polls08$start.date),lubridate::mdy(polls12$start.date),lubridate::ymd(polls16$start.date),lubridate::mdy(polls20$start.date)),
                   end.date = c(lubridate::mdy(polls08$end.date),lubridate::mdy(polls12$end.date),lubridate::ymd(polls16$end.date),lubridate::mdy(polls20$end.date)),
                   elec.date = c(rep("2008-11-04", nrow(polls08)),rep("2012-11-06",nrow(polls12)),rep("2016-11-08",nrow(polls16)),rep("2020-11-03",nrow(polls20))),
                   pop = c(polls08$population,polls12$population,polls16$population,polls20$population)
                   )

polls = transform(polls, race_id = paste(year,as.numeric(as.factor(polls$state)),sep =""))
polls$pollster = as.factor(polls$pollster)
polls$year = as.factor(polls$year)
# Make copy of full pollster name for reference
polls$pollsterFullName = polls$pollster
#Abbreviate to limit Pollster names to a reasonable length.
levels(polls$pollster) = abbreviate(levels(polls$pollster), minlength = 7, use.classes = TRUE,
                                               dot = FALSE, strict = FALSE)

polls = transform(polls, daysuntil = as.numeric(difftime(elec.date,end.date, units = "days")))
polls = transform(polls, margin = dem - rep)
polls = transform(polls, undecided = 100 - dem - rep)
polls = transform(polls,actualDem = NA)
polls = transform(polls,actualRep = NA)

years = unique(polls$year)
states = unique(polls$state)[1:51]

demList = c("Obama, Barack H.","Obama, Barack H.", "Clinton, Hillary")
repList = c("McCain, John","Romney, Mitt","Trump, Donald J.")
electionDay = c("11/04/08","11/06/12","11/08/2016")

# enter in real results
for (myIndex in 1:3){
  myYear = years[myIndex]
  myDem = demList[myIndex]
  myRep = repList[myIndex]
  
  demPopVote = 0
  repPopVote = 0
  totalPopVote = 0
  for (myState in states){
    yearStateResults = subset(elecResults, year == myYear & state_po == myState)
    demVotes = sum(yearStateResults[yearStateResults$candidate == demList[myIndex],"candidatevotes"])
    repVotes = sum(yearStateResults[yearStateResults$candidate == repList[myIndex],"candidatevotes"])
    totalVotes = yearStateResults[yearStateResults$candidate == demList[myIndex],"totalvotes"][1]
    
    polls[polls$state == myState & polls$year == myYear,"actualDem"] = demVotes/totalVotes*100
    polls[polls$state == myState & polls$year == myYear,"actualRep"] = repVotes/totalVotes*100
  
    demPopVote = demPopVote + demVotes
    repPopVote = repPopVote + repVotes
    totalPopVote = totalPopVote + totalVotes
    
  }
  
  
  polls[polls$state == "--" & polls$year == myYear,"actualDem"] = demPopVote/totalPopVote*100
  polls[polls$state == "--" & polls$year == myYear,"actualRep"] = repPopVote/totalPopVote*100
}

polls = transform(polls,actualMargin = actualDem - actualRep)
polls = transform(polls,demBias = dem - actualDem)
polls = transform(polls,repBias = rep - actualRep)


# make sure all populations match, and only include lv and rv
polls[polls$pop == "Likely Voters","pop"] = "lv"
polls[polls$pop == "lv ","pop"] = "lv"
polls[polls$pop == "Registered Voters","pop"] = "rv"
polls[polls$pop == "rv ","pop"] = "rv"
polls = subset(polls, pop == "lv" | pop == "rv")

polls[polls$mode == "Internet","mode"] = "Online"
polls = subset(polls, mode != "Mail" & mode != "Mixed" & mode != "Text/Online" & mode != "Live Phone/Text")


polls = subset(polls, daysuntil < 8 | year == "2020")

# first set of slides: proof of sampling error / overall bias

myYear = 2016
myState = "PA"

mySub = subset(polls, year == myYear & state == myState)
hist(mySub$dem,xlim = c(35,65), breaks = seq(34.5,65.5,1), main = "Poll Results for Clinton in Final Week (2016, Pennsylvania)",xlab = "Poll Result", col = "skyblue")
stripchart(mySub$dem,method = "jitter", add = TRUE, pch = 23, bg = "blue")
abline(v = mySub$actualDem[1], lwd = 3, col = "black")




myYear = 2012
myState = "OH"

mySub = subset(polls, year == myYear & state == myState)
hist(mySub$dem,xlim = c(35,65), breaks = seq(34.5,65.5,1), main = "Poll Results for Obama in Final Week (2012, Ohio)",xlab = "Poll Result", col = "skyblue")
stripchart(mySub$dem,method = "jitter", add = TRUE, pch = 23, bg = "blue")
abline(v = mySub$actualDem[1], lwd = 3, col = "black")





myYear = 2012
myState = "NH"

mySub = subset(polls, year == myYear & state == myState)
hist(mySub$dem,xlim = c(35,65), breaks = seq(34.5,65.5,1), main = "Poll Results for Obama in Final Week (2012, New Hampshire)",xlab = "Poll Result", col = "skyblue")
stripchart(mySub$dem,method = "jitter", add = TRUE, pch = 23, bg = "blue")
abline(v = mySub$actualDem[1], lwd = 3, col = "black")









# for (myYear in years){
#   for (myState in states){
# mySub = subset(polls, year == myYear & state == myState)
# 
# if(nrow(mySub) > 10){
# #hist(mySub$dem,xlim = c(35,65), breaks = seq(34.5,65.5,1), main = "Poll Results for Obama in Final Week (2012, Ohio)",xlab = "Poll Result", col = "skyblue")
# hist(mySub$dem, main = paste(toString(myYear),myState),xlab = "Poll Result", col = "skyblue")
# abline(v = mySub$actualDem[1], lwd = 3, col = "yellow")
# stripchart(mySub$dem,method = "jitter", add = TRUE, pch = 23, bg = "blue")
# }
# }}



# 2020 polls

myYear = "2020"
myState = "PA"

mySub = subset(polls, year == myYear & state == myState)
hist(mySub$dem,xlim = c(35,65), breaks = seq(34.5,65.5,1), main = "Poll Results for Joe Biden (2020, Pennsylvania)",xlab = "Poll Result", col = "skyblue")
stripchart(mySub$dem,method = "jitter", add = TRUE, pch = 23, bg = "blue")



# Load the relevant model into R's working memory:
source("RScripts/pollsterBiasModel.R")

myPollsters = c("SurveyUSA","Rasmssn","PblcPlc","Quinnpc")

dir.create("Markdown/Figures/Jags-ThompsonSeminar/")
dir.create("Simulations/Jags-ThompsonSeminar/")

fileNameRoot = "Markdown/Figures/Jags-ThompsonSeminar/" 
fileNameRootSim= "Simulations/Jags-ThompsonSeminar/" 
graphFileType = "png"

mydata = droplevels(subset(polls, year != "2020"))

# Generate the MCMC chain:
demPrimarySim = runSimulation( mydata , party1 = "demBias", party2 = "repBias",
                               numSavedSteps=100000 , thinSteps=1 , saveName=fileNameRootSim) 
plotDiagnostics(demPrimarySim)


#Get some plotting done
plotParty1PosteriorPredictive( demPrimarySim, p1Name = "demBias",
                               datFrm=mydata, whichPollsters = myPollsters, saveName=fileNameRoot , saveType=graphFileType)

plotMarginalDistributions( demPrimarySim, datFrm = mydata, whichPollsters = myPollsters,
                           p1Name = "demBias", saveName=fileNameRoot , saveType=graphFileType )







# weighted average

mydata = mydata[!is.na(mydata$n),]
mydata = droplevels(mydata)

newdata <- mydata[order(mydata$race_id),]
newdata2 = data.frame(newdata)

#create small data frame for reference
races = unique(newdata$race_id)
raceFullName = rep(NA,length(races))
for(i in 1:length(races)){
  mysub = subset(newdata, race_id == races[i])
  mysub[mysub$state == "--","state"] = "National"
  raceFullName[i] = paste(mysub[1,"year"],mysub[1,"state"])
}
onlyUnique=(data.table::setDT(newdata2)[,.SD[which.max(actualMargin)],keyby=race_id])
actual = onlyUnique$actualMargin
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

predictorsframe = newdata


myDataFrame=predictorsframe

#which race indexing:
whichrace=match(unique(predictorsframe$race_id), predictorsframe$race_id)
whichrace=whichrace-1
whichrace=c(whichrace,nrow(newdata))


fileNameRootSim = "Simulations/thompsonSeminarWeightedandCorrected" 
#fileNameRoot = "Markdown/Figures/WeightedModels/WeightedandCorrected/V05" 
fileNameRoot = "Markdown/Figures/" 
graphFileType = "png" 

myDataFrame$n =myDataFrame $ n/1000
#myDataFrame$samplesize = myDataFrame$log(samplesize)

#MarginOfError = sqrt(.25/myDataFrame$samplesize)*100
#myDataFrame$samplesize =MarginOfError
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
#source("Jags-Ymet-Xnom1met1-MnormalHom.R")
source("RScripts/ThompsonWeightedAverageScript.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( refdataframe ,myDataFrame, pollName="dem" , #daysuntilName="daysuntil", 
                    raceIDName="races", actualName="actual",
                    modeName = "mode", popName = "pop",
                    samplesizeName ="n", whichrace, yearName="year", pollsterName="pollster",
                    biasName="demBias",numSavedSteps=50000 , thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) 
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in c("actualSpread",   
                   "modeImpact[1]", "samplesizeImpact" , "mu[1]","LVImpact[1]") ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}

#------------------------------------------------------------------------------- 
#plot the posterior predictive distrubtions

plotPosteriorPredictive(codaSample=mcmcCoda, refFrame=refdataframe, datFrmPredictor = myDataFrame, 
                        pollName = "margin", raceIDName="races", raceplots=c(1,38,76,3,39,78), whichrace=whichrace, saveName=fileNameRoot , 
                        saveType=graphFileType)

meanPosterior(codaSample=mcmcCoda, refFrame=refdataframe, datFrmPredictor = myDataFrame, 
              pollName = "margin", raceIDName="races", raceplots=1:2, whichrace=whichrace, saveName=fileNameRoot , 
              saveType=graphFileType)
#plot the LV distribution
plotLVPosterior(codaSample=mcmcCoda,  datFrmPredictor = myDataFrame , saveName=fileNameRoot , 
                saveType=graphFileType)

plotModePosterior(codaSample=mcmcCoda,  datFrmPredictor = myDataFrame , saveName=fileNameRoot , 
                  saveType=graphFileType)

#plot the samplesize 
plotSampleSizePosterior(mcmcCoda, datFrm=myDataFrame,  saveName=fileNameRoot , 
                        saveType=graphFileType, title="Sample Size Impact ")
