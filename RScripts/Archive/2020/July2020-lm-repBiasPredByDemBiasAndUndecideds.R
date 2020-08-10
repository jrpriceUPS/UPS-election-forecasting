# Small Fisherian exploration of predicting rep bias using dem bias and undecideds
# provided inspiration for JAGS Party Error Predictor scripts
#
# Runs!
# Notes by Jake 8/10/20

#Fisherian Analyis of undecideds

#create empty data frame, while maintaining all columns from the mydata structure
#Pick min. number of acceptable polls
mydata=read.csv("Data/raw-polls_538_cleaned.csv")

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

Undecided = 100 - myDataFrame$cand1_pct - myDataFrame$cand2_pct
myDataFrame$undecided = Undecided


Model1 = lm(repBias ~ demBias*year + undecided, data=myDataFrame)
summary(Model1)


Model2 = lm(repBias ~ demBias*year + demBias*pollster + undecided, data=myDataFrame)
summary(Model2)

Model3 = lm(repBias ~ demBias*year*pollster + undecided, data=myDataFrame)
summary(Model3)
