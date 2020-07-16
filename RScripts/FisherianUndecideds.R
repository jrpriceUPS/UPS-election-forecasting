#Fisherian Analyis of undecideds

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

myDataFrame = myDataMyPollsters

myDataFrame$pollster = factor( myDataFrame$pollster)

Undecided = 100 - myDataFrame$cand1_pct - myDataFrame$cand2_pct
myDataFrame$undecided = Undecided
