#Linear Modeling of Polling Data 


#load the data
#in this case from 538
mydata <- read.csv("raw-polls_538.csv")

#Establish which years you wish to use data from (between 1998 and 2020):
earlyYear=1998
lateYear=2020
mydata=(subset(mydata, mydata$year>=earlyYear))
mydata=(subset(mydata, mydata$year<=lateYear))

#use lubridate to change dates from character to date data type for functionality. 
mydata$electiondate = lubridate::mdy(mydata$electiondate)
mydata$polldate = lubridate::mdy(mydata$polldate)

#subset to make sure it is a democratic v. republican race
mydata=subset(mydata,cand1_party=="DEM")
mydata=subset(mydata,cand2_party=="REP")



#create empty data frame, while maintaining all columns from the mydata structure
myDataMyPollsters=  mydata[0,]

#run a loop to fill this data frame with every pollster
for(myPollster in unique(mydata$pollster))
{
  #subset to a dataset with just each pollster
  subpoll=subset(mydata, pollster==myPollster)
  
  #if that subset has more than thirty entries, use it:
  if(nrow(subpoll)>200){
  #Use only the most recent poll for each election:
  #Use setDT function from data.table package to get a subset from mydata with just the max. value of the date element for each race (grouped with the keyby function). Call this new subset of data onlyRecentData.
  onlyRecentData1=(data.table::setDT(subpoll)[,.SD[which.max(polldate)],keyby=race_id])
  #combine each of these subsets together
  myDataMyPollsters=rbind(myDataMyPollsters, onlyRecentData1)
  }
}

model = lm(bias~ year + pollster + type_simple + partisan , data = myDataMyPollsters)
