#Clean Up Data
#Author: Haley Reed
#Date: June 22nd 2020

#Read the Data - FiveThirtyEight Data
mydata=read.csv("Data/raw-polls_538.csv")

#Limit to just president
preferredType="Pres-G"
mydata=subset(mydata, type_simple==preferredType)

#Treat pollster as a factor.
mydata$pollster = as.factor(mydata$pollster)
#Rename pollsters with special characters
levels(mydata$pollster) = stringr::str_replace_all(levels(mydata$pollster),"/","-")


#Make copy of full pollster name for reference
mydata$pollsterFull=mydata$pollster
#Abbreviate to limit Pollster names to a reasonable length.
levels(mydata$pollster)=abbreviate(levels(mydata$pollster), minlength = 7, use.classes = TRUE,
           dot = FALSE, strict = FALSE)


#use lubridate to change dates from character to date data type for functionality. 
mydata$electiondate = lubridate::mdy(mydata$electiondate)
mydata$polldate = lubridate::mdy(mydata$polldate)

library(magrittr)
#create daysuntil column to keep track of the number of days between a poll and the actual election.
daysuntil=lubridate::interval(mydata$polldate,mydata$electiondate)
daysuntil=lubridate::as.period(daysuntil) %>% lubridate::day()
mydata$daysuntil=daysuntil                   

#subset to make sure it is a democratic v. republican race
mydata=subset(mydata,cand1_party=="DEM")
mydata=subset(mydata,cand2_party=="REP")

#Treat year as a factor
mydata$year = as.factor(mydata$year)


#Change so that bias is candiate specific, rather than in terms of margins.

demBias = mydata[,15] - mydata[,22]
mydata$demBias=demBias

repBias = mydata[,18] - mydata [,23]
mydata$repBias=repBias









#Add pollster credibility stuff!!

pollsterCred=read.csv("Data/pollster-ratings.csv")

  #start with transparency!
myDataMyPollsters=  mydata[0,]
for(myID in unique(mydata$pollster_rating_id))
{
  #find the transparency of that pollster
  subcred = subset(pollsterCred, Pollster.Rating.ID==myID)
  transparency = subcred$NCPP...AAPOR...Roper
  
  #subset to a dataset with just each pollster
  subpoll=subset(mydata, pollster_rating_id==myID)
  transparencyList = subpoll [ ,0 ]
  transparencyList=rep(transparency,nrow(subpoll))
  subpoll$transparency = transparencyList

    myDataMyPollsters=rbind(myDataMyPollsters, subpoll)
  
}
myDataFrame = myDataMyPollsters
myDataFrame$pollster = factor( myDataFrame$pollster)
mydata = myDataFrame




    #next to delMode 
myDataMyPollsters=  mydata[0,]
for(myID in unique(mydata$pollster_rating_id))
{
  #find the mode of that pollster
  submethod = subset(pollsterCred, Pollster.Rating.ID==myID)
  method = submethod$Methodology
  
  #subset to a dataset with just each pollster
  subpoll=subset(mydata, pollster_rating_id==myID)
  methodList = subpoll [ ,0 ]
  methodList=rep(method,nrow(subpoll))
  subpoll$delMode = methodList
  
  myDataMyPollsters=rbind(myDataMyPollsters, subpoll)
  
}
myDataFrame = myDataMyPollsters
myDataFrame$pollster = factor( myDataFrame$pollster)
mydata = myDataFrame




  #next is LV
myDataMyPollsters=  mydata[0,]
for(myPoll in unique(mydata$poll_id))
{
  #subset to a dataset with just each pollster
  subpoll=subset(mydata, poll_id==myPoll)
  comment1=subpoll[1,28]
  
  #if that subset has more than thirty entries, use it:
  if(grepl("register", comment1, ignore.case = TRUE)){
  LV=FALSE
  LVList = subpoll [ ,0 ]
  LVList=rep(LV,nrow(subpoll))
  subpoll$LV = LVList
  
  myDataMyPollsters=rbind(myDataMyPollsters, subpoll)
  
  }
  else{
    LV=TRUE
    LVList = subpoll [ ,0 ]
    LVList=rep(LV,nrow(subpoll))
    subpoll$LV = LVList
    
    myDataMyPollsters=rbind(myDataMyPollsters, subpoll)
  }
}
myDataFrame = myDataMyPollsters
myDataFrame$pollster = factor( myDataFrame$pollster)

#add an undecided voter column
Undecided = 100 - myDataFrame$cand1_pct - myDataFrame$cand2_pct
myDataFrame$undecided = Undecided
mydata = myDataFrame


#Save the csv file.
write.csv(mydata,'Data/raw-polls_538_cleaned.csv')


#Make a version of the data that includes polls exclusively in the week prior to the election
recentdata = mydata[0,]
for(i in unique(mydata$daysuntil))
  {
    #subset to a dataset with just each day until
    subpoll=subset(mydata, daysuntil==i)

    #if that subset has more than thirty entries, use it:
    if(i<8){
      #combine each of these subsets together
      recentdata=rbind(recentdata, subpoll)
    }
  }
write.csv(recentdata,'Data/raw-polls_538_weekprior.csv')

