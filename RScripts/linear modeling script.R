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

mydata$year = as.factor(mydata$year)

#create empty data frame, while maintaining all columns from the mydata structure
myDataMyPollsters=  mydata[0,]

#run a loop to fill this data frame with every pollster
for(myPollster in unique(mydata$pollster))
{
  #subset to a dataset with just each pollster
  subpoll=subset(mydata, pollster==myPollster)
  
  #if that subset has more than thirty entries, use it:
  if(nrow(subpoll)>300){
  #Use only the most recent poll for each election:
  #Use setDT function from data.table package to get a subset from mydata with just the max. value of the date element for each race (grouped with the keyby function). Call this new subset of data onlyRecentData.
  onlyRecentData1=(data.table::setDT(subpoll)[,.SD[which.max(polldate)],keyby=race_id])
  #combine each of these subsets together
  myDataMyPollsters=rbind(myDataMyPollsters, onlyRecentData1)
  }
}

model = lm(bias~ year + pollster + type_simple + partisan + pollster:year , data = myDataMyPollsters)
summary(model)
# R^2 value = 19.9% - not bad!

#SVAs
car::qqp(rstandard(model))
lattice::densityplot(rstandard(model))
# the residual density plot looks cleary centered on zero and mostly normal with outliers
# robust linear modeling will probably do a very good job here when we construct our Bayesian model

model = lm(bias~ year + pollster + type_simple + partisan , data = myDataMyPollsters)
summary(model)



#plot years and pollster
# Convert cyl column from a numeric to a factor variable
myDataMyPollsters$pollster <- as.factor(myDataMyPollsters$pollster)

#plot by year and pollster
library(ggplot2)
ggplot(myDataMyPollsters, aes(x=year, y=bias, color=pollster)) +
  geom_point() +
  geom_smooth(method=lm)
# Remove confidence intervals
# Extend the regression lines


#helpful to interact year and pollster
model1 = lm(bias~ year + pollster + type_simple +  year:pollster + partisan , data = myDataMyPollsters)
summary(model1)



#plot by year and type_simple
library(ggplot2)
ggplot(myDataMyPollsters, aes(x=year, y=bias, color=type_simple)) +
  geom_point() +
  geom_smooth(method=lm, , se=FALSE, fullrange=TRUE)
# Remove confidence intervals
# Extend the regression lines

#interaction term with year and type needed?
model2 = lm(bias~  year:type_simple + year:pollster, data = myDataMyPollsters)
summary(model2)


model5 = lm(bias~ year:simple_type, data = myDataMyPollsters)
summary(model5)


#plot by year and partisian 

library(ggplot2)
ggplot(myDataMyPollsters, aes(x=year, y=bias, color=partisan)) +
  geom_point() +
  geom_smooth(method=lm, , se=FALSE, fullrange=TRUE)
# Remove confidence intervals
# Extend the regression lines

#interaction with year and partisan needed?
model3 = lm(bias~ year + pollster + year:partisan + type_simple + year:type_simple +  year:pollster + partisan , data = myDataMyPollsters)
summary(model3)
#no thats not better

#sample size?
model4 = lm(bias~ year + pollster + year:partisan + type_simple + year:type_simple +  year:pollster + partisan + samplesize , data = myDataMyPollsters)
summary(model4)

#plot by sample size and pollster


library(ggplot2)
ggplot(myDataMyPollsters, aes(x=samplesize, y=bias, color=pollster)) +
  geom_point() +
  geom_smooth(method=lm)


data5 = transform(myDataMyPollsters, logconc = log(samplesize))
#sample size doesn't make it better

ggplot(data5, aes(x=logconc, y=bias, color=pollster)) +
  geom_point() +
  geom_smooth(method=lm)

model5 = lm(bias ~ year*pollster + year*type_simple + year + type_simple + pollster + logconc , data = data5)
summary(model5)


modela <- lm(bias ~ year*pollster + year*type_simple + year + type_simple + pollster + samplesize, data = myDataMyPollsters)
summary(modela)

