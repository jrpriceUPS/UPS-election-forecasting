#Clean Up Data
#Author: Haley Reed
#Date: June 22nd 2020

#Read the Data - FiveThirtyEight Data
mydata=read.csv("Data/raw-polls_538.csv")

#Limit to just president
preferredType="Pres-G"

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

#Save the csv file.
write.csv(mydata,'Data/raw-polls_538_cleaned.csv')


