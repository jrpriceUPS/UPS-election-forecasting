# Updated data cleaning script to convert raw-polls_538.csv to appropriate form
# Date: August 11, 2020
# Author: Jake Price and Haley Reed

# original version:
#Clean Up Data
#Author: Haley Reed
#Date: June 22nd 2020

# import magrittr for piping ( %>% )
library(magrittr)


# Read the data - FiveThirtyEight Data
pollingData = read.csv("Data/raw-polls_538.csv") # should update this to pull from github repository at some point

# Only general election presidential polls
genElectionPolls = subset(pollingData, type_simple=="Pres-G")


# Rename pollsters with special characters
genElectionPolls$pollster = as.factor(genElectionPolls$pollster)
levels(genElectionPolls$pollster) = stringr::str_replace_all(levels(genElectionPolls$pollster),"/","-")


# Make copy of full pollster name for reference
genElectionPolls$pollsterFullName = genElectionPolls$pollster
#Abbreviate to limit Pollster names to a reasonable length.
levels(genElectionPolls$pollster) = abbreviate(levels(genElectionPolls$pollster), minlength = 7, use.classes = TRUE,
           dot = FALSE, strict = FALSE)


#use lubridate to change dates from character to date data type for functionality. 
genElectionPolls$electiondate = lubridate::mdy(genElectionPolls$electiondate)
genElectionPolls$polldate = lubridate::mdy(genElectionPolls$polldate)

#create daysuntil column to keep track of the number of days between a poll and the actual election.
genElectionPolls$daysuntil = as.numeric(difftime(genElectionPolls$electiondate,genElectionPolls$polldate, units = "days"))        

#subset to make sure it is a democratic v. republican race
genElectionPolls=subset(genElectionPolls,cand1_party=="DEM")
genElectionPolls=subset(genElectionPolls,cand2_party=="REP")

#Treat year as a factor
genElectionPolls$year = as.factor(genElectionPolls$year)


#Change so that bias is candidate specific, rather than in terms of margins.
genElectionPolls$demBias = genElectionPolls[,"cand1_pct"] - genElectionPolls[,"cand1_actual"]
genElectionPolls$repBias = genElectionPolls[,"cand2_pct"] - genElectionPolls[,"cand2_actual"]




#Add pollster credibility information!

# HUGE NOTE - this assumes that a given pollster *always* uses the same delivery mode and LV / RV model
# that is almost certainly not true - it would be best to set up a better database for polling that is not reliant on 538
# I've reached out to 538 to see if they can get us the data we need
# another option: Elliott Morris and the Economist's data: https://github.com/TheEconomist/us-potus-model/tree/master/data

pollsterCred=read.csv("Data/pollster-ratings.csv")

# append a column coding transparency, delivery mode, and population
genElectionPolls = transform(genElectionPolls, transparency = 0, IVR = FALSE, online = FALSE, text = FALSE, live = FALSE, LV = 0)

for(myID in unique(genElectionPolls$pollster_rating_id))
{
  #find the transparency and delivery mode of that pollster
  isTransparent = subset(pollsterCred, Pollster.Rating.ID==myID)$NCPP...AAPOR...Roper
  delMode = subset(pollsterCred, Pollster.Rating.ID==myID)$Methodology
  
  # relabel transparency for all rows corresponding to that pollster
  genElectionPolls[genElectionPolls$pollster_rating_id == myID,"transparency"] = isTransparent
  
  # code delivery mode for that pollster
  genElectionPolls[genElectionPolls$pollster_rating_id == myID,"online"] = (delMode=="IVR/Online"||delMode=="IVR/Online/Live"||delMode=="IVR/Online/Live/Text"||delMode=="IVR/Online/Text"||delMode=="Online"||delMode=="Online/Live"||delMode =="Online/Live/Text"||delMode=="Online/Text")
  genElectionPolls[genElectionPolls$pollster_rating_id == myID,"live"] = (delMode == "IVR/Live"||delMode=="IVR/Live/Text"||delMode=="IVR/Online/Live"||delMode=="IVR/Online/Live/Text"||delMode=="Live"||delMode=="Live*"||delMode=="Live/Text"||delMode=="Online/Live"||delMode=="Online/Live/Text")
  genElectionPolls[genElectionPolls$pollster_rating_id == myID,"IVR"] = (delMode=="IVR"||delMode=="IVR/Live"||delMode=="IVR/Live/Text"||delMode=="IVR/Online"||delMode=="IVR/Online/Live"||delMode=="IVR/Online/Live/Text"||delMode=="IVR/Online/Text"||delMode=="IVR/Text")
  genElectionPolls[genElectionPolls$pollster_rating_id == myID,"text"] = (delMode=="IVR/Live/Text"||delMode=="IVR/Online/Live/Text"||delMode=="IVR/Online/Text"||delMode=="IVR/Text"||delMode=="Live/Text"||delMode=="Online/Live/Text"||delMode=="Online/Text")
  
}


# next is LV
# This should also absolutely be updated!

genElectionPolls$LV = !grepl("register", genElectionPolls$comment, ignore.case = TRUE)


#add an undecided voter column
genElectionPolls$undecided = 100 - genElectionPolls$cand1_pct - genElectionPolls$cand2_pct


#Save the csv file.
write.csv(genElectionPolls,'Data/raw-polls_538_cleaned.csv')


# Make a version of the data that includes polls exclusively in the week prior to the election
finalPolls = subset(genElectionPolls,daysuntil < 8)
write.csv(finalPolls,'Data/raw-polls_538_weekprior.csv')