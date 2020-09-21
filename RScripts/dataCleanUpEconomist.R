#Data Clean Up file to Use the Economist data

#Orginal Version 
#Date: September 13th 2020
#Author: Haley Reed

# Read the data - Economist Data
pollingData2008 = read.csv("Data/all_polls_2008.csv") 
pollingData2012 = read.csv("Data/all_polls_2012.csv") 
pollingData2016 = read.csv("Data/all_polls_2016.csv") 

# Rename pollsters with special characters
pollingData2008$pollster = as.factor(pollingData2008$pollster)
pollingData2012$pollster = as.factor(pollingData2012$pollster)
pollingData2016$pollster = as.factor(pollingData2016$pollster)


# Make copy of full pollster name for reference
pollingData2008$pollsterFullName = pollingData2008$pollster
pollingData2012$pollsterFullName = pollingData2012$pollster
pollingData2016$pollsterFullName = pollingData2016$pollster

#Condense to one dataframe.
pollingData2008$year = "2008"
pollingData2012$year = "2012"
pollingData2016$year = "2016"

#create a 2016 dataframe that has the same variables as 2008 and 2012.
pollingData2016Con <- data.frame(pollingData2016$state, pollingData2016$pollster,
                                 pollingData2016$number.of.observations, pollingData2016$start.date,
                                 pollingData2016$end.date, pollingData2016$trump, pollingData2016$clinton,
                                 pollingData2016$other,pollingData2016$undecided, pollingData2016$mode,
                                 pollingData2016$population,pollingData2016$pollsterFullName,
                                 pollingData2016$year)

#correct the names of the variables
pollingData2016Con$state=pollingData2016Con$pollingData2016.state
pollingData2016Con$pollster=pollingData2016Con$pollingData2016.pollster
pollingData2016Con$number.of.observations=pollingData2016Con$pollingData2016.number.of.observations
pollingData2016Con$start.date=pollingData2016Con$pollingData2016.start.date
pollingData2016Con$end.date=pollingData2016Con$pollingData2016.end.date
pollingData2016Con$other=pollingData2016Con$pollingData2016.other
pollingData2016Con$undecided=pollingData2016Con$pollingData2016.undecided
pollingData2016Con$mode=pollingData2016Con$pollingData2016.mode
pollingData2016Con$population=pollingData2016Con$pollingData2016.population
pollingData2016Con$pollsterFullName=pollingData2016Con$pollingData2016.pollsterFullName
pollingData2016Con$year=pollingData2016Con$pollingData2016.year
pollingData2016Con$democrat=pollingData2016Con$pollingData2016.clinton
pollingData2016Con$republican=pollingData2016Con$pollingData2016.trump

pollingData2016Con = pollingData2016Con[,-(1:13)]  

#standardize the vote share variables - for 2008
pollingData2008$democrat=pollingData2008$obama           
pollingData2008$republican=pollingData2008$mccain           
pollingData2008 <- pollingData2008[,-(6:7)]  

#for 2012...
pollingData2012$democrat=pollingData2012$obama           
pollingData2012$republican=pollingData2012$romney           
pollingData2012 <- pollingData2012[,-(6:7)]  


newFrame = rbind(pollingData2008,pollingData2012,pollingData2016Con)
newFrame$race_id <- seq.int(nrow(newFrame))

write.csv(newFrame,'Data/economist_cleaned.csv')
