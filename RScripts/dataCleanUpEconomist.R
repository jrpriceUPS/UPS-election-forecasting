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

pollingData2016Con <- pollingData2016[state,pollster,number.of.observations,start.date,
                                      end.date,trump,clinton,other,undecided,mode,
                                      population,pollsterFullName,year]

newFrame = rbind(pollingData2008,pollingData2012,pollingData2016)
