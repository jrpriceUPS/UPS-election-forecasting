#Comparing Pollsters
#While looking at one pollster may provide us with some relevant information, it is much more interesting to look at multiple pollsters and compare them. Choose your group below:


#load all the data
mydata <- read.csv("raw-polls_538.csv")

#Establish which years you wish to use data from (between 1998 and 2020):
earlyYear=1998
lateYear=2020

#use lubridate to change dates from character to date data type for functionality. 
mydata$electiondate = lubridate::mdy(mydata$electiondate)
mydata$polldate = lubridate::mdy(mydata$polldate)

#Subset for year constraints given above.
mydata=(subset(mydata, mydata$year>=earlyYear))
mydata=(subset(mydata, mydata$year<=lateYear))


#list all pollsters you wish to compare.
pollsters= c("Rasmussen Reports/Pulse Opinion Research","Monmouth University",
             "Marist College", "ABC News/The Washington Post","SurveyUSA", "Quinnipiac University")

#create empty data frame, while maintaining all columns from the mydata structure
myDataMyPollsters=  mydata[0,]

#run a loop to fill this data frame with every pollster specified above
for(myPollster in unique(pollsters))
{
  #subset to a dataset with just each pollster
  subpoll=subset(mydata, pollster==myPollster)
  #Use only the most recent poll for each election:
  #Use setDT function from data.table package to get a subset from mydata with just the max. value of the date element for each race (grouped with the keyby function). Call this new subset of data onlyRecentData.
  onlyRecentData1=(data.table::setDT(subpoll)[,.SD[which.max(polldate)],keyby=race_id])
  #combine each of these subsets together
  myDataMyPollsters=rbind(myDataMyPollsters, onlyRecentData1)
}



#create new data frame.
df1 <- data.frame(matrix(ncol = 7, nrow = 0))
#First two columns will just be Race and Election Date
x <- c("Race","Election Date")


#define x as the column names for our new dataframe.
colnames(df1) <- x

#loop through pollster filtered data to fill the data frame
for(myRace in unique(myDataMyPollsters$race_id))
{
  #subset the data to only one race at a time
  subrace=subset(myDataMyPollsters, race_id==myRace)
  #Choose a year by picking year of first data point in the frame
  myYear=subrace$year[1]
  
  
  #create new data frame to store info from each race
  newline=data.frame("Race"=myRace, "Election Date"=myYear)
  
  #Loop through each pollster to find their bias
  for(myPollster in unique(pollsters))
  {
    #get bias for that pollster
    pbias=subset(subrace, pollster==myPollster)$bias
    if(length(pbias)==0){
      pbias=NA
    }
    #add that pollster's bias to individual race's data frame
    newline = cbind(newline, pbias)
  }
  
  
  #add the individual race data frames together 
  df1=rbind(df1, newline)
}

#loop through list of desired pollsters for the rest of the columns
for(myPollster in unique(pollsters))
{
  #Append pollster names to column name list
  x=c(x,myPollster)
}

#define x as the column names for our new dataframe.
colnames(df1) <- x



#scatterplot to visualize correlation
lattice::splom(subset(df1, select=pollsters))


#alternative scatterplot
# ggplot2::ggplot(df1) +
#    ggplot2::geom_point(aes(x = .panel_x, y = .panel_y)) +
#    ggforce::facet_matrix(vars(pollsters))

#density plot

lattice::densityplot(~bias,data=myDataMyPollsters,
            groups=pollster,
            xlab="Bias",
            main="Pollster Bias",
            plot.points=TRUE,
            auto.key=TRUE)



#model = lm(bias~ year + pollster + type_simple + partisan + year:pollster + pollster:type_simple, data = mydata)
#model = lm(bias~ year + pollster + type_simple + partisan , data = mydata)
