mydata <- read.csv("raw-polls_538.csv")
pollsters= c("Rasmussen Reports/Pulse Opinion Research","Monmouth University",
             "Marist College", "ABC News/The Washington Post","SurveyUSA")

pollsterOfChoice1="Rasmussen Reports/Pulse Opinion Research"
pollsterOfChoice2="Monmouth University"
pollsterOfChoice3="Marist College"
pollsterOfChoice4="ABC News/The Washington Post"
pollsterOfChoice5="SurveyUSA"




#use lubridate to change dates from character to date data type for functionality. 
mydata$electiondate = lubridate::mdy(mydata$electiondate)
mydata$polldate = lubridate::mdy(mydata$polldate)

#Subset for constraints given above.
mydata=(subset(mydata, mydata$year>=earlyYear))
mydata=(subset(mydata, mydata$year<=lateYear))


mydata1=(subset(mydata,pollster==pollsterOfChoice1))
mydata2=(subset(mydata,pollster==pollsterOfChoice2))
mydata3=(subset(mydata,pollster==pollsterOfChoice3))
mydata4=(subset(mydata,pollster==pollsterOfChoice4))
mydata5=(subset(mydata,pollster==pollsterOfChoice5))


#Use only the most recent poll for each election:
#Use setDT function from data.table package to get a subset from mydata with just the max. value of the date element for each race (grouped with the keyby function). Call this new subset of data onlyRecentData.
onlyRecentData1=(data.table::setDT(mydata1)[,.SD[which.max(polldate)],keyby=race_id])
onlyRecentData2=(data.table::setDT(mydata2)[,.SD[which.max(polldate)],keyby=race_id])
onlyRecentData3=(data.table::setDT(mydata3)[,.SD[which.max(polldate)],keyby=race_id])
onlyRecentData4=(data.table::setDT(mydata4)[,.SD[which.max(polldate)],keyby=race_id])
onlyRecentData5=(data.table::setDT(mydata5)[,.SD[which.max(polldate)],keyby=race_id])
#combine all the rows
total = rbind(onlyRecentData1,onlyRecentData2, onlyRecentData3, onlyRecentData4, onlyRecentData5)

df1 <- data.frame(matrix(ncol = 7, nrow = 0))
x <- c("Race","Election Date",pollsterOfChoice1, pollsterOfChoice2, pollsterOfChoice3,pollsterOfChoice4, pollsterOfChoice5)
colnames(df1) <- x


for(myRace in unique(total$race_id))
{

subrace=subset(total, race_id==myRace)
head(subrace)


myYear=subrace$year[1]


p1bias=subset(subrace, pollster==pollsterOfChoice1)$bias
if(length(p1bias)==0){
  p1bias=NA
}
p2bias=subset(subrace, pollster==pollsterOfChoice2)$bias
if(length(p2bias)==0){
  p2bias=NA
}
p3bias=subset(subrace, pollster==pollsterOfChoice3)$bias
if(length(p3bias)==0){
  p3bias=NA
}
p4bias=subset(subrace, pollster==pollsterOfChoice4)$bias
if(length(p4bias)==0){
  p4bias=NA
}
p5bias=subset(subrace, pollster==pollsterOfChoice5)$bias
if(length(p5bias)==0){
  p5bias=NA
}

newline=data.frame("Race"=myRace, "Election Date"=myYear,
                   pollsterOfChoice1=p1bias, 
                   pollsterOfChoice2=p2bias, 
                   pollsterOfChoice3=p3bias, 
                   pollsterOfChoice4=p4bias,
                   pollsterOfChoice5=p5bias)

df1=rbind(df1, newline)
}


lattice::splom(subset(df1, select=c(pollsterOfChoice1, pollsterOfChoice2,
                                    pollsterOfChoice3, pollsterOfChoice4, pollsterOfChoice5)))



