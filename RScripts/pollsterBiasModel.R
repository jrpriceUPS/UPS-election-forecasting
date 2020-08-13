# Author: Jake Price
# Date: 8/12/20
#
# A script that models the bias of pollster in a given year ("house effect")
#  * Models demBias and repBias separately 
#     - one is "primary" and other is "secondary"
#     - secondary is modeled on primary and undecided voters in poll
#  * Overall year lean pulled from parent distribution of year leans
#  * Pollster bias in a given year comes from parent distribution for that pollster
#     - Pollster parent dist. comes from larger parent distribution of pollster centers

# Let us predict rebBias with demBias. 


# source summary functions (look at this later and pull in only the functions we need)
source("DBDA2E-utilities.R")

#===============================================================================
runSimulation = function( datFrm , party1 = "demBias", party2 = "repBias",
                          numSavedSteps=50000 , thinSteps=1 , saveName=NULL,
                          runjagsMethod=runjagsMethodDefault , nChains=nChainsDefault ) 
  { 
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to lists for delivery to JAGS
  undecided = as.numeric(datFrm[,"undecided"])
  p1bias = as.numeric(datFrm[,party1])
  p2bias = as.numeric(datFrm[,party2])
  pollster = as.numeric(as.factor(datFrm[,"pollster"]))
  PollsterLevels = levels(as.factor(datFrm[,"pollster"]))
  year = as.numeric(as.factor(datFrm[,"year"]))
  YearLevels = levels(as.factor(datFrm[,"year"]))
  
  # indexing
  PollsTotal = length(p1bias)
  PollsterLevelsTotal = length(unique(pollster))
  YearLevelsTotal = length(unique(year))
  
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # For prior on baseline, etc.:
  p1mean = mean(p1bias)
  p1SD = sd(p1bias)
  p2Mean = mean(p2bias)
  p2SD = sd(p2bias)
  undecidedMean=mean(undecided)
  undecidedSD=sd(undecided)
  
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    p1bias = p1bias,
    p2bias = p2bias,
    pollster = pollster,
    undecided=undecided,
    year = year,
    
    # indexing
    PollsTotal = PollsTotal ,
    PollsterLevelsTotal = PollsterLevelsTotal ,
    YearLevelsTotal =  YearLevelsTotal ,
    
    # data properties for scaling the prior:
    p1mean = p1mean,
    p1SD = p1SD,
    p2mean = p2Mean,
    p2SD = p2SD,
    undecidedMean=undecidedMean,
    undecidedSD=undecidedSD
  )
  
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  model {
  
    #Bottom Level (individual poll bias):
    
    for (poll in 1:PollsTotal){
    
      # primary party is modeled as t-distribution with center dictated by year lean + house effects
      mu[poll] <- yearLean[year[poll]] + pollsterBias[pollster[poll],year[poll]] 
      p1bias[poll] ~ dt(mu[poll] , (1/p1PollSpread^2) , nu1 )
    
      # secondary party is modeled on primary party and undecideds
      p2bias[poll] ~ dnorm((responseFirstPartyBias*p1bias[poll] + responseUndecideds*undecided[poll]), 1/p2PollSpread^2)
    
    }
    
    # vague prior for normality of party 1 poll bias
    nu1 ~  dexp(1/30.0) 
    
    # vague priors on appropriate scale for spread of poll biases
    p1PollSpread ~ dunif( p1SD/100 , p1SD*10 )
    p2PollSpread ~ dunif( p2SD/100 , p2SD*10 )
    
    
    responseFirstPartyBias  ~ dnorm( 0 , 1/(10)^2 ) # extremely vague prior for response between parties
    responseUndecideds  ~ dnorm( 0 , 1/(10)^2 ) # vague prior for response to undecideds

   
    # Year lean branch of hierarchy
    for ( year in 1:YearLevelsTotal ) { yearLean[year] ~ dnorm( 0.0 , 1/yearSpread^2) }
    yearSpread ~ dunif( 1/1000 , 100 )
    
    # House effects for each pollster
     for ( pollster in 1:PollsterLevelsTotal ) {
     
        pollsterCenter[pollster] ~ dnorm(0.0, 1/10^2)
        pollsterSpread[pollster] ~ dunif( 1/1000, 100)
        
        for ( year in 1:YearLevelsTotal )
        { 
        
          pollsterBias[pollster,year] ~ dnorm(pollsterCenter[pollster] , 1/pollsterSpread[pollster]^2) 
          
        }
     }
  }
  
  
  " # close quote for modelstring
  writeLines(modelstring,con="TEMPmodel.txt")
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  
  # Let JAGS do parameters automatically...
  #
  #------------------------------------------------------------------------------
  # RUN THE CHAINS
  require(rjags)
  parameters = c("pollsterBias","pollsterCenter","pollsterSpread","yearSpread","yearLean","responseUndecideds","responseFirstPartyBias",
                 "p1PollSpread","p2PollSpread","nu1")
  adaptSteps = 500 
  burnInSteps = 1000 
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
}
#===============================================================================

smryMCMC = function(  codaSamples , datFrm=NULL , dembiasName = "demBias" , pollsterName = "pollster" , yearName = "year",
                      undecidedName="undecided",
                      contrasts=NULL , saveName=NULL ) {
  # All single parameters:
  parameterNames = varnames(codaSamples) 
  if ( !is.null(datFrm) & !is.null(pollsterName) & !is.null(yearName) & !is.null(undecidedName) ) {
    PollsterLevels = levels(as.factor(datFrm[,pollsterName]))
    YearLevels = levels(as.factor(datFrm[,yearName]))
  }
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName] ) )
    thisRowName = parName
    if ( !is.null(datFrm) & !is.null(pollsterName) & !is.null(yearName) & !is.null(undecidedName)) {
      # For row name, extract numeric digits from parameter name. E.g., if
      # parameter name is "b1b2[12,34]" then pull out b1b2, 12 and 34:
      strparts = unlist( strsplit( parName , "\\[|,|\\]"  ) )
      # if there are only the param name and a single index:
      if ( length(strparts)==2 ) {
        # if param name refers to factor 1:
        if ( substr(strparts[1],nchar(strparts[1]),nchar(strparts[1]))=="1" ) {
          thisRowName = paste( thisRowName , PollsterLevels[as.numeric(strparts[2])] )
        }
        # if param name refers to factor 2:
        if ( substr(strparts[1],nchar(strparts[1]),nchar(strparts[1]))=="2" ) {
          thisRowName = paste( thisRowName , YearLevels[as.numeric(strparts[2])] )
        }
      }
      # if there are the param name and two indices:
      if ( length(strparts)==3 ) {
        thisRowName = paste( thisRowName , PollsterLevels[as.numeric(strparts[2])],
                             YearLevels[as.numeric(strparts[3])] )
      }
    }
    rownames(summaryInfo)[NROW(summaryInfo)] = thisRowName
  }
  
  summaryInfo = rbind( summaryInfo , 
                       "repResponsetoDemBias" = summarizePost( mcmcMat[,"repResponsetoDemBias"]  
                                                               
                                                               #,compVal=compValrepResponsetoDemBias , 
                                                               #ROPE=roperepResponsetoDemBias 
                       ) )
  
  summaryInfo = rbind( summaryInfo , 
                       "UndecidedResponse" = summarizePost( mcmcMat[,"UndecidedResponse"] 
                                                            # ,compVal=compValrepResponsetoDemBias , 
                                                            #ROPE=roperepResponsetoDemBias 
                       ) )
  
  # Save results:
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================
c("pollsterBias","pollsterCenter","pollsterSpread","yearSpread","yearLean","responseUndecideds","responseFirstPartyBias",
  "p1PollSpread","p2PollSpread","nu1","mu")
plotDiagnostics= function( ){
  for ( parName in c("p1PollSpread",  "nu1" , "pollsterSpread" , "yearSpread" , "yearLean[1]", "pollsterBias[1,1]"  ) ) {
    diagMCMC( codaObject=mcmcCodaV03 , parName=parName , 
              saveName=fileNameRoot , saveType=graphFileType )
  }
}
#===============================================================================

#===============================================================================

plotParty1PosteriorPredictive = function( codaSamples , 
                                       p1Name = "demBias",
                                       datFrm , saveName=NULL , 
                                       saveType="jpg",
                                       showCurve = FALSE) {
  if(p1Name == "demBias"){
    myCol = "blue"
    myLabel = "Democratic Bias"}
  if(p1Name == "repBias"){
    myCol = "red"
    myLabel = "Republican Bias"}
  
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  p1bias = datFrm[,p1Name]
  pollster = as.numeric(as.factor(datFrm[,"pollster"]))
  PollsterLevels = levels(as.factor(datFrm[,"pollster"]))
  year = as.numeric(as.factor(datFrm[,"year"]))
  YearLevels = levels(as.factor(datFrm[,"year"]))
  
  # Display data with posterior predictive distributions
  for ( Yearidx in 1:length(YearLevels) ) {
    openGraph(width=2*length(PollsterLevels),height=5)
    par( mar=c(4,4,2,1) , mgp=c(3,1,0) )
    
    limits = c(min(p1bias)-0.2*(max(p1bias)-min(p1bias)),max(p1bias)+0.2*(max(p1bias)-min(p1bias))) 
    
    plot(-10,-10,
         xlim=c(0.2,length(PollsterLevels)+0.1) ,
         xlab=paste(pollsterName,yearName,sep="\n") ,
         xaxt="n" , ylab = myLabel ,
         
         ylim=limits ,
         main="Data with Post. Pred.")
    axis( 1 , at=1:length(PollsterLevels) , tick=FALSE ,
          lab=paste( PollsterLevels , YearLevels[Yearidx] , sep="\n" ) )
    for ( Pollsteridx in 1:length(PollsterLevels) ) {
      xPlotVal = Pollsteridx 
      y1Vals = p1bias[ pollster==Pollsteridx & year==Yearidx ]
      points( rep(xPlotVal,length(y1Vals))+runif(length(y1Vals),-0.05,0.05) ,
              y1Vals , pch=1 , cex=1.5 , col=myCol )
      chainSub = round(seq(1,chainLength,length=20))
      for ( chnIdx in chainSub ) {
        m = mcmcMat[chnIdx,paste("pollsterBias[",Pollsteridx,",", Yearidx,"]",sep="")]   # pollster bias
        +         mcmcMat[chnIdx,paste("yearLean[",Yearidx,"]",sep="")]  # year lean
        s = mcmcMat[chnIdx,"p1PollSpread"] # spread
        nu = mcmcMat[chnIdx,"nu1"]# normality


        tlim = qt( c(0.025,0.975) , df=nu )
        yl = m+tlim[1]*s
        yh = m+tlim[2]*s
        ycomb=seq(yl,yh,length=201)
        
        yt = dt( (ycomb-m)/s , df=nu )
        yt = 0.67*yt/max(yt)
        lines( xPlotVal-yt , ycomb , col=myCol )



      }
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PostPred-Dem-",YearLevels[Yearidx]), type=saveType)
    }
  }# end for Yearidx
}




#===============================================================================  

#===============================================================================


plotMarginalDistributions = function( codaSamples , data , 
                                          p1Name="demBias" ,
                                          p2Name = "repBias",
                                          datFrm,
                                          showCurve=FALSE ,  
                                          saveName=NULL , saveType="jpg" ) {

  

  # Marginal histograms for regression coefficients
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  responseFirstPartyBias = mcmcMat[,"responseFirstPartyBias"]
  responseUndecideds  = mcmcMat[,"responseUndecideds"]
  
  openGraph(width=10,height=4)
  par( mar=c(4,1,2,1) )
  layout( matrix( 1:2 , nrow=1 ) )
  
  histInfo = plotPost( responseFirstPartyBias , cex.lab = 1.75 , showCurve=showCurve ,
            xlab="Response to Party 1" , main="Coefficient" )
  histInfo = plotPost( responseUndecideds , cex.lab = 1.75 , showCurve=showCurve ,
            xlab="Response to Undecideds" , main="Coefficient" )
  
  
  # marginal histograms for bias spread
  p1PollSpread = mcmcMat[,"p1PollSpread"]
  p2PollSpread = mcmcMat[,"p2PollSpread"]
  
  xlims = c(floor(min(p1PollSpread,p2PollSpread)),ceiling(max(p1PollSpread,p2PollSpread)))
  
  openGraph(width=4,height=10)
  par( mar=c(4,1,2,1) )
  layout( matrix( 1:2 , nrow=2 ) )
  
  histInfo = plotPost( p1PollSpread , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="Party 1 Spread" , main="Scale", xlim = xlims )
  histInfo = plotPost( p2PollSpread , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="Party 2 Spread" , main="Scale", xlim = xlims )
  
  
  
  # marginal distributions for years
  
  YearLevels = levels(as.factor(datFrm[,"year"]))
  numYears = length(YearLevels)
  
  
  
  openGraph(width=3,height=2)
  par( mar=c(4,1,2,1) )
  layout( matrix( 1:8 , nrow=2 ) )
  
  for (i in 1:numYears){
    histInfo = plotPost( mcmcMat[,paste("yearLean[",i,"]",sep = "")], cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=paste(YearLevels[i]," Lean",sep = "") , xlim = xlims )
  }
  histInfo = plotPost( mcmcMat[,"yearSpread"], cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="Year Spread" , xlim = xlims )
  
  
}
