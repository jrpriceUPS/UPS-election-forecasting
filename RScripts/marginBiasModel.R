# Author: Jake Price
# Date: 8/25/20
#
# A set of functions for running the margin bias of a pollster in a given year ("house effect")
# simulation and making plots
#  * Overall year lean pulled from parent distribution of year leans
#  * Pollster bias in a given year comes from parent distribution for that pollster
#     - Pollster parent dist. comes from larger parent distribution of pollster centers


# source summary functions (look at this later and pull in only the functions we need)
source("DBDA2E-utilities.R")

#===============================================================================
runSimulation = function( datFrm , numSavedSteps=50000 , thinSteps=1 , saveName=NULL,
                          runjagsMethod=runjagsMethodDefault , nChains=nChainsDefault ) 
{ 
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to lists for delivery to JAGS
 
  marginBias = as.numeric(datFrm[,"bias"])
  
  pollster = as.numeric(as.factor(datFrm[,"pollster"]))
  PollsterLevels = levels(as.factor(datFrm[,"pollster"]))
  year = as.numeric(as.factor(datFrm[,"year"]))
  YearLevels = levels(as.factor(datFrm[,"year"]))
  
  # indexing
  PollsTotal = length(marginBias)
  PollsterLevelsTotal = length(unique(pollster))
  YearLevelsTotal = length(unique(year))
  
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # For prior on baseline, etc.:
  marginMean = mean(marginBias)
  marginSD = sd(marginBias)
  
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    marginBias = marginBias,
    pollster = pollster,
    year = year,
    
    # indexing
    PollsTotal = PollsTotal ,
    PollsterLevelsTotal = PollsterLevelsTotal ,
    YearLevelsTotal =  YearLevelsTotal ,
    
    # data properties for scaling the prior:
    marginMean = marginMean,
    marginSD = marginSD
  )
  
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  model {
  
    #Bottom Level (individual poll bias):
    
    for (poll in 1:PollsTotal){
    
      # margin is modeled as t-distribution with center dictated by year lean + house effects
      mu[poll] <- yearLean[year[poll]] + pollsterBias[pollster[poll],year[poll]] 
      marginBias[poll] ~ dt(mu[poll] , (1/marginSpread^2) , nu )
    
    }
    
    # vague prior for normality of margin bias
    nu ~  dexp(1/30.0) 
    
    # vague priors on appropriate scale for spread of poll biases
    marginSpread ~ dunif( marginSD/100 , marginSD*10 )
    

   
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
  parameters = c("pollsterBias","pollsterCenter","pollsterSpread","yearSpread","yearLean",
                 "marginSpread","nu")
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
plotDiagnostics= function( codaObject ){
  for ( parName in c("marginSpread",  "nu"  , "yearSpread" , "yearLean[1]", "pollsterBias[1,1]", "pollsterCenter[1]", "pollsterSpread[1]"  ) ) {
    diagMCMC( codaObject=codaObject , parName=parName)
  }
}
#===============================================================================

#===============================================================================

plotPosteriorPredictive = function( codaSamples , 
                                          datFrm ,
                                          whichPollsters,
                                          saveName=NULL , 
                                          saveType="jpg",
                                          showCurve = FALSE) {
  myCol = "purple"
  myLabel = "Democratic Margin Bias (Dem - Rep)"
  
  
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  marginBias = datFrm[,"bias"]
  year = as.numeric(as.factor(datFrm[,"year"]))
  YearLevels = levels(as.factor(datFrm[,"year"]))
  
  # Display data with posterior predictive distributions
  for ( Yearidx in 1:length(YearLevels) ) {
    openGraph(width=2*length(whichPollsters),height=5)
    par( mar=c(4,4,2,1) , mgp=c(3,1,0) )
    
    limits = c(min(marginBias)-0.2*(max(marginBias)-min(marginBias)),max(marginBias)+0.2*(max(marginBias)-min(marginBias))) 
    
    plot(-10,-10,
         xlim=c(0.2,length(whichPollsters)+0.1), xlab = "",
         xaxt="n" , ylab = myLabel ,
         
         ylim=limits ,
         main="Data with Post. Pred.")
    axis( 1 , at=1:length(whichPollsters) , tick=FALSE ,
          lab=paste( whichPollsters , YearLevels[Yearidx] , sep="\n" ) )
    for ( i in 1:length(whichPollsters) ) {
      pollsterIdx = which(levels(as.factor(datFrm$pollster))==whichPollsters[i])
      
      
      xPlotVal = i
      y1Vals = marginBias[ as.numeric(as.factor(datFrm$pollster))==pollsterIdx & year==Yearidx ]
      points( rep(xPlotVal,length(y1Vals))+runif(length(y1Vals),-0.05,0.05) ,
              y1Vals , pch=1 , cex=1.5 , col=myCol )
      chainSub = round(seq(1,chainLength,length=20))
      for ( chnIdx in chainSub ) {
        m = mcmcMat[chnIdx,paste("pollsterBias[",pollsterIdx,",", Yearidx,"]",sep="")]   # pollster bias
        +         mcmcMat[chnIdx,paste("yearLean[",Yearidx,"]",sep="")]  # year lean
        s = mcmcMat[chnIdx,"marginSpread"] # spread
        nu = mcmcMat[chnIdx,"nu"]# normality
        
        
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
      saveGraph( file=paste0(saveName,"PostPred-margin-",YearLevels[Yearidx]), type=saveType)
    }
  }# end for Yearidx
}




#===============================================================================  

#===============================================================================


plotMarginalDistributions = function( codaSamples , data , 
                                      datFrm,
                                      whichPollsters,
                                      showCurve=FALSE ,  
                                      saveName=NULL , saveType="jpg" ) {
  
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  
  
  # Marginal histograms for regression coefficients
  coeffMarginal(mcmcMat, showCurve, saveName, saveType)
  
  # marginal distributions for years
  yearMarginal(mcmcMat,datFrm, showCurve, saveName, saveType)
  
  
  
  
  # marginal distributions for pollsters
  pollsterMarginal( mcmcMat, datFrm, whichPollsters, showCurve, saveName, saveType)
  
}


yearMarginal = function( mcmcMat, datFrm, showCurve, saveName, saveType) {
  # marginal distributions for years
  
  YearLevels = levels(as.factor(datFrm[,"year"]))
  numYears = length(YearLevels)
  
  openGraph(width=16,height=8)
  par( mar=c(4,1,2,1) )
  layout( matrix( 1:6 , nrow=2, byrow = TRUE ) )
  
  for (i in 1:numYears){
    histInfo = plotPost( mcmcMat[,paste("yearLean[",i,"]",sep = "")], cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=paste(YearLevels[i]," Lean",sep = ""))
  }
  histInfo = plotPost( mcmcMat[,"yearSpread"], cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="Year Spread" )
  
  saveGraph( file=paste0(saveName,"YearMarginals-"), type=saveType)
  
}

coeffMarginal = function( mcmcMat, showCurve, saveName, saveType ) {
  
  # marginal histograms for bias spread
  marginSpread = mcmcMat[,"marginSpread"]
  nu = mcmcMat[,"nu"]
  
  
  openGraph(width=10,height=4)
  par( mar=c(6,1,2,1) )
  layout( matrix( 1:2 , nrow=1 ) )
  
  histInfo = plotPost( marginSpread , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="Margin Bias (Dem - Rep) Spread" , main="Scale" )
  histInfo = plotPost( nu , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="Margin Bias Normality" , main="nu")
  
  saveGraph( file=paste0(saveName,"SpreadNormalityMarginals-"), type=saveType)
  
}

pollsterMarginal = function( mcmcMat, datFrm, whichPollsters, showCurve, saveName, saveType){
  # create marginal distributions for each pollster requested
  
  YearLevels = levels(as.factor(datFrm[,"year"]))
  numYears = length(YearLevels)
  
  
  
  for (i in 1:length(whichPollsters)){
    
    openGraph(width=20,height=8)
    par( mar=c(4,1,2,1) )
    layout( matrix( 1:8 , nrow=2, byrow = TRUE ) )
    
    pollsterIdx = which(levels(as.factor(datFrm$pollster))==whichPollsters[i])
    fullName = datFrm$pollsterFullName[which(datFrm$pollster==whichPollsters[i])[1]]
    
    pollsterCenter = mcmcMat[,paste("pollsterCenter[",pollsterIdx,"]",sep = "")]
    pollsterSpread = mcmcMat[,paste("pollsterSpread[",pollsterIdx,"]",sep = "")]
    
    histInfo = plotPost( pollsterCenter, cex.lab = 1.75 , showCurve=showCurve ,
                         xlab="Average House Effect", main = fullName )
    histInfo = plotPost( pollsterSpread, cex.lab = 1.75 , showCurve=showCurve ,
                         xlab="House Effect Variability", main = fullName )
    
    for (yearIdx in 1:numYears) {
      pollsterBias = mcmcMat[,paste("pollsterBias[",pollsterIdx,",",yearIdx,"]",sep = "")] +
        mcmcMat[,paste("yearLean[",yearIdx,"]",sep = "")]
      
      histInfo = plotPost( pollsterBias, cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=paste(YearLevels[yearIdx]," House Effect + Year Lean",sep=""), main = fullName )
    }
    
    saveGraph( file=paste0(saveName,whichPollsters[i],"-"), type=saveType)
  }
  
  
  
  
}