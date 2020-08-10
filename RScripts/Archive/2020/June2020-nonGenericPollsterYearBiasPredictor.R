# non-generic version of JAGS script to predict bias using pollster and year
# Notes by Jake 8/10/20


#JAGS-2Factor-V03
#06-29-2020
#Modification of JAGS-2Factor-V02.R 
  #Change Parameter names to be non-generic

source("DBDA2E-utilities.R")

#===============================================================================
genMCMC = function( datFrm , biasName = "bias" , pollsterName = "pollster" , yearName = "year",
                    numSavedSteps=50000 , thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to generic x,y variable names for model:
  bias = as.numeric(datFrm[,biasName])
  pollster = as.numeric(as.factor(datFrm[,pollsterName]))
  PollsterLevels = levels(as.factor(datFrm[,pollsterName]))
  year = as.numeric(as.factor(datFrm[,yearName]))
  YearLevels = levels(as.factor(datFrm[,yearName]))
  PollsTotal = length(bias)
  PollsterLevelsTotal = length(unique(pollster))
  YearLevelsTotal = length(unique(year))
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # For prior on baseline, etc.:
  biasMean = mean(bias)
  biasSD = sd(bias)
  
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    bias = bias ,
    pollster = pollster ,
    year = year ,
    PollsTotal = PollsTotal ,
    PollsterLevelsTotal = PollsterLevelsTotal ,
    YearLevelsTotal =  YearLevelsTotal ,
    # data properties for scaling the prior:
    biasMean = biasMean ,
    biasSD = biasSD
  )
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  model {
  
    # Bottom level of the hierarchy (individual poll bias)
    for ( poll in 1:PollsTotal) {
      bias[poll] ~ dt(mu[poll], 1/biasSpread^2, nuY )
      mu[poll] <- yearLean[year[poll]] + pollsterBias[pollster[poll],year[poll]] 
    }
    nuY ~  dexp(1/30.0) 
    biasSpread ~ dunif( biasSD/100 , biasSD*10 )
  
  
    # Middle level of the hierarchy (year lean)
    for ( year in 1:YearLevelsTotal ) { yearLean[year] ~ dnorm( 0.0 , 1/yearSpread^2) }
    yearSpread ~ dunif( biasSD/100 , biasSD*10 )
    
    # Middle level of the hierarchy (pollster bias with year)
     for ( pollster in 1:PollsterLevelsTotal ) { for ( year in 1:YearLevelsTotal ) { 
     pollsterBias[pollster,year] ~ dnorm(0.0 , 1/pollsterSpread^2) }
     
    
     }
     pollsterSpread ~ dunif( biasSD/100 , biasSD*10 )
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
  parameters = c( "biasSpread" , "nuY" , "pollsterSpread", "yearSpread", "yearLean", "pollsterBias" )
  adaptSteps = 500 
  burnInSteps = 1000 
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          #inits=initsList , 
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

smryMCMC = function(  codaSamples , datFrm=NULL , pollsterName=NULL , yearName = NULL,
                      contrasts=NULL , saveName=NULL ) {
  # All single parameters:
  parameterNames = varnames(codaSamples) 
  if ( !is.null(datFrm) & !is.null(pollsterName) & !is.null(yearName) ) {
    PollsterLevels = levels(as.factor(datFrm[,pollsterName]))
    YearLevels = levels(as.factor(datFrm[,yearName]))
  }
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName] ) )
    thisRowName = parName
    if ( !is.null(datFrm) & !is.null(pollsterName) & !is.null(yearName) ) {
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
  
  # All contrasts:
  # if ( !is.null(contrasts) ) {
  #   if ( is.null(datFrm) | is.null(pollsterName) ) {
  #     show(" *** YOU MUST SPECIFY THE DATA FILE AND FACTOR NAMES TO DO CONTRASTS. ***\n")
  #   } else {
  # contrasts:
  #     if ( !is.null(contrasts) ) {
  #       for ( cIdx in 1:length(contrasts) ) {
  #         thisContrast = contrasts[[cIdx]]
  #         left = right = rep(FALSE,length(xlevels))
  #         for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
  #           left = left | xlevels==thisContrast[[1]][nIdx]
  #         }
  #         left = normalize(left)
  #         for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
  #           right = right | xlevels==thisContrast[[2]][nIdx]
  #         }
  #         right = normalize(right)
  #         contrastCoef = matrix( left-right , ncol=1 )
  #         postContrast = ( mcmcMat[,paste("b[",1:length(xlevels),"]",sep="")] 
  #                          %*% contrastCoef )
  #         summaryInfo = rbind( summaryInfo , 
  #                              summarizePost( postContrast ,
  #                                             compVal=thisContrast$compVal ,
  #                                             ROPE=thisContrast$ROPE ) )
  #         rownames(summaryInfo)[NROW(summaryInfo)] = (
  #           paste( paste(thisContrast[[1]],collapse=""), ".v.",
  #                  paste(thisContrast[[2]],collapse=""),sep="") )
  #       }
  #     }
  #   }
  # }
  
  # Save results:
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotDiagnostics= function( ){
for ( parName in c("biasSpread",  "nuY" , "pollsterSpread" , "yearSpread" , "yearLean[1]",
                   "pollsterBias[1,1]"  ) ) {
  diagMCMC( codaObject=mcmcCodaV03 , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
}
#===============================================================================

plotPosteriorPredictive = function( codaSamples , 
                     datFrm , biasName=NULL , 
                     saveName=NULL , saveType="jpg",
                     showCurve = FALSE) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  pollsterName="pollster"
  yearName="year"
  bias = datFrm[,biasName]
  pollster = as.numeric(as.factor(datFrm[,pollsterName]))
  PollsterLevels = levels(as.factor(datFrm[,pollsterName]))
  year = as.numeric(as.factor(datFrm[,yearName]))
  YearLevels = levels(as.factor(datFrm[,yearName]))
 # Display data with posterior predictive distributions
    for ( Yearidx in 1:length(YearLevels) ) {
      openGraph(width=2*length(PollsterLevels),height=5)
      par( mar=c(4,4,2,1) , mgp=c(3,1,0) )
      plot(-10,-10,
           xlim=c(0.2,length(PollsterLevels)+0.1) ,
           xlab=paste(pollsterName,yearName,sep="\n") ,
           xaxt="n" , ylab=biasName ,
           ylim=c(min(bias)-0.2*(max(bias)-min(bias)),max(bias)+0.2*(max(bias)-min(bias))) ,
           main="Data with Post. Pred.")
      axis( 1 , at=1:length(PollsterLevels) , tick=FALSE ,
            lab=paste( PollsterLevels , YearLevels[Yearidx] , sep="\n" ) )
      for ( Pollsteridx in 1:length(PollsterLevels) ) {
        xPlotVal = Pollsteridx #+ (Yearidx-1)*length(PollsterLevels)
        yVals = bias[ pollster==Pollsteridx & year==Yearidx ]
        points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) ,
                yVals , pch=1 , cex=1.5 , col="red" )
        chainSub = round(seq(1,chainLength,length=20))
        for ( chnIdx in chainSub ) {
          m = mcmcMat[chnIdx,paste("pollsterBias[",Pollsteridx,",", Yearidx,"]",sep="")]   # pollster bias
  +         mcmcMat[chnIdx,paste("yearLean[",Yearidx,"]",sep="")]  # year lean
          s = mcmcMat[chnIdx,"biasSpread"] # spread
          nu = mcmcMat[chnIdx,"nuY"]# normality


          tlim = qt( c(0.025,0.975) , df=nu )
          yl = m+tlim[1]*s
          yh = m+tlim[2]*s
          ycomb=seq(yl,yh,length=201)
          #ynorm = dnorm(ycomb,mean=m,sd=s)
          #ynorm = 0.67*ynorm/max(ynorm)
          yt = dt( (ycomb-m)/s , df=nu )
          yt = 0.67*yt/max(yt)
          lines( xPlotVal-yt , ycomb , col="skyblue" )



        }
      }
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,"PostPred-",YearLevels[Yearidx]), type=saveType)
      }
    }# end for Yearidx

}
  
  
  
  #===============================================================================  
  
plotYearPosterior = function( codaSamples , 
                                    datFrm , biasName=NULL , 
                                    saveName=NULL , saveType="jpg",
                                    showCurve = FALSE) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  pollsterName="pollster"
  yearName="year"
  bias = datFrm[,biasName]
  pollster = as.numeric(as.factor(datFrm[,pollsterName]))
  PollsterLevels = levels(as.factor(datFrm[,pollsterName]))
  year = as.numeric(as.factor(datFrm[,yearName]))
  YearLevels = levels(as.factor(datFrm[,yearName]))
   #plot each year
  for ( Yearidx in 1:length(YearLevels) ) {
    openGraph(width=8,height=8)
    
    # posterior of the mean for that pollster
    histInfo = plotPost( mcmcMat[,paste("yearLean[",Yearidx,"]",sep="")] , cex.lab = 1.75 , showCurve=showCurve ,
                         #compVal=compValMu , ROPE=ropeMu ,
                         xlab=YearLevels[Yearidx] , main=paste("Mean") ,
                         col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"YearLean-",YearLevels[Yearidx]), type=saveType)
    }
  }
  
}
#===============================================================================
plotPollsterPosterior = function( codaSamples , 
                                    datFrm , biasName=NULL ,
                                    saveName=NULL , saveType="jpg",
                                    showCurve = FALSE) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  pollsterName="pollster"
  yearName="year"
  bias = datFrm[,biasName]
  pollster = as.numeric(as.factor(datFrm[,pollsterName]))
  PollsterLevels = levels(as.factor(datFrm[,pollsterName]))
  year = as.numeric(as.factor(datFrm[,yearName]))
  YearLevels = levels(as.factor(datFrm[,yearName]))
 #layout(matrix(1:ceiling(length(YearLevels)/4)*4,nrow=4))
  #plot each pollster
  for (Pollsteridx in 1:length(PollsterLevels)) { 
    for (Yearidx in 1:length(YearLevels) )     {
    
      openGraph(width=8,height=8)
    
    # posterior of the mean for that pollster

    histInfo = plotPost(mcmcMat[,paste("pollsterBias[",Pollsteridx,",", Yearidx,"]",sep="")] , cex.lab = 1.75 , showCurve=showCurve ,
                        #compVal=compValMu , ROPE=ropeMu ,
                        xlab= PollsterLevels[Pollsteridx]  , main= YearLevels[Yearidx]  ,
                        col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PollsterBias-",PollsterLevels[Pollsteridx],YearLevels[Yearidx]), type=saveType)
    }
  }
 }
}
