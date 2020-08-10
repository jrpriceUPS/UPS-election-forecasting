#Weighted with Just sample size and transparency predictors
# Notes by Jake 8/10/20


source("DBDA2E-utilities.R")

#===============================================================================
genMCMC = function( refFrame ,datFrmPredictor, pollName="poll" , #daysuntilName="daysuntil", 
                    raceIDName="raceID", actualName="actual",
                    transparencyName="transparency",
                    samplesizeName ="samplesize", whichrace,
                    
                    numSavedSteps=50000 , thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to generic xNom,y variable names for model:
  
  poll = as.numeric(datFrmPredictor[,pollName])
  actual = as.numeric(refFrame[,actualName])
  raceID = as.numeric(as.factor(refFrame[,raceIDName]))
  
  
  transparency = as.numeric(as.factor(datFrmPredictor[,transparencyName]))
  transparencylevels = levels(as.factor(datFrmPredictor[,transparencyName]))
  
  
  
  
  samplesize = as.numeric(datFrmPredictor[,samplesizeName])
  
  
  pollTotal = length(poll)
  NraceIDLvl = length(raceID)
  NtransparencyLvl = length(unique(transparency))
  
  
  
  
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # lmInfo = lm( datFrm[,scoreName] ~ datFrm[,samplesizeName] + 
  #                datFrm[,transparencyName] + datFrm[,transparencyName] + datFrm[,LVName] + datFrm[,transparencyName])
  # residSD = sqrt(mean(lmInfo$residuals^2)) # residual root mean squared deviation
  # For hyper-prior on deflections:
  agammaShRa = unlist( gammaShRaFromModeSD( mode=sd(actual)/2 , sd=2*sd(actual) ) )
  agammaShRasamplesizeImpact = unlist( gammaShRaFromModeSD( mode=sd(samplesize)/2 , sd=2*sd(samplesize) ) )
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    actual=actual,
    poll=poll ,
    whichrace=whichrace,
    raceID=raceID,
    transparency = transparency,
    
    
    samplesize = samplesize,
    
    NtransparencyLvl = NtransparencyLvl ,
    
    
    NraceIDLvl=NraceIDLvl,
    # data properties for scaling the prior:
    samplesizeSD = sd(samplesize) ,
    actualSD = sd(actual) ,
    agammaShRa = agammaShRa,
    agammaShRasamplesizeImpact = agammaShRasamplesizeImpact
    
  )
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  model {
  
  for ( race1 in 1:NraceIDLvl ){
    actual[race1] ~ dnorm(mu[race1], 1/actualSpread^2)
    
   for(myPoll in (whichrace[race1]+1):whichrace[race1+1]){
     weight[myPoll]=transparencyImpact[transparency[myPoll]]+
    samplesizeImpact*samplesize[myPoll]
   
   }

    
    mu[race1] <- sum(nWeight[(whichrace[race1]+1):whichrace[race1+1]]*poll[(whichrace[race1]+1):whichrace[race1+1]])
    summedWeights[race1] = sum(weight[(whichrace[race1]+1):whichrace[race1+1]])
    
    for(myPoll1 in (whichrace[race1]+1):whichrace[race1+1]){
      nWeight[myPoll1]=weight[myPoll1]/summedWeights[race1]
    }
  }


  
    for ( mytransparency in 1:NtransparencyLvl ) { transparencyImpact[mytransparency] ~ dnorm( 0.0 , 1/transparencySpread^2 ) }
    transparencySpread ~ dgamma( agammaShRa[1] , agammaShRa[2] ) 
    
   
    
    
    samplesizeImpact ~ dgamma( agammaShRasamplesizeImpact[1] , agammaShRasamplesizeImpact[2] ) 
    actualSpread ~ dunif(actualSD/100, actualSD*10)
    
    
  }
  " # close quote for modelstring
  writeLines(modelstring,con="TEMPmodel.txt")
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # initsList = list(
  #   a = unname( c(
  #     lmInfo$coef["(Intercept)"] 
  #     + mean(datFrm[,xMetName]) * lmInfo$coef["datFrm[, xMetName]"] ,
  #     lmInfo$coef["(Intercept)"] 
  #     + mean(datFrm[,xMetName]) * lmInfo$coef["datFrm[, xMetName]"] 
  #     + lmInfo$coef[grep( "xNomName" , names(lmInfo$coef) )] ) ) - mean(y) ,
  #   aSigma = sd( c(
  #     lmInfo$coef["(Intercept)"] 
  #     + mean(datFrm[,xMetName]) * lmInfo$coef["datFrm[, xMetName]"] ,
  #     lmInfo$coef["(Intercept)"] 
  #     + mean(datFrm[,xMetName]) * lmInfo$coef["datFrm[, xMetName]"] 
  #     + lmInfo$coef[grep( "xNomName" , names(lmInfo$coef) )] ) ) ,
  #   ySigma = residSD ,
  #   aMet = unname(lmInfo$coefficients["datFrm[, xMetName]"]) ,
  #   nu = 1
  #   # Let JAGS do other parameters automatically...
  # )
  #show( initsList ) 
  #------------------------------------------------------------------------------
  # RUN THE CHAINS
  
  parameters = c(  "transparencyImpact" , "samplesizeImpact"  ,  "actualSpread", "mu"  )
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
  
  #   nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
  #   # Create, initialize, and adapt the model:
  #   jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , inits=initsList , 
  #                           n.chains=nChains , n.adapt=adaptSteps )
  #   # Burn-in:
  #   cat( "Burning in the MCMC chain...\n" )
  #   update( jagsModel , n.iter=burnInSteps )
  #   # The saved MCMC chain:
  #   cat( "Sampling final MCMC chain...\n" )
  #   codaSamples = coda.samples( jagsModel , variable.names=parameters , 
  #                               n.iter=nIter , thin=thinSteps )
  
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
}

#===============================================================================


smryMCMC = function(  codaSamples , datFrm=NULL , transparencyName="transparency" , LVName="LV" , 
                      
                      samplesizeName ="samplesize",
                      #contrasts=NULL ,
                      saveName=NULL ) {
  # All single parameters:
  parameterNames = varnames(codaSamples) 
  if ( !is.null(datFrm) & !is.null(transparencyName) ) {
    transparencylevels = levels(as.factor(datFrm[,transparencyName]))
  }
  
  
  if ( !is.null(datFrm) & !is.null(samplesizeName) ) {
    samplesizelevels = levels(as.factor(datFrm[,samplesizeName]))
  }
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName] ) )
    thisRowName = parName
    if ( !is.null(datFrm) & !is.null(transparencyName) ) {
      # For row name, extract numeric digits from parameter name. E.g., if
      # parameter name is "beta[12,34]" then pull out 12 and 34:
      levelVal = as.numeric( 
        grep( "^[1-9]" , # grep only substrings that begin with digits.
              # Return sll substrings split by "[" or "," or "]":
              unlist( strsplit( parName , "\\[|,|\\]"  ) ) , 
              value=TRUE ) )
      if ( length(levelVal) > 0 ) { 
        
        # Assumes there is only a single factor, i.e., levelVal has only entry: 
        # thisRowName = paste(thisRowName,transparencylevels[levelVal]) 
      }
    }
    rownames(summaryInfo)[NROW(summaryInfo)] = thisRowName
  }
  # # All contrasts:
  # if ( !is.null(contrasts) ) {
  #   if ( is.null(datFrm) | is.null(xNomName) ) {
  #     show(" *** YOU MUST SPECIFY THE DATA FILE AND FACTOR NAMES TO DO CONTRASTS. ***\n")
  #   } else {
  #     # contrasts:
  #     if ( !is.null(contrasts) ) {
  #       for ( cIdx in 1:length(contrasts) ) {
  #         thisContrast = contrasts[[cIdx]]
  #         left = right = rep(FALSE,length(xNomlevels))
  #         for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
  #           left = left | xNomlevels==thisContrast[[1]][nIdx]
  #         }
  #         left = normalize(left)
  #         for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
  #           right = right | xNomlevels==thisContrast[[2]][nIdx]
  #         }
  #         right = normalize(right)
  #         contrastCoef = matrix( left-right , ncol=1 )
  #         postContrast = ( mcmcMat[,paste("aMet[",1:length(xNomlevels),"]",sep="")] 
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
  
  # remove the extra columns
  summaryInfo= summaryInfo[,1:7]
  # Save results:
  
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotPosteriorPredictive = function( codaSamples, refFrame ,datFrmPredictor, pollName="cand1_pct" , 
                                    raceIDName="races", actualName="actual",
                                    whichrace,
                                    saveName=NULL , saveType="jpg", raceplots,
                                    showCurve = FALSE) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  
  actual = refFrame[,actualName]
  polls= datFrmPredictor[,pollName]
  raceID = as.numeric(as.factor(refFrame[,raceIDName]))
  NraceIDLvl = length(raceID)
  
  
  for ( raceidx in raceplots ){
    openGraph(width=8,height=8)
    raceNameidx=refFrame[raceidx,2]
    pollresults = polls[(whichrace[raceidx]+1):whichrace[raceidx+1]]
    plot(actual[raceidx],-.1, xlim = c(floor(min(c(pollresults,actual[raceidx]))/10)*10,ceiling(max(c(pollresults,actual[raceidx]))/10)*10),
         ylim=c(-2,3), cex=2, pch=8, col="steelblue", 
         ylab="Posterior Density for Actual Result", xlab="Dem. Voting Share", main=raceNameidx)
    abline(a=0,b=0)
    
    points(pollresults, runif(length(pollresults))-1.5, col="blue")
    chainSub = round(seq(1,chainLength,length=20))
    
    
    
    for ( chnIdx in chainSub ) {
      m = mcmcMat[chnIdx,paste("mu[",raceidx,"]",sep="")]
      
      s = mcmcMat[chnIdx,"actualSpread"] # spread
      
      
      nlim = qnorm( c(0.01,0.99) )
      yl = m+nlim[1]*s
      yh = m+nlim[2]*s
      ycomb=seq(yl,yh,length=201)
      ynorm = dnorm(ycomb,mean=m,sd=s)
      ynorm = 2.75*ynorm/max(ynorm)
      
      lines( ycomb , ynorm , col="skyblue" )
      
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PostPred-",raceNameidx), type=saveType)
    }
  }
}

meanPosterior = function (codaSamples, refFrame ,datFrmPredictor, pollName="cand1_pct" , 
                          raceIDName="races", actualName="actual",
                          whichrace,
                          saveName=NULL , saveType="jpg", raceplots,
                          showCurve = FALSE) {
  
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  
  actual = refFrame[,actualName]
  polls= datFrmPredictor[,pollName]
  raceID = as.numeric(as.factor(refFrame[,raceIDName]))
  NraceIDLvl = length(raceID)
  
  
  for ( raceidx in raceplots ){
    openGraph(width=13,height=8)
    #Arrange plots to 2 rows and 1 column. 
    layout(matrix(c(1,2), ncol=1))
    raceNameidx=refFrame[raceidx,2]
    pollresults = polls[(whichrace[raceidx]+1):whichrace[raceidx+1]]
    
    #plot the mean
    plotPost( mcmcMat[,paste("mu[",raceidx,"]",sep="")], cex.lab = 1.75 , showCurve=showCurve ,
              xlab=bquote(LVImpact) , main="Mean",xlim = c(floor(min(c(pollresults,actual[raceidx]))/10)*10,ceiling(max(c(pollresults,actual[raceidx]))/10)*10)
    )
    
    #plot post-predictive
    plot(actual[raceidx],-.1, xlim = c(floor(min(c(pollresults,actual[raceidx]))/10)*10,ceiling(max(c(pollresults,actual[raceidx]))/10)*10),
         ylim=c(-2,3), cex=2, pch=8, col="steelblue", 
         ylab="Posterior Density for Actual Result", xlab="Dem. Voting Share", main=raceNameidx)
    abline(a=0,b=0)
    
    points(pollresults, runif(length(pollresults))-1.5, col="blue")
    chainSub = round(seq(1,chainLength,length=20))
    
    
    
    for ( chnIdx in chainSub ) {
      m = mcmcMat[chnIdx,paste("mu[",raceidx,"]",sep="")]
      
      s = mcmcMat[chnIdx,"actualSpread"] # spread
      
      
      nlim = qnorm( c(0.01,0.99) )
      yl = m+nlim[1]*s
      yh = m+nlim[2]*s
      ycomb=seq(yl,yh,length=201)
      ynorm = dnorm(ycomb,mean=m,sd=s)
      ynorm = 2.75*ynorm/max(ynorm)
      
      lines( ycomb , ynorm , col="skyblue" )
      
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PostMean-",raceNameidx), type=saveType)
    }
  }
}


plotSampleSizePosterior = function( codaSamples , 
                                    datFrm  , 
                                    saveName=NULL , saveType="jpg",
                                    showCurve = FALSE, title="Sample Size Impact") {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  samplesizeName ="samplesize"
  pollName ="poll"
  
  openGraph(width=8,height=8)
  
  # posterior of the mean for sample size distrubtion 
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  samplesizeImpact= mcmcMat[,"samplesizeImpact"]
  plotPost( samplesizeImpact , cex.lab = 1.75 , showCurve=showCurve ,
            xlab=bquote(samplesizeImpact) , main=title )
  print(mean(samplesizeImpact))
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,"samplesizeImpact", type=saveType))
    
  }
}



plotTransparencyPosterior = function( codaSamples , 
                                      datFrmPredictor  , 
                                      saveName=NULL , saveType="jpg",
                                      showCurve = FALSE) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  transparencyName ="transparency"
  pollName ="poll"
  
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  transparency = as.numeric(as.factor(datFrmPredictor[,transparencyName]))
  transparencyLevels = levels(as.factor(datFrmPredictor[,transparencyName]))
  NtransparencyLvl = length(unique(transparency))
  for ( transparencyidx in 1:length(transparencyLevels)) {
    openGraph(width=4,height=4)
    
    # posterior of the mean for sample size distrubtion 
    
    #give better titles - using delModLevels
    if(transparencyidx==1){titleT="Untransparent Impact"}
    if(transparencyidx==2){titleT="Transparent Impact"}
    
    
    
    
    plotPost( mcmcMat[,paste("transparencyImpact[",transparencyidx,"]",sep="")], cex.lab = 1.75 , showCurve=showCurve ,
              xlab=bquote(transparencyImpact) , main=titleT)
    
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,titleT, type=saveType))
      
    }
  }
}





