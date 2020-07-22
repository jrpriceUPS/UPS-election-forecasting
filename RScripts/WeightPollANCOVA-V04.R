#ANCOVA for Weight of Poll V03
#06/30/2020


# Accompanies the book:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================
genMCMC = function( refFrame ,datFrmPredictor, pollName="poll" , #daysuntilName="daysuntil", 
                    raceIDName="raceID", actualName="actual",
                    delModeName="delMode" , LVName="LV" , transparencyName="transparency", 
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
  
  
  delMode = as.numeric(as.factor(datFrmPredictor[,delModeName]))
  delModelevels = levels(as.factor(datFrmPredictor[,delModeName]))
  LV = as.numeric(as.factor(datFrmPredictor[,LVName]))
  LVlevels = levels(as.factor(datFrmPredictor[,LVName]))
  transparency = as.numeric(as.factor(datFrmPredictor[,transparencyName]))
  transparencylevels = levels(as.factor(datFrmPredictor[,transparencyName]))  
  
  samplesize = as.numeric(datFrmPredictor[,samplesizeName])
  
  
  pollTotal = length(poll)
  NraceIDLvl = length(raceID)
  NdelModeLvl = length(unique(delMode))
  NLVLvl = length(unique(LV))
  NtransparencyLvl = length(unique(transparency))
  
  
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # lmInfo = lm( datFrm[,scoreName] ~ datFrm[,samplesizeName] + 
  #                datFrm[,delModeName] + datFrm[,delModeName] + datFrm[,LVName] + datFrm[,transparencyName])
  # residSD = sqrt(mean(lmInfo$residuals^2)) # residual root mean squared deviation
  # For hyper-prior on deflections:
  agammaShRa = unlist( gammaShRaFromModeSD( mode=sd(actual)/2 , sd=2*sd(actual) ) )
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    actual=actual,
    poll=poll ,
    whichrace=whichrace,
    raceID=raceID,
    delMode = delMode,
    LV = LV,
    transparency = transparency,
    samplesize = samplesize,
    
    NdelModeLvl = NdelModeLvl ,
    NLVLvl = NLVLvl,
    NtransparencyLvl = NtransparencyLvl ,
    NraceIDLvl=NraceIDLvl,
    # data properties for scaling the prior:
    samplesizeSD = sd(samplesize) ,
    actualSD = sd(actual) ,
    agammaShRa = agammaShRa
    
  )
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  model {
  
  for ( race1 in 1:NraceIDLvl ){
    actual[race1] ~ dnorm(mu[race1], 1/actualSpread^2)
    
   for(myPoll in (whichrace[race1]+1):whichrace[race1+1]){
     weight[myPoll]=delModeImpact[delMode[myPoll]]+LVImpact[LV[myPoll]]+
     transparencyImpact[transparency[myPoll]]+samplesizeImpact*samplesize[myPoll]
   
   }

    
    mu[race1] <- sum(nWeight[(whichrace[race1]+1):whichrace[race1+1]]*poll[(whichrace[race1]+1):whichrace[race1+1]])
    summedWeights[race1] ~ sum(weight[])
    for(myPoll1 in (whichrace[race1]+1):whichrace[race1+1]){
      nWeight[myPoll1]=weight[myPoll1]/summedWeights
    }
  }

  
    for ( mydelMode in 1:NdelModeLvl ) { delModeImpact[mydelMode] ~ dnorm( 0.0 , 1/delModeSpread^2 ) 
                               }
   
   
    delModeSpread ~ dgamma( agammaShRa[1] , agammaShRa[2] ) 
    
    for ( myLV in 1:NLVLvl ) { LVImpact[myLV] ~ dnorm( 0.0 , 1/LVSpread^2 ) 
                               }
    LVSpread ~ dgamma( agammaShRa[1] , agammaShRa[2] ) 
    
    for ( mytransparency in 1:NtransparencyLvl ) { transparencyImpact[mytransparency] ~ dnorm( 0.0 , 1/transparencySpread^2 ) 
                              }
    transparencySpread ~ dgamma( agammaShRa[1] , agammaShRa[2] ) 
    
    samplesizeImpact ~ dnorm( 0 , 1/(samplesizeSD*.01)^2 )
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
  
  parameters = c(  "delModeImpact" , "samplesizeImpact" , "LVImpact" , "transparencyImpact", "actualSpread"  )
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

smryMCMC = function(  codaSamples , datFrm=NULL , delModeName="delMode" , LVName="LV" , 
                      transparencyName="transparency", 
                      samplesizeName ="samplesize",
                      #contrasts=NULL ,
                      saveName=NULL ) {
  # All single parameters:
  parameterNames = varnames(codaSamples) 
  if ( !is.null(datFrm) & !is.null(delModeName) ) {
    delModelevels = levels(as.factor(datFrm[,delModeName]))
  }
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName] ) )
    thisRowName = parName
    if ( !is.null(datFrm) & !is.null(delModeName) ) {
      # For row name, extract numeric digits from parameter name. E.g., if
      # parameter name is "beta[12,34]" then pull out 12 and 34:
      levelVal = as.numeric( 
        grep( "^[1-9]" , # grep only substrings that begin with digits.
              # Return sll substrings split by "[" or "," or "]":
              unlist( strsplit( parName , "\\[|,|\\]"  ) ) , 
              value=TRUE ) )
      if ( length(levelVal) > 0 ) { 
        # Assumes there is only a single factor, i.e., levelVal has only entry: 
        thisRowName = paste(thisRowName,delModelevels[levelVal]) 
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
  # Save results:
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , datFrm ,
                     contrasts=NULL , saveName=NULL , saveType="jpg",
                     scoreName="score" , 
                     delModeName="delMode" ,
                     samplesizeName ="samplesize"
                     
) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  score = as.numeric(datFrm[,scoreName])
  
  # daysuntil = as.numeric(datFrm[,daysuntilName])
  delMode = as.numeric(as.factor(datFrm[,delModeName]))
  delModelevels = levels(as.factor(datFrm[,delModeName]))
  # LV = as.numeric(as.factor(datFrm[,LVName]))
  #  LVlevels = levels(as.factor(datFrm[,LVName]))
  # transparency = as.numeric(as.factor(datFrm[,transparencyName]))
  #transparencylevels = levels(as.factor(datFrm[,transparencyName]))  
  
  samplesize = as.numeric(datFrm[,samplesizeName])
  
  
  NdelModeLvl = length(unique(delMode))
  NLVLvl = length(unique(LV))
  NtransparencyLvl = length(unique(transparency))
  
  # xNom = as.numeric(as.factor(datFrm[,xNomName]))
  # xNomlevels = levels(as.factor(datFrm[,xNomName]))
  # xMet = as.numeric(datFrm[,xMetName])
  # Ntotal = length(y)
  # NxNomLvl = length(unique(xNom))
  
  # Display data with posterior predictive distributions:
  for ( xNomLvlIdx in 1:NdelModeLvl ) {
    # Open blank graph with appropriate limits:
    xLim = c( min(samplesize)-0.2*(max(samplesize)-min(samplesize)) ,
              max(samplesize)+0.2*(max(samplesize)-min(samplesize)) )
    yLim = c( min(score)-0.2*(max(score)-min(score)) , 
              max(score)+0.2*(max(score)-min(score)) )
    openGraph(width=4,height=5)
    par(mar=c(3,3,3,0.5)) # number of margin lines: bottom,left,top,right
    par(mgp=c(1.75,0.5,0)) # which margin lines to use for labels
    plot(2*max(xLim),2*max(yLim), # point out of range not seen 
         xlab=samplesizeName , xlim=xLim , ylab=yName , ylim=yLim , 
         main=paste(NLVLvl[xNomLvlIdx],"Data\nwith Post. Pred. Distrib.") ) 
    # plot credible regression lines and noise profiles:
    nSlice = 3
    curveXpos = seq(min(samplesize),max(samplesize),length=nSlice)
    curveWidth = (max(samplesize)-min(samplesize))/(nSlice+2)
    nPredCurves=30
    for ( i in floor(seq(from=1,to=nrow(mcmcMat),length=nPredCurves)) ) {
      intercept = mcmcMat[i,paste0("modeImpact[",xNomLvlIdx,"]")]
      slope = mcmcMat[i,paste0("modeImpact[delMode[",xNomLvlIdx,"]]")]
      noise = mcmcMat[i,"scoreSpread"]
      abline( a=intercept , b=slope , col="skyblue" )
      for ( j in 1:nSlice ) {
        hdiLo = intercept+slope*curveXpos[j] - 1.96*noise
        hdiHi = intercept+slope*curveXpos[j] + 1.96*noise
        yComb = seq( hdiLo , hdiHi , length=75 )
        xVals = dnorm( yComb , mean=intercept+slope*curveXpos[j] , sd=noise ) 
        xVals = curveWidth * xVals / max(xVals)
        lines( curveXpos[j] - xVals , yComb , col="skyblue" )
        lines( curveXpos[j] - 0*xVals , yComb , col="skyblue" , lwd=2 )
      }
    }
    # plot data points:
    includeVec = ( delMode == xNomLvlIdx )
    xVals = samplesize[includeVec]
    yVals = score[includeVec]
    points( xVals , yVals , pch=1 , cex=1.5 , col="red" )
    
    
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PostPred-",xNomlevels[xNomLvlIdx]), 
                 type=saveType)
    }
  }
  # 
  # # Display contrast posterior distributions:
  # if ( !is.null(contrasts) ) {
  #   if ( is.null(datFrm) | is.null(xNomName) ) {
  #     show(" *** YOU MUST SPECIFY THE DATA FILE AND FACTOR NAMES TO DO CONTRASTS. ***\n")
  #   } else {
  #     for ( cIdx in 1:length(contrasts) ) {
  #       thisContrast = contrasts[[cIdx]]
  #       left = right = rep(FALSE,length(xNomlevels))
  #       for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
  #         left = left | xNomlevels==thisContrast[[1]][nIdx]
  #       }
  #       left = normalize(left)
  #       for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
  #         right = right | xNomlevels==thisContrast[[2]][nIdx]
  #       }
  #       right = normalize(right)
  #       contrastCoef = matrix( left-right , ncol=1 )
  #       postContrast = ( mcmcMat[,paste("aMet[",1:length(xNomlevels),"]",sep="")] 
  #                        %*% contrastCoef )
  #       openGraph(height=8,width=4)
  #       layout(matrix(1:2,ncol=1))
  #       plotPost( postContrast , xlab="Difference" ,
  #                 main=paste0( 
  #                   paste(thisContrast[[1]],collapse="."), 
  #                   "\nvs\n",
  #                   paste(thisContrast[[2]],collapse=".") ) ,
  #                 compVal=thisContrast$compVal , ROPE=thisContrast$ROPE )
  #       plotPost( postContrast/mcmcMat[,"ySigma"] , 
  #                 xlab="Effect Size" ,
  #                 main=paste0( 
  #                   paste(thisContrast[[1]],collapse="."), 
  #                   "\nvs\n",
  #                   paste(thisContrast[[2]],collapse=".") ) ,
  #                 compVal=0.0 , 
  #                 ROPE=c(-0.1,0.1) )
  #       
  #       if ( !is.null(saveName) ) {
  #         saveGraph( file=paste0(saveName, paste0( 
  #           paste(thisContrast[[1]],collapse=""), 
  #           ".v.",
  #           paste(thisContrast[[2]],collapse="") ) ), 
  #           type=saveType )
  #       }
  #     }
  #   }
  # } # end if ( !is.null(contrasts) )
}
