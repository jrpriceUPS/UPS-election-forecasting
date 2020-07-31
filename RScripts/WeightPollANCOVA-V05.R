#ANCOVA for Weight of Poll V04
#06/30/2020


# Accompanies the book:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================
genMCMC = function( refFrame ,datFrmPredictor, pollName="poll" , #daysuntilName="daysuntil", 
                    raceIDName="raceID", actualName="actual",
                    IVRName="IVR" , onlineName="online", liveName="live", textName="text",
                    LVName="LV" , transparencyName="transparency", 
                    samplesizeName ="samplesize", whichrace, yearName="year", pollsterName="pollster",
                    biasName="bias",
                    
                    numSavedSteps=50000 , thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to generic xNom,y variable names for model:
  
  poll = as.numeric(datFrmPredictor[,pollName])
  actual = as.numeric(refFrame[,actualName])
  raceID = as.numeric(as.factor(refFrame[,raceIDName]))
  
  
  IVR = as.numeric(as.factor(datFrmPredictor[,IVRName]))
  IVRlevels = levels(as.factor(datFrmPredictor[,IVRName]))
  online = as.numeric(as.factor(datFrmPredictor[,onlineName]))
  onlinelevels = levels(as.factor(datFrmPredictor[,onlineName]))
  live = as.numeric(as.factor(datFrmPredictor[,liveName]))
  livelevels = levels(as.factor(datFrmPredictor[,liveName]))
  text = as.numeric(as.factor(datFrmPredictor[,textName]))
  textlevels = levels(as.factor(datFrmPredictor[,textName]))
  
  
  LV = as.numeric(as.factor(datFrmPredictor[,LVName]))
  LVlevels = levels(as.factor(datFrmPredictor[,LVName]))
  transparency = as.numeric(as.factor(datFrmPredictor[,transparencyName]))
  transparencylevels = levels(as.factor(datFrmPredictor[,transparencyName]))  
  
  samplesize = as.numeric(datFrmPredictor[,samplesizeName])
  
  
  pollTotal = length(poll)
  NraceIDLvl = length(raceID)
  NIVRLvl = length(unique(IVR))
  NonlineLvl = length(unique(online))
  NliveLvl = length(unique(live))
  NtextLvl = length(unique(text))
  NLVLvl = length(unique(LV))
  NtransparencyLvl = length(unique(transparency))
  
 #data for bias part of the model
  
  bias = as.numeric(datFrmPredictor[,biasName])
  pollster = as.numeric(as.factor(datFrmPredictor[,pollsterName]))
  PollsterLevels = levels(as.factor(datFrmPredictor[,pollsterName]))
  year = as.numeric(as.factor(datFrmPredictor[,yearName]))
  YearLevels = levels(as.factor(datFrmPredictor[,yearName]))
  PollsTotal = length(bias)
  PollsterLevelsTotal = length(unique(pollster))
  YearLevelsTotal = length(unique(year))
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # For prior on baseline, etc.:
  biasMean = mean(bias)
  biasSD = sd(bias)
  
  
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # lmInfo = lm( datFrm[,scoreName] ~ datFrm[,samplesizeName] + 
  #                datFrm[,delModeName] + datFrm[,delModeName] + datFrm[,LVName] + datFrm[,transparencyName])
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
    IVR=IVR,
    online=online,
    live=live,
    text=text,
    LV = LV,
    transparency = transparency,
    samplesize = samplesize,
    
    # NIVRLvl=NIVRLvl,
    # NonlineLvl=NonlineLvl,
    # NliveLvl=NliveLvl,
    # NtextLvl=NtextLvl,
    
    NLVLvl = NLVLvl,
    NtransparencyLvl = NtransparencyLvl ,
    NraceIDLvl=NraceIDLvl,
    # data properties for scaling the prior:
    samplesizeSD = sd(samplesize) ,
    actualSD = sd(actual) ,
    agammaShRa = agammaShRa,
    agammaShRasamplesizeImpact = agammaShRasamplesizeImpact,
    
    
    #data for the bias part of the model
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
  
  print(online)
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  model {
  
  for ( race1 in 1:NraceIDLvl ){
    actual[race1] ~ dnorm(mu[race1], 1/actualSpread^2)
    
   for(myPoll in (whichrace[race1]+1):whichrace[race1+1]){
     weight[myPoll]=IVRImpact*(IVR[myPoll]-1) + onlineImpact*(online[myPoll]-1) +
     liveImpact*(live[myPoll]-1) + textImpact*(text[myPoll]-1)+
     LVImpact*(LV[myPoll]-1)+
     transparencyImpact*(transparency[myPoll]-1)
     +samplesizeImpact*samplesize[myPoll]
     
     
      bias[myPoll] ~ dt(mubias[myPoll], 1/biasSpread^2, nuY )
      mubias[myPoll] <- yearLean[year[myPoll]] + pollsterBias[pollster[myPoll],year[myPoll]] 
   
   }

    poll[(whichrace[race1]+1):whichrace[race1+1]] <- poll[(whichrace[race1]+1):whichrace[race1+1]] + bias[(whichrace[race1]+1):whichrace[race1+1]] 
    mu[race1] <- sum(nWeight[(whichrace[race1]+1):whichrace[race1+1]]*poll[(whichrace[race1]+1):whichrace[race1+1]])
    summedWeights[race1] = sum(weight[(whichrace[race1]+1):whichrace[race1+1]])
    
    for(myPoll1 in (whichrace[race1]+1):whichrace[race1+1]){
      nWeight[myPoll1]=weight[myPoll1]/summedWeights[race1]
    }
  }

IVRImpact ~ dgamma( agammaShRa[1] , agammaShRa[2] ) 
onlineImpact ~ dgamma( agammaShRa[1] , agammaShRa[2] ) 
liveImpact ~ dgamma( agammaShRa[1] , agammaShRa[2] ) 
textImpact ~ dgamma( agammaShRa[1] , agammaShRa[2] ) 

LVImpact ~ dgamma( agammaShRa[1] , agammaShRa[2] )  
   
    
transparencyImpact ~ dgamma( agammaShRa[1] , agammaShRa[2] )  
   
    
    samplesizeImpact ~ dgamma( agammaShRasamplesizeImpact[1] , agammaShRasamplesizeImpact[2] ) 
    actualSpread ~ dunif(actualSD/100, actualSD*10)
 

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
  
  parameters = c( "IVRImpact" , "onlineImpact", "liveImpact", "textImpact",
                  "samplesizeImpact" , "LVImpact" , "transparencyImpact", "actualSpread", "mu" 
  )
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

smryMCMC = function(  codaSamples , datFrm=NULL ,  LVName="LV" , 
                      IVRName="IVR", onlineName="online", liveName="live", textName="text",
                      transparencyName="transparency", 
                      samplesizeName ="samplesize",
                      #contrasts=NULL ,
                      saveName=NULL ) {
  # All single parameters:
  parameterNames = varnames(codaSamples) 
  if ( !is.null(datFrm) & !is.null(IVRName) ) {
    IVRlevels = levels(as.factor(datFrm[,IVRName]))
  }
  if ( !is.null(datFrm) & !is.null(onlineName) ) {
    onlinelevels = levels(as.factor(datFrm[,onlineName]))
  }
  if ( !is.null(datFrm) & !is.null(liveName) ) {
    livelevels = levels(as.factor(datFrm[,liveName]))
  }
  if ( !is.null(datFrm) & !is.null(textName) ) {
    textlevels = levels(as.factor(datFrm[,textName]))
  }
  
  
  
  if ( !is.null(datFrm) & !is.null(transparencyName) ) {
    transparencylevels = levels(as.factor(datFrm[,transparencyName]))
  }
  if ( !is.null(datFrm) & !is.null(LVName) ) {
    LVlevels = levels(as.factor(datFrm[,LVName]))
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
        # thisRowName = paste(thisRowName,delModelevels[levelVal]) 
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

plotLVPosterior = function( codaSamples , 
                            datFrmPredictor  , 
                            saveName=NULL , saveType="jpg",
                            showCurve = FALSE) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  LVName ="LV"
  pollName ="poll"
  
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  LV = as.numeric(as.factor(datFrmPredictor[,LVName]))
  LVLevels = levels(as.factor(datFrmPredictor[,LVName]))
  NLVLvl = length(unique(LV))
  
  openGraph(width=8,height=8)
  
  # posterior of the mean for sample size distrubtion 
  
  #give better titles - using LVLevels[1]=FALSE, LVLevels[2]=TRUE
  title="Likely-Voter Turnout Model Impact"
  
  plotPost( mcmcMat[,paste("LVImpact",sep="")], cex.lab = 1.75 , showCurve=showCurve ,
            xlab=bquote(LVImpact) , main=title)
  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,title, type=saveType))
    
  }
}


plotIVRPosterior = function( codaSamples , 
                             datFrmPredictor  , 
                             saveName=NULL , saveType="jpg",
                             showCurve = FALSE) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  IVRName ="IVR"
  pollName ="poll"
  
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  IVR = as.numeric(as.factor(datFrmPredictor[,IVRName]))
  IVRLevels = levels(as.factor(datFrmPredictor[,IVRName]))
  NIVRLvl = length(unique(IVR))
  
  openGraph(width=4,height=4)
  
  # posterior of the mean for sample size distrubtion 
  
  #give better titles - using delModLevels
  titleIVR="IVR Impact"
  
  
  plotPost( mcmcMat[,paste("IVRImpact",sep="")], cex.lab = 1.75 , showCurve=showCurve ,
            xlab=bquote(IVRImpact) , main=titleIVR)
  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,titleIVR, type=saveType))
    
  }
}


plotonlinePosterior = function( codaSamples , 
                                datFrmPredictor  , 
                                saveName=NULL , saveType="jpg",
                                showCurve = FALSE) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  onlineName ="online"
  pollName ="poll"
  
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  online = as.numeric(as.factor(datFrmPredictor[,onlineName]))
  onlineLevels = levels(as.factor(datFrmPredictor[,onlineName]))
  NonlineLvl = length(unique(online))
  
  openGraph(width=4,height=4)
  
  # posterior of the mean for sample size distrubtion 
  
  #give better titles - using delModLevels
  titleOnline="Online Impact"
  
  
  
  
  
  plotPost( mcmcMat[,paste("onlineImpact",sep="")], cex.lab = 1.75 , showCurve=showCurve ,
            xlab=bquote(onlineImpact) , main=titleOnline)
  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,titleOnline, type=saveType))
    
  }
}

plottextPosterior = function( codaSamples , 
                              datFrmPredictor  , 
                              saveName=NULL , saveType="jpg",
                              showCurve = FALSE) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  textName ="text"
  pollName ="poll"
  
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  text = as.numeric(as.factor(datFrmPredictor[,textName]))
  textLevels = levels(as.factor(datFrmPredictor[,textName]))
  NtextLvl = length(unique(text))
  
  openGraph(width=4,height=4)
  
  # posterior of the mean for sample size distrubtion 
  
  #give better titles - using delModLevels
  
  titleText="Text Impact"
  
  
  
  plotPost( mcmcMat[,paste("textImpact",sep="")], cex.lab = 1.75 , showCurve=showCurve ,
            xlab=bquote(textImpact) , main=titleText)
  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,titleText, type=saveType))
    
  }
}

plotlivePosterior = function( codaSamples , 
                              datFrmPredictor  , 
                              saveName=NULL , saveType="jpg",
                              showCurve = FALSE) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  liveName ="live"
  pollName ="poll"
  
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  live = as.numeric(as.factor(datFrmPredictor[,liveName]))
  liveLevels = levels(as.factor(datFrmPredictor[,liveName]))
  NliveLvl = length(unique(live))
  
  openGraph(width=4,height=4)
  
  # posterior of the mean for sample size distrubtion 
  
  
  
  titleLive="Live Impact"
  
  
  
  plotPost( mcmcMat[,paste("liveImpact",sep="")], cex.lab = 1.75 , showCurve=showCurve ,
            xlab=bquote(liveImpact) , main=titleLive)
  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,titleLive, type=saveType))
    
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
  
  openGraph(width=4,height=4)
  
  # posterior of the mean for sample size distrubtion 
  
  #give better titles - using delModLevels
  
  titleT="Transparent Impact"
  
  
  
  
  plotPost( mcmcMat[,paste("transparencyImpact",sep="")], cex.lab = 1.75 , showCurve=showCurve ,
            xlab=bquote(transparencyImpact) , main=titleT)
  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,titleT, type=saveType))
    
  }
}




