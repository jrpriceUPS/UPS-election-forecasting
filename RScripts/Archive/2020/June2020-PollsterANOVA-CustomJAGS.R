# modification of code from Kruschke's DBDA2E book
# does 2-factor ANOVA to predict bias using pollster and year
# includes more customized model definition and plotting
# called by "June2020-PollsterANOVA-CustomJAGS-Example"
#
# Runs!
# Notes by Jake, 8/6/20

source("DBDA2E-utilities.R")

#===============================================================================
genMCMC = function( datFrm , yName="y" , xName="x" ,
                    numSavedSteps=50000 , thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to generic x,y variable names for model:
  y = as.numeric(datFrm[,yName]) # extract metric predicted as a numeric vector
  x = as.numeric(as.factor(datFrm[,xName])) # extract nominal predictor as numeric vector
  xlevels = levels(as.factor(datFrm[,xName])) # make a reference sheet of nominal levels (which pollster corresponds to 1, 2, etc.)
  Ntotal = length(y)
  NxLvl = length(unique(x))
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # For prior on baseline, etc.:
  yMean = mean(y)
  ySD = sd(y)
  
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    y = y ,
    x = x ,
    Ntotal = Ntotal ,
    NxLvl = NxLvl ,
    # data properties for scaling the prior:
    yMean = yMean ,
    ySD = ySD
  )
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  model {
  
    # Bottom level of the hierarchy (individual poll bias)
    for ( i in 1:Ntotal ) {
      y[i] ~ dt( a[x[i]], 1/ySigma^2, nuY )
    }
    nuY ~  dexp(1/30.0) 
    ySigma ~ dunif( ySD/100 , ySD*10 )
    
    
    # Middle level of the hierarchy (pollster bias)
    for ( j in 1:NxLvl ) { a[j] ~ dnorm( a0 , 1/aSigma^2) }
    a0 ~ dnorm( 0 , 1/(ySD*5)^2 )
    aSigma ~ dunif( ySD/100 , ySD*10 )
    
  }
  " # close quote for modelstring
  writeLines(modelstring,con="TEMPmodel.txt")
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  #initsList = list(
  #  a0 = yMean ,
  #  a = aggregate( y , list( x ) , mean )[,2] - yMean ,
  #  ySigma = mean( aggregate( y , list( x ) , sd )[,2] )
    # Let JAGS do other parameters automatically...
  #)
  #------------------------------------------------------------------------------
  # RUN THE CHAINS
  
  require(runjags)
  parameters = c( "a" ,  "ySigma" , "nuY" , "a0", "aSigma" )
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

smryMCMC = function(  codaSamples , datFrm=NULL , xName=NULL ,
                      contrasts=NULL , saveName=NULL ) {
  # All single parameters:
  parameterNames = varnames(codaSamples) 
  if ( !is.null(datFrm) & !is.null(xName) ) {
    xlevels = levels(as.factor(datFrm[,xName]))
  }
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName] ) )
    thisRowName = parName
    if ( !is.null(datFrm) & !is.null(xName) ) {
      # For row name, extract numeric digits from parameter name. E.g., if
      # parameter name is "beta[12,34]" then pull out 12 and 34:
      levelVal = as.numeric( 
        grep( "^[1-9]" , # grep only substrings that begin with digits.
              # Return sll substrings split by "[" or "," or "]":
              unlist( strsplit( parName , "\\[|,|\\]"  ) ) , 
              value=TRUE ) )
      if ( length(levelVal) > 0 ) { 
        # Assumes there is only a single factor, i.e., levelVal has only entry: 
        thisRowName = paste(thisRowName,xlevels[levelVal]) 
      }
    }
    rownames(summaryInfo)[NROW(summaryInfo)] = thisRowName
  }
  # All contrasts:
  # if ( !is.null(contrasts) ) {
  #   if ( is.null(datFrm) | is.null(xName) ) {
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

plotMCMC = function( codaSamples , 
                     datFrm , yName="y" , xName="x" , contrasts=NULL ,
                     saveName=NULL , saveType="png", showCurve = FALSE  ) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  y = datFrm[,yName]
  x = as.numeric(as.factor(datFrm[,xName]))
  xlevels = levels(as.factor(datFrm[,xName]))
  # Display data with posterior predictive distributions
  # openGraph(width=min(14,4.5*length(xlevels)),height=6)
  # par(mar=c(3,3,2,0.5)) # number of margin lines: bottom,left,top,right
  # par(mgp=c(1.75,0.5,0)) # which margin lines to use for labels
  # plot(-1,0, 
  #      xlim=c(0.1,length(xlevels)+0.1) , 
  #      xlab=xName , cex=.25, xaxt="n" , ylab=yName ,
  #      ylim=c(min(y)-0.2*(max(y)-min(y)),max(y)+0.2*(max(y)-min(y))) , 
  #      main="Data with Posterior Predictive Distrib.")
  # axis( 1 , at=1:length(xlevels) , tick=FALSE , lab=xlevels )
  # for ( xidx in 1:length(xlevels) ) {
  #   xPlotVal = xidx  
  #   yVals = y[ x==xidx ]
  #   points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) , 
  #           yVals , pch=1 , cex=1.5 , col="red" )
  #   chainSub = round(seq(1,chainLength,length=20))
  #   for ( chnIdx in chainSub ) {
  #     m = mcmcMat[chnIdx,paste("a[",xidx,"]",sep="")] # grab the center for the current pollster at that step in chain
  #     s = mcmcMat[chnIdx,paste("ySigma",sep="")] # grab the standard deviation for the current pollster at that step in chain
  #     nu = mcmcMat[chnIdx,"nuY"] # grab the normality parameter at that step
  #     tlim = qt( c(0.025,0.975) , df=nu )
  #     yl = m+tlim[1]*s
  #     yh = m+tlim[2]*s
  #     ycomb=seq(yl,yh,length=201)
  #     #ynorm = dnorm(ycomb,mean=m,sd=s)
  #     #ynorm = 0.67*ynorm/max(ynorm)
  #     yt = dt( (ycomb-m)/s , df=nu )
  #     yt = 0.67*yt/max(yt)
  #     lines( xPlotVal-yt , ycomb , col="skyblue" ) 
  #   }
  # }
   if ( !is.null(saveName) ) {
     saveGraph( file=paste(saveName,"PostPred",sep=""), type=saveType)
  }
  
  # plot posteriors for each pollster
  for ( xidx in 1:length(xlevels) ) {
    openGraph(width=3,height=8)
    layout(matrix(1:2,ncol=1))
    
    
    # Compute limits for plots of data with posterior pred. distributions
    y = datFrm[,yName]
    y = y[x==xidx]
    xLim = c( min(y)-0.1*(max(y)-min(y)) , max(y)+0.1*(max(y)-min(y)) )
    xBreaks = seq( xLim[1] , xLim[2] , 
                   length=ceiling((xLim[2]-xLim[1])/(sd(y)/4)) )
    histInfo = hist(y,breaks=xBreaks,plot=FALSE)
    yMax = 1.2 * max( histInfo$density )
    xVec = seq( xLim[1] , xLim[2] , length=501 )
    #-----------------------------------------------------------------------------
    # Plot data y and smattering of posterior predictive curves:
    histInfo = hist( y , prob=TRUE , xlim=xLim , ylim=c(0,yMax) , breaks=xBreaks,
                     col="red2" , border="white" , xlab="y" , ylab="" , 
                     yaxt="n" , cex.lab=1.5 , main="Data w. Post. Pred." )
    chainSub = round(seq(1,chainLength,length=20))
    for ( chnIdx in chainSub ) {
      m = mcmcMat[chnIdx,paste("a[",xidx,"]",sep="")] # grab the center for the current pollster at that step in chain
      s = mcmcMat[chnIdx,paste("ySigma",sep="")] # grab the standard deviation for the current pollster at that step in chain
      nu = mcmcMat[chnIdx,"nuY"] # grab the normality parameter at that step
      
      lines(xVec, dt( (xVec-m)/s, 
                      df=nu)/s, 
            type="l" , col="skyblue" , lwd=1 )
    }
    text( max(xVec) , yMax , bquote(N==.(length(y))) , adj=c(1.1,1.1) )
    
    # posterior of the mean for that pollster
    histInfo = plotPost( mcmcMat[,paste("a[",xidx,"]",sep="")] , cex.lab = 1.75 , showCurve=showCurve ,
                         #compVal=compValMu , ROPE=ropeMu ,
                         xlab=xlevels[xidx] , main=paste("Mean") , 
                         col="skyblue" )
    
    

    
  }
  
  openGraph(width=10,height=3)
  layout(matrix(1:2,ncol=2))
  histInfo = plotPost( mcmcMat[,"ySigma"] , cex.lab = 1.75 , showCurve=showCurve ,
                       #compVal=compValMu , ROPE=ropeMu ,
                       xlab="Pollster Standard Deviation" , main=paste("Standard Deviation") , 
                       col="skyblue" )
  
  histInfo = plotPost( mcmcMat[,"nuY"] , cex.lab = 1.75 , showCurve=showCurve ,
                       #compVal=compValMu , ROPE=ropeMu ,
                       xlab="Pollster Normality" , main=paste("Normality") , 
                       col="skyblue" )
  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Normality&SD",sep=""), type=saveType)
  }
  
  
  openGraph(width=10,height=3)
  layout(matrix(1:2,ncol=2))
  histInfo = plotPost( mcmcMat[,"a0"] , cex.lab = 1.75 , showCurve=showCurve ,
                       #compVal=compValMu , ROPE=ropeMu ,
                       xlab="Average Pollster Bias" , main=paste("Average Bias") , 
                       col="skyblue" )
  
  histInfo = plotPost( mcmcMat[,"aSigma"] , cex.lab = 1.75 , showCurve=showCurve ,
                       #compVal=compValMu , ROPE=ropeMu ,
                       xlab="Spread of Pollster Bias" , main=paste("Spread") , 
                       col="skyblue" )

  
  # if ( !is.null(contrasts) ) {
  #   if ( is.null(datFrm) | is.null(xName) ) {
  #     show(" *** YOU MUST SPECIFY THE DATA FILE AND FACTOR NAMES TO DO CONTRASTS. ***\n")
  #   } else {
  #     for ( cIdx in 1:length(contrasts) ) {
  #       thisContrast = contrasts[[cIdx]]
  #       left = right = rep(FALSE,length(xlevels))
  #       for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
  #         left = left | xlevels==thisContrast[[1]][nIdx]
  #       }
  #       left = normalize(left)
  #       for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
  #         right = right | xlevels==thisContrast[[2]][nIdx]
  #       }
  #       right = normalize(right)
  #       contrastCoef = matrix( left-right , ncol=1 )
  #       postContrast = ( mcmcMat[,paste("b[",1:length(xlevels),"]",sep="")] 
  #                        %*% contrastCoef )
  #       openGraph(height=8,width=4)
  #       layout(matrix(1:2,ncol=1))
  #       plotPost( postContrast , xlab="Difference" ,
  #                 main=paste0( 
  #                   paste(thisContrast[[1]],collapse="."), 
  #                   "\nvs\n",
  #                   paste(thisContrast[[2]],collapse=".") ) ,
  #                 compVal=thisContrast$compVal , ROPE=thisContrast$ROPE )
  #       plotPost( postContrast/mcmcMat[,"ySigma"] , xlab="Effect Size" ,
  #                 main=paste0( 
  #                   paste(thisContrast[[1]],collapse="."), 
  #                   "\nvs\n",
  #                   paste(thisContrast[[2]],collapse=".") ) ,
  #                 compVal=0.0 , 
  #                 ROPE=c(-0.1,0.1) )
  #     if ( !is.null(saveName) ) {
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

