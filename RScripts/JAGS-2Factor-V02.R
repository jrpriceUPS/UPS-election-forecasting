#JAGS-2Factor-V02
#06-22-2020
#Modification of JAGS-2Factor.R 

source("DBDA2E-utilities.R")

#===============================================================================
genMCMC = function( datFrm , yName="y" , x1Name="x1" , x2Name="x2",
                    numSavedSteps=50000 , thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to generic x,y variable names for model:
  y = as.numeric(datFrm[,yName])
  x1 = as.numeric(as.factor(datFrm[,x1Name]))
  x1levels = levels(as.factor(datFrm[,x1Name]))
  x2 = as.numeric(as.factor(datFrm[,x2Name]))
  x2levels = levels(as.factor(datFrm[,x2Name]))
  Ntotal = length(y)
  Nx1Lvl = length(unique(x1))
  Nx2Lvl = length(unique(x2))
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # For prior on baseline, etc.:
  yMean = mean(y)
  ySD = sd(y)
  
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    y = y ,
    x1 = x1 ,
    x2 = x2 ,
    Ntotal = Ntotal ,
    Nx1Lvl = Nx1Lvl ,
    Nx2Lvl = Nx2Lvl ,
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
      y[i] ~ dt(mu[i], 1/ySigma^2, nuY )
      mu[i] <- a2[x2[i]] + a1[x1[i],x2[i]] 
    }
    nuY ~  dexp(1/30.0) 
    ySigma ~ dunif( ySD/100 , ySD*10 )
  
  
    # Middle level of the hierarchy (year lean)
    for ( j2 in 1:Nx2Lvl ) { a2[j2] ~ dnorm( 0.0 , 1/a2Sigma^2) }
    a2Sigma ~ dunif( ySD/100 , ySD*10 )
    
    # Middle level of the hierarchy (pollster bias with year)
     for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) { 
     a1[j1,j2] ~ dnorm(0.0 , 1/a2Sigma^2) }
     
    
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
  parameters = c( "ySigma" , "nuY" , "a2Sigma", "a1", "a2" )
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

smryMCMC = function(  codaSamples , datFrm=NULL , x1Name=NULL , x2Name = NULL,
                      contrasts=NULL , saveName=NULL ) {
  # All single parameters:
  parameterNames = varnames(codaSamples) 
  if ( !is.null(datFrm) & !is.null(x1Name) & !is.null(x2Name) ) {
    x1levels = levels(as.factor(datFrm[,x1Name]))
    x2levels = levels(as.factor(datFrm[,x2Name]))
  }
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName] ) )
    thisRowName = parName
    if ( !is.null(datFrm) & !is.null(x1Name) & !is.null(x2Name) ) {
      # For row name, extract numeric digits from parameter name. E.g., if
      # parameter name is "b1b2[12,34]" then pull out b1b2, 12 and 34:
      strparts = unlist( strsplit( parName , "\\[|,|\\]"  ) )
      # if there are only the param name and a single index:
      if ( length(strparts)==2 ) { 
        # if param name refers to factor 1:
        if ( substr(strparts[1],nchar(strparts[1]),nchar(strparts[1]))=="1" ) { 
          thisRowName = paste( thisRowName , x1levels[as.numeric(strparts[2])] )
        }
        # if param name refers to factor 2:
        if ( substr(strparts[1],nchar(strparts[1]),nchar(strparts[1]))=="2" ) { 
          thisRowName = paste( thisRowName , x2levels[as.numeric(strparts[2])] )
        }
      }
      # if there are the param name and two indices:
      if ( length(strparts)==3 ) { 
        thisRowName = paste( thisRowName , x1levels[as.numeric(strparts[2])], 
                             x2levels[as.numeric(strparts[3])] )
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
                     datFrm , yName="y" , x1Name="x1" , x2Name="x2" ,
                     x1contrasts=NULL , 
                     x2contrasts=NULL , 
                     x1x2contrasts=NULL ,
                     saveName=NULL , saveType="jpg",
                     showCurve = FALSE) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  y = datFrm[,yName]
  x1 = as.numeric(as.factor(datFrm[,x1Name]))
  x1levels = levels(as.factor(datFrm[,x1Name]))
  x2 = as.numeric(as.factor(datFrm[,x2Name]))
  x2levels = levels(as.factor(datFrm[,x2Name]))
# #Display data with posterior predictive distributions
#   for ( x2idx in 1:length(x2levels) ) {
#     openGraph(width=2*length(x1levels),height=5)
#     par( mar=c(4,4,2,1) , mgp=c(3,1,0) )
#     plot(-10,-10,
#          xlim=c(0.2,length(x1levels)+0.1) ,
#          xlab=paste(x1Name,x2Name,sep="\n") ,
#          xaxt="n" , ylab=yName ,
#          ylim=c(min(y)-0.2*(max(y)-min(y)),max(y)+0.2*(max(y)-min(y))) ,
#          main="Data with Post. Pred.")
#     axis( 1 , at=1:length(x1levels) , tick=FALSE ,
#           lab=paste( x1levels , x2levels[x2idx] , sep="\n" ) )
#     for ( x1idx in 1:length(x1levels) ) {
#       xPlotVal = x1idx #+ (x2idx-1)*length(x1levels)
#       yVals = y[ x1==x1idx & x2==x2idx ]
#       points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) ,
#               yVals , pch=1 , cex=1.5 , col="red" )
#       chainSub = round(seq(1,chainLength,length=20))
#       for ( chnIdx in chainSub ) {
#         m = mcmcMat[chnIdx,paste("a1[",x1idx,",", x2idx,"]",sep="")]   # pollster bias
# +         mcmcMat[chnIdx,paste("a2[",x2idx,"]",sep="")]  # year lean
#         s = mcmcMat[chnIdx,"ySigma"] # spread
#         nu = mcmcMat[chnIdx,"nuY"]# normality
# 
# 
#         tlim = qt( c(0.025,0.975) , df=nu )
#         yl = m+tlim[1]*s
#         yh = m+tlim[2]*s
#         ycomb=seq(yl,yh,length=201)
#         #ynorm = dnorm(ycomb,mean=m,sd=s)
#         #ynorm = 0.67*ynorm/max(ynorm)
#         yt = dt( (ycomb-m)/s , df=nu )
#         yt = 0.67*yt/max(yt)
#         lines( xPlotVal-yt , ycomb , col="skyblue" )
# 
# 
# 
#       }
#     }
#     if ( !is.null(saveName) ) {
#       saveGraph( file=paste0(saveName,"PostPred-",x2levels[x2idx]), type=saveType)
#     }
#   }# end for x2idx

  #openGraph(width=8,height=8)
  #layout(matrix(1:4,nrow=2))
  #layout(matrix(4:4))
  #layout(matrix(1:ceiling(length(x2levels)/4)*4,nrow=4))
  for ( x2idx in 1:length(x2levels) ) {
    openGraph(width=8,height=8)
    
    # posterior of the mean for that pollster
    histInfo = plotPost( mcmcMat[,paste("a2[",x2idx,"]",sep="")] , cex.lab = 1.75 , showCurve=showCurve ,
                         #compVal=compValMu , ROPE=ropeMu ,
                         xlab=x2levels[x2idx] , main=paste("Mean") ,
                         col="skyblue" )
    if ( !is.null(saveName) ) {
            saveGraph( file=paste0(saveName,"YearLean-",x2levels[x2idx]), type=saveType)
          }
  }
  
  #openGraph(width=8,height=8)
  #layout(matrix(1:4,nrow=2))
  
  #layout(matrix(4:4))
  #layout(matrix(1:ceiling(length(x2levels)/4)*4,nrow=4))
  for ( x1idx in 1:length(x1levels) ) {
    openGraph(width=8,height=8)
    
    # posterior of the mean for that pollster
    histInfo = plotPost(mcmcMat[,paste("a1[",x1idx,",", x2idx,"]",sep="")] , cex.lab = 1.75 , showCurve=showCurve ,
                         #compVal=compValMu , ROPE=ropeMu ,
                         xlab=x1levels[x1idx] , main=paste("Mean") ,
                         col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PollsterBias-",x1levels[x1idx]), type=saveType)
    }
  }
  
}
