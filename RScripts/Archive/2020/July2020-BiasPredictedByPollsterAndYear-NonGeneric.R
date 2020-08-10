# first attempt at including dem error, rep. error, and margin error all in one posterior
# This runs only a single example of the above and also includes some plots - similar to PollsterYearPredBias-Interactions
# non-generic naming scheme

# Runs! Sources "RScripts/Jags-2Factor-V03.R" at the end so I'll need to come back and rename that eventually
# Notes by Jake, 8/6/20


#JAGS-2Factor-V04
#3PartPosterior

#06-29-2020
#Modification of JAGS-2Factor-V03.R 
#Change so that it produces a posterior graph with 3Posterior distribution lines drawn


source("DBDA2E-utilities.R")

#===============================================================================
genMCMC = function( datFrm , biasName = "bias" , 
                    pollsterName = "pollster" , yearName = "year",
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

##############################################

plotPosteriorPredictive = function( codaSamplesbias , codaSamplesrepBias, codaSamplesdemBias,
                                    datFrm , biasName=NULL , 
                                    saveName=NULL , saveType="jpg",
                                    showCurve = FALSE, pollsterName="pollster" ) {
  mcmcMat = as.matrix(codaSamplesbias,chains=TRUE)
  chainLength = NROW( mcmcMat )
  #pollsterName="pollster"
  yearName="year"
  demBias = datFrm[,"repBias"]
  repBias = datFrm[,"demBias"]
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
      
      yVals = demBias[ pollster==Pollsteridx & year==Yearidx ]
      points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) ,
              yVals , pch=1 , cex=1.5 , col="blue" )
      
      yVals = repBias[ pollster==Pollsteridx & year==Yearidx ]
      points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) ,
              yVals , pch=1 , cex=1.5 , col="red" )
      
      yVals = bias[ pollster==Pollsteridx & year==Yearidx ]
      points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) ,
              yVals , pch=1 , cex=1.5 , col="black" )
      
      
      
      chainSub = round(seq(1,chainLength,length=10))
      
      
      #Marginal Bias
      mcmcMat = as.matrix(codaSamplesbias,chains=TRUE)
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
        lines( xPlotVal-yt , ycomb , col="black" )
      }
      
      
    
      
      
      
      #Dem Bias
      mcmcMatD = as.matrix(codaSamplesdemBias,chains=TRUE)
      for ( chnIdx in chainSub ) {
        m = mcmcMatD[chnIdx,paste("pollsterBias[",Pollsteridx,",", Yearidx,"]",sep="")]   # pollster bias
        +         mcmcMatD[chnIdx,paste("yearLean[",Yearidx,"]",sep="")]  # year lean
        s = mcmcMatD[chnIdx,"biasSpread"] # spread
        nu = mcmcMatD[chnIdx,"nuY"]# normality
        
        
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
      
      
      #Rep Bias
      mcmcMatR = as.matrix(codaSamplesrepBias,chains=TRUE)
      for ( chnIdx in chainSub ) {
        m = mcmcMatR[chnIdx,paste("pollsterBias[",Pollsteridx,",", Yearidx,"]",sep="")]   # pollster bias
        +         mcmcMat[chnIdx,paste("yearLean[",Yearidx,"]",sep="")]  # year lean
        s = mcmcMatR[chnIdx,"biasSpread"] # spread
        nu = mcmcMatR[chnIdx,"nuY"]# normality
        
        
        tlim = qt( c(0.025,0.975) , df=nu )
        yl = m+tlim[1]*s
        yh = m+tlim[2]*s
        ycomb=seq(yl,yh,length=201)
        #ynorm = dnorm(ycomb,mean=m,sd=s)
        #ynorm = 0.67*ynorm/max(ynorm)
        yt = dt( (ycomb-m)/s , df=nu )
        yt = 0.67*yt/max(yt)
        lines( xPlotVal-yt , ycomb , col="red" )
      }
      
     
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PostPred-",YearLevels[Yearidx]), type=saveType)
    }
  }# end for Yearidx
  
}





#This version is for creating visuals for the presentation.
plotPosteriorPredictiveV02 = function( codaSamplesbias , codaSamplesrepBias, codaSamplesdemBias,
                                    datFrm , biasName=NULL , 
                                    saveName=NULL , saveType="jpg",
                                    showCurve = FALSE, pollsterName="pollster" ) {
  mcmcMat = as.matrix(codaSamplesbias,chains=TRUE)
  chainLength = NROW( mcmcMat )
  #pollsterName="pollster"
  yearName="year"
  demBias = datFrm[,"demBias"]
  repBias = datFrm[,"repBias"]
  bias = datFrm[,biasName]
  pollster = as.numeric(as.factor(datFrm[,pollsterName]))
  PollsterLevels = levels(as.factor(datFrm[,pollsterName]))
  year = as.numeric(as.factor(datFrm[,yearName]))
  YearLevels = levels(as.factor(datFrm[,yearName]))
  # Display data with posterior predictive distributions
  for ( Yearidx in 1:length(YearLevels) ) {
    openGraph(width=5*length(PollsterLevels),height=7)
    par( mar=c(4,4,2,1) , mgp=c(3,1,0) )
    plot(-10,-10,
         xlim=c(.2,length(PollsterLevels)+0.1) ,
         xlab="Pollster Name & Year" ,
         xaxt="n" , ylab="Bias" ,
         ylim=c(min(bias)-0.2*(max(bias)-min(bias)),max(bias)+0.2*(max(bias)-min(bias))) ,
         main=paste("Posterior Predictions for Margin of Bias - ", YearLevels[Yearidx])
         
         )
    axis( 1 , at=1:length(PollsterLevels) , tick=FALSE ,
          lab=paste( PollsterLevels , YearLevels[Yearidx] , sep="\n" ) )
    for ( Pollsteridx in 1:length(PollsterLevels) ) {
      xPlotVal = Pollsteridx #+ (Yearidx-1)*length(PollsterLevels)
      

      
      yVals = bias[ pollster==Pollsteridx & year==Yearidx ]
      points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) ,
              yVals , pch=1 , cex=1.5 , col="black" )
      
      
      
      chainSub = round(seq(1,chainLength,length=20))
      
      
      #Marginal Bias
      mcmcMat = as.matrix(codaSamplesbias,chains=TRUE)
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
        lines( xPlotVal-yt , ycomb , col="black" )
      }
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PostPred-JustMargin",YearLevels[Yearidx]), type=saveType)
    }
  }
      
  
  
  #open new graph
  for ( Yearidx in 1:length(YearLevels) ) {
    openGraph(width=5*length(PollsterLevels),height=7)
    par( mar=c(4,4,2,1) , mgp=c(3,1,0) )
    plot(-10,-10,
         xlim=c(.2,length(PollsterLevels)+0.1) ,
         xlab="Pollster Name & Year" ,
         xaxt="n" , ylab="Bias" ,
         ylim=c(min(bias)-0.2*(max(bias)-min(bias)),max(bias)+0.2*(max(bias)-min(bias))) ,
         main=paste("Comparing Posterior Predictions for Different Measures of Bias - ", YearLevels[Yearidx])
         
    )
    axis( 1 , at=1:length(PollsterLevels) , tick=FALSE ,
          lab=paste( PollsterLevels , YearLevels[Yearidx] , sep="\n" ) )
    for ( Pollsteridx in 1:length(PollsterLevels) ) {
      xPlotVal = Pollsteridx #+ (Yearidx-1)*length(PollsterLevels)
      
      yVals = demBias[ pollster==Pollsteridx & year==Yearidx ]
      points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) ,
              yVals , pch=1 , cex=1.5 , col="blue" )
      
      yVals = repBias[ pollster==Pollsteridx & year==Yearidx ]
      points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) ,
              yVals , pch=1 , cex=1.5 , col="red" )
      
      chainSub = round(seq(1,chainLength,length=10))
      
      #Dem Bias
      mcmcMat = as.matrix(codaSamplesdemBias,chains=TRUE)
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
      
      
      #Rep Bias
      mcmcMat = as.matrix(codaSamplesrepBias,chains=TRUE)
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
        lines( xPlotVal-yt , ycomb , col="red" )
      }
      
      legend( x = c(.12, .67), y = c(-12, -17), legend=c("Democratic Error", "Republican Error"), col=c("skyblue" ,"red"), lty=1:2, cex=0.9)
      abline(a=0,b=0, col="ivory4")
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PostPred-",YearLevels[Yearidx]), type=saveType)
    }
  }# end for Yearidx
  
}

