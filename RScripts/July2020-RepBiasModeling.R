# First script used to model repBias with demBias, undecideds using JAGS
# demBias modeled with t-distribution, center determined by pollster + year (with interactions)
# repBias modeled with ANCOVA response to demBias and undecideds
#
# Notes by Jake 8/10/20

# Let us predict rebBias with demBias. 



source("DBDA2E-utilities.R")

#===============================================================================
genMCMC = function( datFrm , dembiasName = "demBias" , pollsterName = "pollster" , yearName = "year",
                    repBiasName="repBias", undecidedName="undecided",
                    numSavedSteps=50000 , thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to generic x,y variable names for model:
  undecided = as.numeric(datFrm[,undecidedName])
  dembias = as.numeric(datFrm[,dembiasName])
  repBias = as.numeric(datFrm[,repBiasName])
  pollster = as.numeric(as.factor(datFrm[,pollsterName]))
  PollsterLevels = levels(as.factor(datFrm[,pollsterName]))
  year = as.numeric(as.factor(datFrm[,yearName]))
  YearLevels = levels(as.factor(datFrm[,yearName]))
  PollsTotal = length(dembias)
  PollsterLevelsTotal = length(unique(pollster))
  YearLevelsTotal = length(unique(year))
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # For prior on baseline, etc.:
  dembiasMean = mean(dembias)
  dembiasSD = sd(dembias)
  repMean = mean(repBias)
  repSD = sd(repBias)
  undecidedMean=mean(undecided)
  undecidedSD=sd(undecided)
  
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    dembias = dembias,
    repBias=repBias,
    pollster = pollster ,
    undecided=undecided,
    undecidedMean=undecidedMean,
    undecidedSD=undecidedSD,
    year = year ,
    PollsTotal = PollsTotal ,
    PollsterLevelsTotal = PollsterLevelsTotal ,
    YearLevelsTotal =  YearLevelsTotal ,
    # data properties for scaling the prior:
    repMean = repMean,
    repSD = repSD,
    dembiasMean = dembiasMean ,
    dembiasSD = dembiasSD
  )

  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  model {
  
    #Bottom Level (individual poll rep bias):
    for (poll in 1:PollsTotal){
    repBias[poll] ~ dnorm((repResponsetoDemBias*dembias[poll] + UndecidedResponse*undecided[poll]), 1/repBiasSpread^2)
    dembias[poll] ~ dt(mu[poll] , (1/dembiasSpread^2) , nuY )
    mu[poll] <- yearLean[year[poll]] + pollsterBias[pollster[poll],year[poll]] 
    }
    nuY ~  dexp(1/30.0) 
    repBiasSpread ~ dunif( repSD/100 , repSD*10 )
    dembiasSpread ~ dunif( dembiasSD/100 , dembiasSD*10 )
    
    
    repResponsetoDemBias  ~ dnorm( 0 , 1/(10)^2 )
    UndecidedResponse  ~ dnorm( 0 , 1/(10)^2 )

   
    # Middle level of the hierarchy (year lean)
    for ( year in 1:YearLevelsTotal ) { yearLean[year] ~ dnorm( 0.0 , 1/yearSpread^2) }
    yearSpread ~ dunif( dembiasSD/100 , dembiasSD*10 )
    
    # Middle level of the hierarchy (pollster dem bias with year)
     for ( pollster in 1:PollsterLevelsTotal ) { for ( year in 1:YearLevelsTotal ) { 
     pollsterBias[pollster,year] ~ dnorm(0.0 , 1/pollsterSpread^2) }
     
    
     }
     pollsterSpread ~ dunif( dembiasSD/100 , dembiasSD*10 )
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
  parameters = c(  "repResponsetoDemBias", "UndecidedResponse","repBiasSpread","dembiasSpread" ,
                  "nuY" , "pollsterSpread", "yearSpread", "yearLean", "pollsterBias" )
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

plotDiagnostics= function( ){
  for ( parName in c("dembiasSpread",  "nuY" , "pollsterSpread" , "yearSpread" , "yearLean[1]", "pollsterBias[1,1]"  ) ) {
    diagMCMC( codaObject=mcmcCodaV03 , parName=parName , 
              saveName=fileNameRoot , saveType=graphFileType )
  }
}
#===============================================================================

#===============================================================================

plotPosteriorPredictiveDem = function( codaSamples , 
                                    datFrm , dembiasName=NULL , 
                                    saveName=NULL , saveType="jpg",
                                    showCurve = FALSE) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  dembiasName = "demBias" 
  pollsterName = "pollster" 
  repBiasName="repBias"
  yearName = "year"
  undecidedName="undecided"
  demBias = datFrm[,dembiasName]
  repBias= datFrm[,repBiasName]
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
         xaxt="n" , ylab=dembiasName ,
         ylim=c(min(demBias)-0.2*(max(demBias)-min(demBias)),max(demBias)+0.2*(max(demBias)-min(demBias))) ,
         main="Data with Post. Pred.")
    axis( 1 , at=1:length(PollsterLevels) , tick=FALSE ,
          lab=paste( PollsterLevels , YearLevels[Yearidx] , sep="\n" ) )
    for ( Pollsteridx in 1:length(PollsterLevels) ) {
      xPlotVal = Pollsteridx #+ (Yearidx-1)*length(PollsterLevels)
      yVals = demBias[ pollster==Pollsteridx & year==Yearidx ]
      points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) ,
              yVals , pch=1 , cex=1.5 , col="red" )
      chainSub = round(seq(1,chainLength,length=20))
      for ( chnIdx in chainSub ) {
        m = mcmcMat[chnIdx,paste("pollsterBias[",Pollsteridx,",", Yearidx,"]",sep="")]   # pollster bias
        +         mcmcMat[chnIdx,paste("yearLean[",Yearidx,"]",sep="")]  # year lean
        s = mcmcMat[chnIdx,"dembiasSpread"] # spread
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
      saveGraph( file=paste0(saveName,"PostPred-Dem-",YearLevels[Yearidx]), type=saveType)
    }
  }# end for Yearidx


# 
#   for ( Yearidx in 1:length(YearLevels) ) {
#     openGraph(width=2*length(PollsterLevels),height=5)
#     par( mar=c(4,4,2,1) , mgp=c(3,1,0) )
#     plot(-10,-10,
#          xlim=c(0.2,length(PollsterLevels)+0.1) ,
#          xlab=paste(pollsterName,yearName,sep="\n") ,
#          xaxt="n" , ylab=repBiasName ,
#          ylim=c(min(repBias)-0.2*(max(repBias)-min(repBias)),max(repBias)+0.2*(max(repBias)-min(repBias))) ,
#          main="Data with Post. Pred.")
#     axis( 1 , at=1:length(PollsterLevels) , tick=FALSE ,
#           lab=paste( PollsterLevels , YearLevels[Yearidx] , sep="\n" ) )
#     for ( Pollsteridx in 1:length(PollsterLevels) ) {
#       xPlotVal = Pollsteridx #+ (Yearidx-1)*length(PollsterLevels)
#       yVals = repBias[ pollster==Pollsteridx & year==Yearidx ]
#       points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) ,
#               yVals , pch=1 , cex=1.5 , col="red" )
#      
#      # repBias[poll] ~ dnorm((repResponsetoDemBias*dembias[poll] + UndecidedResponse*undecided[poll]), 1/repBiasSpread^2)
#       
#        chainSub = round(seq(1,chainLength,length=20))
#       for ( chnIdx in chainSub ) {
#         m = mcmcMat[chnIdx,paste("repResponsetoDemBias")]*(mcmcMat[chnIdx,paste("demBias[",Pollsteridx,",", Yearidx,"]",sep="")]   # pollster bias
#         +         mcmcMat[chnIdx,paste("yearLean[",Yearidx,"]",sep="")])  # year lean
#         s = mcmcMat[chnIdx,"repBiasSpread"] # spread
# 
# 
# 
#         normlim = qnorm( c(0.025,0.975) )
#         yl = m+normlim[1]*s
#         yh = m+normlim[2]*s
#         ycomb=seq(yl,yh,length=201)
#         ynorm = dnorm(ycomb,mean=m,sd=s)
#         ynorm = 0.67*ynorm/max(ynorm)
#         
#         lines( xPlotVal-ynorm , ycomb , col="skyblue" )
# 
# 
# 
#       }
#     }
#     if ( !is.null(saveName) ) {
#       saveGraph( file=paste0(saveName,"PostPred-Rep-",YearLevels[Yearidx]), type=saveType)
#     }
# 
# }
}


#===============================================================================  

#===============================================================================

plotMCMCrep = function( codaSamples , data , xName="undecided" , yName="repBias" ,
                    
                     showCurve=FALSE ,  
                     saveName=NULL , saveType="jpg" ) {
  # Marginal histograms:
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  repResponsetoDemBias = mcmcMat[,"repResponsetoDemBias"]
  UndecidedResponse  = mcmcMat[,"UndecidedResponse"]
  repBiasSpread = mcmcMat[,"repBiasSpread"]
  decideOpenGraph = function( panelCount , saveName , finished=FALSE , 
                              nRow=2 , nCol=3 ) {
    
    # If finishing a set:
    if ( finished==TRUE ) {
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,ceiling((panelCount-1)/(nRow*nCol))), 
                   type=saveType)
      }
      panelCount = 1 # re-set panelCount
      return(panelCount)
    } else {
      # If this is first panel of a graph:
      if ( ( panelCount %% (nRow*nCol) ) == 1 ) {
        # If previous graph was open, save previous one:
        if ( panelCount>1 & !is.null(saveName) ) {
          saveGraph( file=paste0(saveName,(panelCount%/%(nRow*nCol))), 
                     type=saveType)
        }
        # Open new graph
        openGraph(width=nCol*9.0/3,height=nRow*2.0)
        layout( matrix( 1:(nRow*nCol) , nrow=nRow, byrow=TRUE ) )
        par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
      }
      # Increment and return panel count:
      panelCount = panelCount+1
      return(panelCount)
    }
  }
  # Original scale:
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( repResponsetoDemBias , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(repResponsetoDemBias) , main="Republican Response to Democratic Bias" )

    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
    histInfo = plotPost( UndecidedResponse , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(UndecidedResponse) , main="Rebulican Response to Undecided" )
  
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( repBiasSpread , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(repBiasSpread) , main=paste("Scale") )
  
 
  panelCount = decideOpenGraph( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMarg") )
  }
