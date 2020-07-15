# Let us predict rebBias with demBias. 



source("DBDA2E-utilities.R")

#===============================================================================
genMCMC = function( datFrm , biasName = "demBias" , pollsterName = "pollster" , yearName = "year",
                    repBiasName="repBias", undecidedName="undecided",
                    numSavedSteps=50000 , thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to generic x,y variable names for model:
  undecided = as.numeric(datFrm[,undecidedName])
  bias = as.numeric(datFrm[,biasName])
  repBias = as.numeric(datFrm[,repBiasName])
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
  repMean = mean(repBias)
  repSD = sd(repBias)
  undecidedMean=mean(undecided)
  undecidedSD=sd(undecided)
  
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    bias = bias,
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
    biasMean = biasMean ,
    biasSD = biasSD
  )

  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  model {
  
    #Bottom Level (individual poll rep bias):
    for (poll in 1:PollsTotal){
    repBias[poll] ~ dnorm((beta0+beta1*bias[poll] + beta2*undecided[poll]), 1/repBiasSpread^2)
    bias[poll] ~ dt(mu[poll] , (1/biasSpread^2) , nuY )
    mu[poll] <- yearLean[year[poll]] + pollsterBias[pollster[poll],year[poll]] 
    }
    nuY ~  dexp(1/30.0) 
    repBiasSpread ~ dunif( repSD/100 , repSD*10 )
    biasSpread ~ dunif( biasSD/100 , biasSD*10 )
    
    beta0 ~ dnorm( 0 , 1/(10)^2 )
    beta1  ~ dnorm( 0 , 1/(10)^2 )
    beta2  ~ dnorm( 0 , 1/(10)^2 )

   
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
  parameters = c( "beta0", "beta1", "beta2","repBiasSpread","biasSpread" ,
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

smryMCMC = function(  codaSamples , datFrm=NULL , biasName = "demBias" , pollsterName = "pollster" , yearName = "year",
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
                       "beta0" = summarizePost( mcmcMat[,"beta0"] 
                                                #, compVal=compValBeta0 , 
                                                #ROPE=ropeBeta0 
                                                ) )
  summaryInfo = rbind( summaryInfo , 
                       "beta1" = summarizePost( mcmcMat[,"beta1"]  
                                                
                                              #,compVal=compValBeta1 , 
                                                #ROPE=ropeBeta1 
                                              ) )
 
  summaryInfo = rbind( summaryInfo , 
                       "beta2" = summarizePost( mcmcMat[,"beta2"] 
                                               # ,compVal=compValBeta1 , 
                                                #ROPE=ropeBeta1 
                                               ) )
  
  # Save results:
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotDiagnostics= function( ){
  for ( parName in c("biasSpread",  "nuY" , "pollsterSpread" , "yearSpread" , "yearLean[1]", "pollsterBias[1,1]"  ) ) {
    diagMCMC( codaObject=mcmcCodaV03 , parName=parName , 
              saveName=fileNameRoot , saveType=graphFileType )
  }
}
#===============================================================================
