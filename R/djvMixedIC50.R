# http://kbroman.org/pkg_primer/pages/docs.html
# http://cran.r-project.org/web/packages/roxygen2/index.html


#' Estimating IC50 values using a mixed effects model.
#'
#' \pkg{djvMixedIC50} allows you to estimate IC50 and AUC values using a non-linear mixed model. It estimates all the responses simultaneously and capitalizes on the entire set of responses for each cell line to infer its shape. This significantly reduces the ambiguity of the cell line's shape and thus improves the estimate.
#'
#' The package was initially designed for the GDSC compound sensitivity screen produced by the Wellcome Trust Sanger Institute. Their typical assay is a 9 point assay with 2-fold dilutions. However, all screen comprising sets of cell lines exposed to a series of compounds can benefit from this modeling approach.
#' 
#' As a design choice, the IC50 is not directly estimated on the concentration range, but on the dilution steps. In this way, the internal IC50 estimates are comparable in terms of response range. Internally the fold-dilution is constant at 2, other constants can be accommodated by using a custom predictor vector (see also \code{\link{getXfromConcSeries}}). 
#' 
#' After the initialization of the data by \code{\link{initIC50}}, and the curve fitting of the data by \code{\link{fitModel}}, the data frames containing the resultant statistics are constructed by \code{\link{gatherModelStats}}. In this function, the internal IC50 estimates are translated to actual (log) concentration values based on the highest test concentration. The individual data matrices holding these statistics can be accessed the named list. For more information on how to format the data for a custom screen, please consult the relevant sections in \code{\link{initIC50}}.
#' 
#' The resultant matrices hold the natural logarithm of the estimated micro Molar IC50 concentration. The area under the curve (AUC) is calculated using the integral of the estimated model parameters for the dose-response curve and normalized by the total area. The numerically integrated version is designated AUCtrap, calculated using the trapezoid rule. The residual of the model fit is summarized by the root mean square error.
#' 
#' See also: Pharmacogenomics, May 2016, Vol. 17, No. 7, Pages 691-700, or \url{http://www.ncbi.nlm.nih.gov/pubmed/27180993}
#'
#' 
#' @docType package
#' @name djvMixedIC50
#' 
#' @examples
#' gDat      <- initIC50(szFileName='~/myFile.csv')
#' fmMod1    <- fitModel(gDat)
#' outStats  <- gatherModelStats(gDat=gDat,fmMod1=fmMod1)
#' dfIC50    <- outStats$IC50
#' dfRMSE    <- outStats$RMSE
#' dfAUC     <- outStats$AUC
#' dfAUCtrap <- outStats$AUCtrap
#' 
#' The input data looks as follows:
#' 
#' Grouped Data: y ~ x | drug/CL
#'   x          y     CL drug maxc
#' 1 9 0.88079708 MC-CAR    1    2
#' 2 8 0.79139147 MC-CAR    1    2
#' 3 7 0.66075637 MC-CAR    1    2
#' 4 6 0.50000000 MC-CAR    1    2
#' 5 5 0.33924363 MC-CAR    1    2
#' 6 4 0.20860853 MC-CAR    1    2
#' 7 3 0.11920292 MC-CAR    1    2
#' 8 2 0.06496917 MC-CAR    1    2
#' 9 1 0.03444520 MC-CAR    1    2
NULL

#' Preparing data for curve fitting
#'
#' Prepares data for use in \code{\link{fitModel}}. As the model was initially developed for estimating IC50 values from the GDSC, this function serves as a convenience function to process this data into a form suitable for curve fitting. It takes the default GDSC file format and processes the dose-response and control well intensities to obtain the relative viabilities. These viabilities returned, together with the reverse-order dilution steps, the cell line name, the drug, and the maximum screening concentration. The object returned is a \code{\link{groupedData}}. 
#' 
#' Please see the section on Custom Screens for more information how to prepare the data for fitting.
#' 
#' @section Custom screen:
#' 
#' Internally, the fitting considers the viability to be a function of the dilution steps. In this way, the maximum screening concentration always takes the same 'x' value, and the variation in the parameters can be interpreted in terms of the dilution steps. By design, the highest screening concentration is the 9'th dilution step, and the lowest screening concentration is the first dilution step. The 'y' value is the 1-viability. The CL value is the cell line name. The drug value is the drug identifier, and the maxc hold the maximum concentration.
#' 
#' A data frame holding the data needs to contain the following (literal) column names: "x","y","CL","drug", and "maxc". In which x is an integer of value 9 and lower (see also \code{\link{getXfromConcSeries}}), y is 1-viability [0..1], CL is the cell line name (factor), drug is the drug identifier (factor), and finally maxc is the maximum test concentration (float, >0). 
#' 
#' The data frame is then converted to a groupedData object by the following call (dfDat is the data frame holding the data).
#' 
#' gDat <- groupedData(y ~ x | drug/CL, data = dfDat, FUN = mean, labels = list(x = "Concentration", y = "Viability"), units = list(x = "uM/l", y = "percentage killed"))
#'
#' @param szFileName  name to use as data source
#'
#' @return groupedData object

#'
#' @examples
#' gDat      <- initIC50(szFileName='~/myFile.csv')
#' fmMod1    <- fitModel(gDat)
#' outStats  <- gatherModelStats(gDat=gDat,fmMod1=fmMod1)
#' dfIC50    <- outStats$IC50
#' dfRMSE    <- outStats$RMSE
#' dfAUC     <- outStats$AUC
#' dfAUCtrap <- outStats$AUCtrap
#' 
#' The input data looks as follows:
#' 
#' Grouped Data: y ~ x | drug/CL
#'   x          y     CL drug maxc
#' 1 9 0.88079708 MC-CAR    1    2
#' 2 8 0.79139147 MC-CAR    1    2
#' 3 7 0.66075637 MC-CAR    1    2
#' 4 6 0.50000000 MC-CAR    1    2
#' 5 5 0.33924363 MC-CAR    1    2
#' 6 4 0.20860853 MC-CAR    1    2
#' 7 3 0.11920292 MC-CAR    1    2
#' 8 2 0.06496917 MC-CAR    1    2
#' 9 1 0.03444520 MC-CAR    1    2
#' 
#' @export
initIC50 <- function(szFileName='~/pub/gdsc_drug_sensitivity_raw_data_w5.csv'){
  library(foreach)
  library(doMC)
  registerDoMC()
  library(nlme)
  
  myDat <- read.delim(szFileName,sep=',')
  myDat2 <- myDat
  d1 <- as.matrix(myDat[,23:70])
  d2 <- matrix(as.integer(d1),nrow=dim(d1)[1],ncol=dim(d1)[2],byrow=F)
  myDat2[,23:70] <- d2
  myDat2$posctrl <- rowMeans(d2,na.rm=T)
  d1 <- as.matrix(myDat2[,71:102])
  d2 <- matrix(as.integer(d1),nrow=dim(d1)[1],ncol=dim(d1)[2],byrow=F)
  myDat2[,71:102] <- d2
  myDat2$negctrl <- rowMeans(d2,na.rm=T)
  myDat3 <- myDat2[,-c(23:102)]
  
  myDat3$raw_max <- (myDat3$raw_max-myDat3$negctrl)/(myDat3$posctrl-myDat3$negctrl)
  myDat3$raw2 <- (myDat3$raw2-myDat3$negctrl)/(myDat3$posctrl-myDat3$negctrl)
  myDat3$raw3 <- (myDat3$raw3-myDat3$negctrl)/(myDat3$posctrl-myDat3$negctrl)
  myDat3$raw4 <- (myDat3$raw4-myDat3$negctrl)/(myDat3$posctrl-myDat3$negctrl)
  myDat3$raw5 <- (myDat3$raw5-myDat3$negctrl)/(myDat3$posctrl-myDat3$negctrl)
  myDat3$raw6 <- (myDat3$raw6-myDat3$negctrl)/(myDat3$posctrl-myDat3$negctrl)
  myDat3$raw7 <- (myDat3$raw7-myDat3$negctrl)/(myDat3$posctrl-myDat3$negctrl)
  myDat3$raw8 <- (myDat3$raw8-myDat3$negctrl)/(myDat3$posctrl-myDat3$negctrl)
  myDat3$raw9 <- (myDat3$raw9-myDat3$negctrl)/(myDat3$posctrl-myDat3$negctrl)
  
  myDat3$raw_max[myDat3$raw_max>1] <- 1
  myDat3$raw2[myDat3$raw2>1] <- 1
  myDat3$raw3[myDat3$raw3>1] <- 1
  myDat3$raw4[myDat3$raw4>1] <- 1
  myDat3$raw5[myDat3$raw5>1] <- 1
  myDat3$raw6[myDat3$raw6>1] <- 1
  myDat3$raw7[myDat3$raw7>1] <- 1
  myDat3$raw8[myDat3$raw8>1] <- 1
  myDat3$raw9[myDat3$raw9>1] <- 1
  
  myDat3$raw_max[myDat3$raw_max<0] <- 0
  myDat3$raw2[myDat3$raw2<0] <- 0
  myDat3$raw3[myDat3$raw3<0] <- 0
  myDat3$raw4[myDat3$raw4<0] <- 0
  myDat3$raw5[myDat3$raw5<0] <- 0
  myDat3$raw6[myDat3$raw6<0] <- 0
  myDat3$raw7[myDat3$raw7<0] <- 0
  myDat3$raw8[myDat3$raw8<0] <- 0
  myDat3$raw9[myDat3$raw9<0] <- 0
  
  szDrug <- unique(myDat3$DRUG_ID)
  
  i5 <- c(9,7,5,3,1)
  i9 <- c(9,8,7,6,5,4,3,2,1)
  
  newDat <- foreach(i= 1:length(szDrug), .combine=rbind) %dopar%  {
    myLoc <- subset(myDat3, DRUG_ID==szDrug[i])
    dummy <- rep(NA,dim(myLoc)[1]*9)
    dfLoc <- data.frame(x=dummy,y=dummy,CL=dummy,drug=dummy,maxc=dummy)
    iCnt <- 1
    for(j in 1:dim(myLoc[1])){
      if(is.na(myLoc$raw6[j])==TRUE){
        dfLoc[iCnt:(iCnt+4),] <- as.data.frame(cbind(i5,c(myLoc[j,14:18]),as.character(myLoc$CELL_LINE_NAME[j]),myLoc$DRUG_ID[j],myLoc$MAX_CONC[j]))
        iCnt <- iCnt+5
      }else{
        dfLoc[iCnt:(iCnt+8),] <- as.data.frame(cbind(i9,c(myLoc[j,14:22]),as.character(myLoc$CELL_LINE_NAME[j]),myLoc$DRUG_ID[j],myLoc$MAX_CONC[j]))
        iCnt <- iCnt+9
      }
    }
    dfLoc <- dfLoc[1:iCnt-1,]
    return(dfLoc)
  }
  
  dfDat <- data.frame(x=unlist(newDat$x),y=unlist(newDat$y),CL=unlist(newDat$CL),drug=unlist(newDat$drug),maxc=unlist(newDat$maxc))
  dfDat$y <- 1-dfDat$y
  gDat <- groupedData(y ~ x | drug/CL, data = dfDat, FUN = mean, labels = list(x = "Concentration", y = "Viability"), units = list(x = "uM/l", y = "percentage killed"))
  return(gDat)
}


#' Curve fitting the data to the model
#'
#' Curve fitting the data prepared with \code{\link{initIC50}} to the model. The resultant object is further processed by \code{\link{gatherModelStats}}. 
#'
#' @param gDat  A \code{\link{groupedData}} object prepared in \code{\link{initIC50}}
#' @param vStart  The starting values of for the curve fit.
#' @param bLargeScale  Whether the covariance between the position and scale parameter on the cell line level are assumed to be correlated or not. In a large scale screen this is set to true as it further stabilizes the fit. In small bespoke screens this is set to false as the model otherwise struggles to converge.
#' @param bSilent  If set to true (default), suppresses status output. 
#'
#' @return fitted model, a nlme object
#'
#' @examples
#' gDat      <- initIC50(szFileName='~/myFile.csv')
#' fmMod1    <- fitModel(gDat)
#' outStats  <- gatherModelStats(gDat=gDat,fmMod1=fmMod1)
#' dfIC50    <- outStats$IC50
#' dfRMSE    <- outStats$RMSE
#' dfAUC     <- outStats$AUC
#' dfAUCtrap <- outStats$AUCtrap
#'
#' @export
fitModel <- function(gDat, vStart = c(8.886464, 1.495953 ),bLargeScale=TRUE,bSilent=TRUE){
  # start with sanity checks; catch most common mistakes early
  # first, do we have all requires columns?
  if(length(which(colnames(gDat)=='CL'))==0){
    stop('gDat is required to contain the cell lines names in a column with name CL')  
  }
  if(length(which(colnames(gDat)=='x'))==0){
    stop('gDat is required to contain the X-from concentration step names in a column with name x')  
  }
  if(length(which(colnames(gDat)=='y'))==0){
    stop('gDat is required to contain the relative kill (1-viability) in a column with name y')  
  }
  if(length(which(colnames(gDat)=='drug'))==0){
    stop('gDat is required to contain the drug in a column with name drug')  
  }
  if(length(which(colnames(gDat)=='maxc'))==0){
    stop('gDat is required to contain the maxc (maximum concentration) in a column with name maxc')  
  }
  # coding of relative kill (1-viability) correct; mean of gDat$y at gDat$x==9 should be higher than at min(gDat$x)
  tmp.minx <- min(gDat$x)
  tmp.whichxmin <- which(gDat$x == tmp.minx)
  tmp.whichxmax <- which(gDat$x == 9)
  tmp.ymin <- mean(na.omit(gDat$y[tmp.whichxmin]))
  tmp.ymax <- mean(na.omit(gDat$y[tmp.whichxmax]))
  if(tmp.ymin > tmp.ymax){
    stop('Coding of relative viabilities seems incorrect; note that y is defined as 1-viability.')
  }
  library(nlme)
  if(bLargeScale){
    fmv5 <- nlme(y ~ logist3(x, xmid, scal),fixed= xmid+scal~1,random=list(CL = pdSymm(xmid+scal~1),drug=pdDiag(xmid~1)),data = gDat, start=vStart,method='REML')
    if(!bSilent){
      summary(fmv5)
    }
    return(fmv5)
  }else{
    fmv5 <- nlme(y ~ logist3(x, xmid, scal), fixed = xmid + scal ~  1, random = list(CL = pdDiag(xmid + scal ~ 1), drug = pdDiag(xmid ~ 1)), data = gDat, start = vStart, method = "REML",control = nlmeControl(pnlsTol = 0.2, msVerbose = FALSE,tolerance=1e-4,returnObject = T))
    if(!bSilent){
      summary(fmv5)
    }
    return(fmv5)
  }
}

#' Gathering model statistics
#'
#' The data prepared with \code{\link{initIC50}}, and fitted with \code{\link{fitModel}} to the model, is used here to assemble output statiscs. As \code{\link{fitModel}} only parametrizes the model, the relevant statistics are here extracted. These values include the IC50, the RMSE, the model based AUC, and the trapezoid AUC.
#' When the series contains replicates, the normal AUC is define, but the trapezoid AUC is not.
#'
#' @param gDat  A \code{\link{groupedData}} object prepared in \code{\link{initIC50}}
#' @param fmMod1  The model fitted with \code{\link{fitModel}}
#'
#' @return A list with named items, holding the IC50s, AUCs, AUCtrap, RMSE, the gDat, the mcv, and the fitted model.
#'
#' @examples
#' gDat      <- initIC50(szFileName='~/myFile.csv')
#' fmMod1    <- fitModel(gDat)
#' outStats  <- gatherModelStats(gDat=gDat,fmMod1=fmMod1)
#' dfIC50    <- outStats$IC50
#' dfRMSE    <- outStats$RMSE
#' dfAUC     <- outStats$AUC
#' dfAUCtrap <- outStats$AUCtrap
#'
#' @export
gatherModelStats <- function(gDat,fmMod1){
  library(caTools)
  library(nlme)
  # do some sanity checking first
  if(length(which(colnames(gDat)=='maxc'))==0){
    warning('Missing maxc (maximum concentration), it should be a column in gDat')
    return(NULL)
  }
  
  mc1 <- coef(fmMod1)
  mx1 <- matrix(unlist(strsplit(rownames(mc1),split = '/')),nrow=dim(mc1)[1],ncol=2,byrow = T)
  mc1 <- cbind(mc1,mx1)
  
  colnames(mc1)[3:4] <- c('CL','drug')
  szLines <- unique(mc1$CL)
  szDrug <- unique(mc1$drug)
  mxXmid <- matrix(NA,nrow=length(unique(mc1[,3])),ncol=length(unique(mc1[,4])))
  
  colnames(mxXmid) <- szDrug
  rownames(mxXmid) <- szLines
  mxXauc <- mxXmid
  mxXic50 <- mxXauc
  mxXres <- mxXmid
  mxXtrap <- mxXmid
  #get IC50 and AUC from model parameters
  for(i in 1:dim(mxXmid)[2]){
    j <- which(mc1$drug == szDrug[i])
    j2 <- which(gDat$drug == szDrug[i])
    
    vLines <- unique(gDat$CL[j2])
    out <- lapply(1:length(vLines),function(ii){
      j2.1 <- which(gDat$CL[j2] == vLines[ii])
      j2.2 <- which(mc1$CL[j] == vLines[ii])
      v1 <- log(getConcentrationFromFold(xmid = mc1$xmid[j[j2.2]],fHighConc = max(gDat$maxc[j2[j2.1]]),iFold = 2))
      v2 <- 1-(getIntegral(9,mc1$xmid[j[j2.2]],mc1$scal[j[j2.2]])-getIntegral(min(gDat$x[j2[j2.1]]),mc1$xmid[j[j2.2]],mc1$scal[j[j2.2]]))/(max(gDat$x[j2[j2.1]])-min(gDat$x[j2[j2.1]]))
      yhat <- logist3(gDat$x[j2[j2.1]],mc1$xmid[j[j2.2]],mc1$scal[j[j2.2]])
      yres <- gDat$y[j2[j2.1]]-yhat
      v3 <- sqrt(mean(yres^2))
      v4 <- -trapz(gDat$x[j2[j2.1]],gDat$y[j2[j2.1]])/(max(gDat$x[j2[j2.1]])-min(gDat$x[j2[j2.1]]))
      c(v1,v2,v3,v4)
    })
    out2 <- do.call(rbind,out)
    mxXic50[as.character(vLines),i] <- out2[,1]
    mxXauc[as.character(vLines),i] <- out2[,2]
    mxXres[as.character(vLines),i] <- out2[,3]
    mxXtrap[as.character(vLines),i] <- out2[,4]
  }
  
  colnames(mxXic50) <- paste('DrugID',colnames(mxXic50),sep='')
  colnames(mxXauc) <- paste('DrugID',colnames(mxXauc),sep='')
  colnames(mxXres) <- paste('DrugID',colnames(mxXres),sep='')
  colnames(mxXtrap) <- paste('DrugID',colnames(mxXtrap),sep='')
  dfIC50 <- as.data.frame(mxXic50)
  dfAUC <- as.data.frame(mxXauc)
  dfRes <- as.data.frame(mxXres)
  dfAUCtrapz <- as.data.frame(mxXtrap)
  dfCL <- as.data.frame(ranef(fmMod1)$CL)

  lList <- list(IC50=dfIC50,AUC=dfAUC,CL=dfCL,RMSE=dfRes,AUCtrap=dfAUCtrapz,mcv=mc1,gDat=gDat,fmMod=fmMod1)
  return(lList)
}

#' @export
getIntegral <- function(x,xmid,scal){
  a <- xmid
  b <- scal
  return(b*log(exp(a/b)+exp(x/b))-a)
}

#' @export
logistInit3 <- function(mCall, LHS, data)
{   
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 3) {
    stop("Too few distinct input values to fit a logistic")
  }
  xmid <- NLSstClosestX(xy, 0.5 ) 
  scal <- NLSstClosestX(xy, 0.75 ) - xmid
  value <- c(xmid, scal)
  names(value) <- mCall[c("xmid", "scal")]
  value
}

#' @export
logist3 <- selfStart( ~ 1/(1 + exp(-(x - xmid)/(scal))), initial = logistInit3, parameters = c("xmid", "scal"))


#' @export
getConcentrationFromFold <- function(xmid,fHighConc=10,iFold=2,returnLogVal=1,iHighFold=9){
  fFrac <- 1/(iFold^(iHighFold-xmid))
  return(fHighConc*fFrac)
}
#' @export
getFoldFromConcentration <- function(fConc,fHighConc=10,iFold=2,returnLogVal=1,iHighFold=9){
  fFrac <- (fConc/fHighConc)
  return(iHighFold + (log(fFrac)/log(iFold)))
}


#' Plotting a single dose response curve
#'
#' The data prepared with \code{\link{initIC50}}, and fitted with \code{\link{fitModel}} to the model, is used here to assemble output statistics. As \code{\link{fitModel}} only parametrizes the model, the relevant statistics are here extracted. These values include the IC50, the RMSE, the model based AUC, and the trapezoid AUC.
#'
#' 
#'
#' @param gDat  A groupedData object prepared in \code{\link{initIC50}}
#' @param mc1  The internal coefficient of the nonlinear model, obtainable by \code{mc1 <- coef(fmMod1)}, using \code{\link{coef}} from \code{\link{nlme}}
#' @param szCL  The name/identifier of the cell line
#' @param szDrug  The name/identifier of the drug
#' @param dfOld  Optional argument for supplying reference values
#' @param szPrefix  Optional parameter for prefixing the plot title
#' @param szDrugName  Optional parameter for the drug name, printed as part of plot title
#'
#' @return No values are returned
#'
#' @examples
#' plotSingleResponse(gDat,mc1,'MC-CAR', 'DrugID1005',szPrefix='Screen 1',szDrugName='Platinum')
#'
#' @export
plotSingleResponse <- function(gDat,mc1,szCL, szDrug,dfOld=NULL,szPrefix='',szDrugName=''){
  library(scales)
  #no longer a function parameter
  iFold <- 2
  ##
  gDatLL <- subset(gDat, CL== szCL & drug == szDrug)
  if(dim(gDatLL)[1] ==0){
    plot(0,main=sprintf('No Data for Line %s on drug %s',szCL, szDrug))
  }else{
    mc1L <- mc1[sprintf('%s/%s',szCL,szDrug),]
    xx <- seq(from =0 , by=0.01,to=12)
    lxx <- log(getConcentrationFromFold(xx,fHighConc = gDatLL$maxc[1],iFold = iFold))
    lxd <- log(getConcentrationFromFold(gDatLL$x,fHighConc = gDatLL$maxc[1],iFold = iFold))
    yhat <- 1-logist3(xx,mc1L[1,1],mc1L[1,2])
    if(is.null(dfOld) == FALSE){
      j <- which(colnames(dfOld) == sprintf('DrugID%s',szDrug))
      lmax <- max(lxx,dfOld[szCL,j])
    }else{
      lmax <- max(lxx)
    }
    if((szPrefix == '') & (szDrugName=='')){
      szMain <- sprintf('%s on %s',szCL, szDrug)
    }else{
      if(nchar(szDrugName)==0){
        szMain <- sprintf('%s\n%s on %s',szPrefix,szCL, szDrug)
      }else{
        szMain <- sprintf('%s\n%s on %s',szPrefix,szCL, szDrugName)
      }
    }
    plot(lxx,yhat,col=1,type='l',ylab='Relative viability',xlab='log([Conc in uM])',xlim=c(min(lxx),lmax),ylim=c(-0.02,1.02),main=szMain,lwd=3)
    polygon(c(lxx,rev(lxx)),c(rep(0,length(yhat)),rev(yhat)),col=alpha('gray',0.15),border=NA)
    points(lxd,1-gDatLL$y,pch=20)
    grid()
    abline(h=0.5,col=alpha('gray',0.7),lty=5,lwd=3)
    abline(v=log(getConcentrationFromFold(xmid = mc1L[1,1],iFold = iFold,fHighConc = gDatLL$maxc[1])),col=2,lwd=3)
    if(is.null(dfOld) == FALSE){
      j <- which(colnames(dfOld) == sprintf('DrugID%s',szDrug))
      abline(v=dfOld[szCL,j],col=4,lwd=3)
    }
  }
}


#' Convenience function for dilution series coding
#'
#' A convenience function for getting the predictor values for custom dilution series. The function takes as argument a decreasing vector of fractions and converts them to the internal coding scheme of the curve fit, which defines 9 to be the highest screening concentration. The default assumption is a two fold dilution, and this function is there to accommodate other divisors. The series of 1, 1/2, 1/4, 1/8 would result in 9, 8, 7, 6.
#' 
#'
#' @param fDilutionSeries  A series of decreasing fractions.
#'
#' @return Series of predictor values
#'
#' @examples
#' % example of how to get x for a sqrt(2)-fold dilution per step
#' xSeries <- getXfromConcSeries(c(1, 0.7071, 0.5, 0.35355, 0.25))
#'
#' @export
getXfromConcSeries <- function (fDilutionSeries) # redefine func from djvMixedIC50
{
  if (min(fDilutionSeries) <= 0) {
    stop("Dilution series needs to be positive and larger than zero")
  }
  return((log(fDilutionSeries)/log(2)) + 9)
}
