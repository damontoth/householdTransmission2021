rm(list=ls())

QUU <- 1/(1-0.66)-1
QVV <- 1/(1-0.47)-1
QWW <- 1/(1-0.35)-1

QUV <- QVV
QVU <- QVV
QWV <- QVV
QVW <- QVV
QUW <- QWW
QWU <- QWW

sigV <- 0.5  #arbitrary
sigW <- 0.25 #arbitrary

xU <- QUU
xV <- QVU
xW <- QWU

fVV <- QVV/(sigV*xU)
fWW <- QWW/(sigW*xW)
fUV <- QUV/(sigV*xU)
fUW <- QUW/(sigW*xU)
fVW <- QVW/(sigW*xV)
fWV <- QWV/(sigV*xW)

#Assume within-group differences are due to susceptibility only:
getValsSusc <- function(){
	fUU <- 1; fVV <- 1; fWW <- 1
	xU <- 1; xV <- 1; xW <- 1
	sigU <- QUU/(fUU*xU)
	sigV <- QVV/(fVV*xV)
	sigW <- QWW/(fWW*xW)
	fUV <- QUV/(xU*sigV)
	fVU <- QVU/(xV*sigU)
	fUW <- QUW/(xU*sigW)
	fWU <- QWU/(xW*sigU)
	fVW <- QVW/(xV*sigW)
	fWV <- QWV/(xW*sigV)

	list(x = c(U=xU,V=xV,W=xW), sig = c(U=sigU, V=sigV, W=sigW), 
	     f = c(UU=fUU, VV=fVV, WW=fWW, UV=fUV, VU=fVU, UW=fUW, WU=fWU, VW=fVW, WV=fWV))
}

#Assume within-group differences are due to transmissibility only:
getValsTransm <- function(){
	fUU <- 1; fVV <- 1; fWW <- 1
	sigU <- 1; sigV <- 1; sigW <- 1
	xU <- QUU/(fUU*sigU)
	xV <- QVV/(fVV*sigV)
	xW <- QWW/(fWW*sigW)
	fUV <- QUV/(xU*sigV)
	fVU <- QVU/(xV*sigU)
	fUW <- QUW/(xU*sigW)
	fWU <- QWU/(xW*sigU)
	fVW <- QVW/(xV*sigW)
	fWV <- QWV/(xW*sigV)

	list(x = c(U=xU,V=xV,W=xW), sig = c(U=sigU, V=sigV, W=sigW), 
	     f = c(UU=fUU, VV=fVV, WW=fWW, UV=fUV, VU=fVU, UW=fUW, WU=fWU, VW=fVW, WV=fWV))
}

#Assume within-group differences are due to susceptibility & transmissibility in equal proportions:
getValsSuscTransm <- function(){
	fUU <- 1; fVV <- 1; fWW <- 1
	xU <- 1
	sigU <- QUU
	xV <- sqrt(QVV/QUU)
	sigV <- sqrt(QUU*QVV)
	xW <- sqrt(QWW/QUU)
	sigW <- sqrt(QUU*QWW)
	fUV <- QUV/(xU*sigV)
	fVU <- QVU/(xV*sigU)
	fUW <- QUW/(xU*sigW)
	fWU <- QWU/(xW*sigU)
	fVW <- QVW/(xV*sigW)
	fWV <- QWV/(xW*sigV)

	list(x = c(U=xU,V=xV,W=xW), sig = c(U=sigU, V=sigV, W=sigW), 
	     f = c(UU=fUU, VV=fVV, WW=fWW, UV=fUV, VU=fVU, UW=fUW, WU=fWU, VW=fVW, WV=fWV))
}

#Assume within-group differences are due to contact only:
getValsContact <- function(){
	xU <- 1; xV <- 1; xW <- 1
	sigU <- 1; sigV <- 1; sigW <- 1;
	fUU <- QUU
	fVV <- QVV
	fWW <- QWW
	fUV <- QUV
	fVU <- QVU
	fUW <- QUW
	fWU <- QWU
	fVW <- QVW
	fWV <- QWV

	list(x = c(U=xU,V=xV,W=xW), sig = c(U=sigU, V=sigV, W=sigW), 
	     f = c(UU=fUU, VV=fVV, WW=fWW, UV=fUV, VU=fVU, UW=fUW, WU=fWU, VW=fVW, WV=fWV))
}

#Assume within-group differences are due to susceptibility, transmissibility, & contact in equal proportions:
getValsSuscTransmContact <- function(){
	xU <- 1; sigU <- 1; fUU <- QUU
	xV <- (QVV/QUU)^(1/3)
	sigV <- xV
	fVV <- sigV*QUU
	xW <- (QWW/QUU)^(1/3)
	sigW <- xW
	fWW <- sigW*QUU
	fUV <- QUV/(xU*sigV)
	fVU <- QVU/(xV*sigU)
	fUW <- QUW/(xU*sigW)
	fWU <- QWU/(xW*sigU)
	fVW <- QVW/(xV*sigW)
	fWV <- QWV/(xW*sigV)

	list(x = c(U=xU,V=xV,W=xW), sig = c(U=sigU, V=sigV, W=sigW), 
	     f = c(UU=fUU, VV=fVV, WW=fWW, UV=fUV, VU=fVU, UW=fUW, WU=fWU, VW=fVW, WV=fWV))
}

stc <- getValsSuscTransmContact()
st <- getValsSuscTransm()
susc <- getValsSusc()
tran <- getValsTransm()
con <- getValsContact()

mechBarPlot <- function(p,lbl){
	par(cex=0.7, mar=c(3,1,3,1)+0.1)
	par(fig=c(0.05,0.5,0.5,1))
	barplot(p$x/max(p$x),main='Relative transmissibility')
	text(-0.1,1.15,lbl,xpd=TRUE,font=2,cex=1.3)
	par(fig=c(0.55,1,0.5,1), new=TRUE)
	barplot(p$sig/max(p$sig),main='Relative susceptibility')
	par(fig=c(0.05,1,0,0.5), new=TRUE)
	barplot(p$f/max(p$f), main='Relative amount of contact')
	
}

mechBarPlot(con,'A')
x11()
mechBarPlot(st,'B')
x11()
mechBarPlot(susc,'C')
x11()
mechBarPlot(tran,'D')
