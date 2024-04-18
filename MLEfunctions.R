dgammabinommixUVW <- function (yu, sizeu, yv, sizev, yw, sizew, q, k, sigv, sigw){
	if(k==0){
		out <- ifelse(sizeu+sizev+sizew == 0, 1, ifelse(yu==0 & yv==0 & yw==0, q, ifelse(yu==sizeu & yv==sizev & yw==sizew, 1-q, 0)))
	}else if(k==Inf){
		out <- dbinom(yu,sizeu,1-q) * dbinom(yv,sizev,1-q^sigv) * dbinom(yw,sizew,1-q^sigw)
	}else{
		qk <- q^(1/k)
		ju <- rep(0:yu, each=(yv+1)*(yw+1))
		jv <- rep(rep(0:yv, each=yw+1), yu+1)
		jw <- rep(0:yw, (yu+1)*(yv+1))
		if(yu == sizeu & yv == sizev & yw == sizew){
			ju <- ju[-1]; jv <- jv[-1]; jw <- jw[-1]
			out <- 1 + sum((-1)^(ju+jv+jw) * choose(yu,ju) * choose(yv,jv) * choose(yw,jw) * q*(qk+(ju + sigv*jv + sigw*jw)*(1-qk))^(-k))
		}else{
			out <- choose(sizeu,yu) * choose(sizev,yv) * choose(sizew,yw) * 
			sum((-1)^(ju+jv+jw) * choose(yu,ju) * choose(yv,jv) * choose(yw,jw) *
			q*(qk+(sizeu-yu+ju + sigv*(sizev-yv+jv) + sigw*(sizew-yw+jw))*(1-qk))^(-k))
		}
	}
	out
}

dgammabinommixUV <- function (yu, sizeu, yv, sizev, q, k, sigv){
	if(k==0){
		out <- ifelse(sizeu+sizev == 0, 1, ifelse(yu==0 & yv==0, q, ifelse(yu==sizeu & yv==sizev, 1-q, 0)))
	}else if(k==Inf){
		out <- dbinom(yu,sizeu,1-q) * dbinom(yv,sizev,1-q^sigv)
	}else{
		qk <- q^(1/k)
		ju <- rep(0:yu, each=yv+1)
		jv <- rep(0:yv, yu+1)
		if(yu == sizeu & yv == sizev){
			ju <- ju[-1]; jv <- jv[-1]
			out <- 1 + sum((-1)^(ju+jv) * choose(yu,ju) * choose(yv,jv) * q*(qk+(ju + sigv*jv)*(1-qk))^(-k))
		}else{
			out <- choose(sizeu,yu) * choose(sizev,yv) * 
			sum((-1)^(ju+jv) * choose(yu,ju) * choose(yv,jv) *
			q*(qk+(sizeu-yu+ju + sigv*(sizev-yv+jv))*(1-qk))^(-k))
		}
	}
	out
}

dgammabinommixU <- function (yu, sizeu, q, k){
	if(k==0){
		out <- ifelse(sizeu == 0, 1, ifelse(yu==0, q, ifelse(yu==sizeu, 1-q, 0)))
	}else if(k==Inf){
		out <- dbinom(yu,sizeu,1-q)
	}else{
		qk <- q^(1/k)
		ju <- 0:yu
		if(yu == sizeu){
			ju <- ju[-1]
			out <- 1 + sum((-1)^ju * choose(yu,ju) * q*(qk+ju*(1-qk))^(-k))
		}else{
			out <- choose(sizeu,yu) * 
			sum((-1)^(ju) * choose(yu,ju) *
			q*(qk+(sizeu-yu+ju)*(1-qk))^(-k))
		}
	}
	out
}

dgammabinommixUVW <- Vectorize(dgammabinommixUVW)
dgammabinommixUV <- Vectorize(dgammabinommixUV)
dgammabinommixU <- Vectorize(dgammabinommixU)

getH_UVW <- function(parTransm, p){
	ph <- p$phInit

	ph[p$pInd100YZ] <- dgammabinommixUVW(p$yu100,p$zu100,p$yv100,p$zv100,p$yw100,p$zw100,
						    	 1-parTransm['pUU'],parTransm['kU'],parTransm['sigUV'],parTransm['sigUW'])
	ph[p$pInd010YZ] <- dgammabinommixUVW(p$yu010,p$zu010,p$yv010,p$zv010,p$yw010,p$zw010,
						    	 1-parTransm['pVU'],parTransm['kV'],parTransm['sigVV'],parTransm['sigVW'])
	ph[p$pInd001YZ] <- dgammabinommixUVW(p$yu001,p$zu001,p$yv001,p$zv001,p$yw001,p$zw001,
						    	 1-parTransm['pWU'],parTransm['kW'],parTransm['sigWV'],parTransm['sigWW'])

	for(i in seq_along(p$xVal)){

		x <- p$xVal[i]
		j <- p$jStack[i,1]:p$jStack[i,2]
		
		ni <- p$niVals[j]

		uj <- (p$usCum[i]+1):p$usCum[i+1]
		vj <- (p$vsCum[i]+1):p$vsCum[i+1]
		wj <- (p$wsCum[i]+1):p$wsCum[i+1]

		phuSums <- rep(0,length(j))
		if(p$uCum[i+1] > p$uCum[i])
			phuSums[p$xuPos[j]] <- p$xuPosFrac[(p$uCum[i]+1):p$uCum[i+1]] * 
						     tapply(ph[p$indU1[uj]] * ph[p$indU2[uj]], rep(1:sum(p$xuPos[j]),ni[p$xuPos[j]]), sum)

		phvSums <- rep(0,length(j))
		if(p$vCum[i+1] > p$vCum[i])
			phvSums[p$xvPos[j]] <- p$xvPosFrac[(p$vCum[i]+1):p$vCum[i+1]] * 
						     tapply(ph[p$indV1[vj]] * ph[p$indV2[vj]], rep(1:sum(p$xvPos[j]),ni[p$xvPos[j]]), sum)

		phwSums <- rep(0,length(j))
		if(p$wCum[i+1] > p$wCum[i])
			phwSums[p$xwPos[j]] <- p$xwPosFrac[(p$wCum[i]+1):p$wCum[i+1]] * 
						     tapply(ph[p$indW1[wj]] * ph[p$indW2[wj]], rep(1:sum(p$xwPos[j]),ni[p$xwPos[j]]), sum)

		ph[p$phInd[j]] <- phuSums + phvSums + phwSums
	}
	ph
}

getH_UV <- function(parTransm, p){
	ph <- p$phInit

	ph[p$pInd10YZ] <- dgammabinommixUV(p$yu10,p$zu10,p$yv10,p$zv10,1-parTransm['pUU'],parTransm['kU'],parTransm['sigUV'])
	ph[p$pInd01YZ] <- dgammabinommixUV(p$yu01,p$zu01,p$yv01,p$zv01,1-parTransm['pVU'],parTransm['kV'],parTransm['sigVV'])
	
	for(i in seq_along(p$xVal)){

		x <- p$xVal[i]
		j <- p$jStack[i,1]:p$jStack[i,2]
		
		ni <- p$niVals[j]

		uj <- (p$usCum[i]+1):p$usCum[i+1]
		vj <- (p$vsCum[i]+1):p$vsCum[i+1]
		
		phuSums <- rep(0,length(j))
		phuSums[p$xuPos[j]] <- p$xuPosFrac[(p$uCum[i]+1):p$uCum[i+1]] * tapply(ph[p$indU1[uj]] * ph[p$indU2[uj]], rep(1:sum(p$xuPos[j]),ni[p$xuPos[j]]), sum)

		phvSums <- rep(0,length(j))
		phvSums[p$xvPos[j]] <- p$xvPosFrac[(p$vCum[i]+1):p$vCum[i+1]] * tapply(ph[p$indV1[vj]] * ph[p$indV2[vj]], rep(1:sum(p$xvPos[j]),ni[p$xvPos[j]]), sum)

		ph[p$phInd[j]] <- phuSums + phvSums
	}
	ph
}

getH_U <- function(parTransm, p){
	ph <- p$phInit

	ph[p$pInd1YZ] <- dgammabinommixU(p$yu1,p$zu1,1-parTransm['pUU'],parTransm['kU'])
		
	for(i in seq_along(p$xVal)){

		x <- p$xVal[i]
		j <- p$jStack[i,1]:p$jStack[i,2]
		
		ni <- p$niVals[j]

		uj <- (p$usCum[i]+1):p$usCum[i+1]
		
		phuSums <- rep(0,length(j))
		phuSums[p$xuPos[j]] <- p$xuPosFrac[(p$uCum[i]+1):p$uCum[i+1]] * tapply(ph[p$indU1[uj]] * ph[p$indU2[uj]], rep(1:sum(p$xuPos[j]),ni[p$xuPos[j]]), sum)

		ph[p$phInd[j]] <- phuSums
	}
	ph
}

getT <- function(parTransm, p){

	ph <- getH(parTransm, p)
	pt <- ph

	for(i in seq_along(p$yVal)){
		jy <- (p$yCum[i]+1):p$yCum[i+1]
		j <- (p$niCum[i]+1):p$niCum[i+1]

		prods <- ph[p$hInd[j]] * pt[p$tInd[j]]
		sums <- tapply(prods,rep(1:p$yCount[i],p$niValsT[jy]),sum)

		pt[p$ptInd[jy]] <- sums
	}
	pt
}

getPhsig <- function(pH, sig, k)
	as.numeric(ifelse(k < Inf, 1-(1+sig*((1-pH)^(-1/k)-1))^(-k), 1-(1-pH)^sig))

getSigXY <- function(pHouseXref, pHouseXY, k = 1)
	as.numeric(((1-pHouseXY)^(-1/k)-1) / ((1-pHouseXref)^(-1/k)-1))

zprod <- function(x,y) ifelse(y==0,0,x*y)

getM <- function(parImport, parTransm, d, p){

	pts <- getT(parTransm, p)[p$ptInds] * p$ptCoef

	dbProds <- d$ldbCoef + colSums( zprod(log(parImport), t(d$infGrp)) + zprod(log(1-parImport), t(d$uninfGrp)))
	
	tapply(exp(dbProds) * pts, d$mTerms, sum)
}

getParTransmHouseSize <- function(parTransm, houseSizeGroup){

	parT <- parTransm[setdiff(names(parTransm),c('sig34','sig5plus'))]

	names(parT)[names(parT) == 'pUU2'] <- 'pUU'
	names(parT)[names(parT) == 'pVU2'] <- 'pVU'
	names(parT)[names(parT) == 'pWU2'] <- 'pWU'	
	
	if(houseSizeGroup > 1){
		sigSize <- ifelse(houseSizeGroup == 2, parTransm['sig34'], parTransm['sig5plus'])
		parT['pUU'] <- getPhsig(parTransm['pUU2'], sigSize, parTransm['kU'])
		if(numTransmGroups > 1) parT['pVU'] <- getPhsig(parTransm['pVU2'], sigSize, parTransm['kV'])
		if(numTransmGroups > 2) parT['pWU'] <- getPhsig(parTransm['pWU2'], sigSize, parTransm['kW'])
	}
	parT
}

getM1 <- function(parImport, parTransm)
	getM(parImport, getParTransmHouseSize(parTransm, 1), d1, p1)

getM2 <- function(parImport, parTransm)
	getM(parImport, getParTransmHouseSize(parTransm, 2), d2, p2)

getM3 <- function(parImport, parTransm)
	getM(parImport, getParTransmHouseSize(parTransm, 3), d3, p3)

epiFun1 <- function(x) do.call(getM1, getEpiPar(x))
epiFun2 <- function(x) do.call(getM2, getEpiPar(x))
epiFun3 <- function(x) do.call(getM3, getEpiPar(x))

getSol <- function(initpar,fn) optim(initpar,fn,control=list(abstol=1e-10,maxit=10000))

getPtests <- function(m,phiGrp,piGrp,d){
	p <- 0
	if(all(c(phiGrp,piGrp) > 0) & all(c(phiGrp,piGrp) <= 1)){
		pr <- 1
		for(r in 1:ncol(phiGrp)){
			pr <- pr * t(phiGrp[,r] ^ t(d$truePos[,,r])) * t((1-phiGrp[,r]) ^ t(d$falseNeg[,,r])) *
				     t(piGrp[,r] ^ t(d$trueNeg[,,r])) * t((1-piGrp[,r]) ^ t(d$falsePos[,,r]))
		}
		x <- m * d$coef * apply(pr, 1, prod)
		p <- tapply(x,rep(seq_along(d$numTerms),d$numTerms),sum) * d$adjustCoef
	}
	p
}

getLLtests <- function(m,phiGrp,piGrp,cf,d) sum(cf$freq * log(getPtests(m,phiGrp,piGrp,d)))

testsFunBase <- function(m, x, cf, d) getLLtests(m, matrix(x[phiInds],4,byrow=TRUE),
	 					     	      matrix(x[piInds],4,byrow=TRUE),cf,d)
testsFun <- function(m1, m2, m3, x)
	testsFunBase(m1,x,cf1,d1) + testsFunBase(m2,x,cf2,d2) + testsFunBase(m3,x,cf3,d3)

getSolIter <- function(epiParInit,epiFun1,epiFun2,epiFun3,testsParInit,testsFun){

	epiPar <- epiParInit
	testsPar <- testsParInit

	parIter <- c(epiPar,testsPar)
	m1 <- epiFun1(epiPar); m2 <- epiFun2(epiPar); m3 <- epiFun3(epiPar)
	LLiter <- testsFun(m1, m2, m3, testsPar)
	stillIterating <- TRUE
	optValue <- Inf
	outMessage <- 'it worked!'

	while(stillIterating){
			
		solTests <- getSol(testsPar, function(x) -testsFun(m1, m2, m3, x))
		testsPar <- solTests$par
		if(solTests$value > optValue){
			stillIterating <- FALSE
			outMessage <- 'error: tests iteration got worse'
		}else{
			solEpi <- getSol(epiPar, function(x) ifelse(all(x>0), -testsFun(epiFun1(x), epiFun2(x), epiFun3(x), testsPar), Inf))
			epiPar <- solEpi$par
			parIter <- cbind(parIter,c(epiPar,testsPar))
			LLiter <- c(LLiter,solEpi$value)
			if(solEpi$value > solTests$value){
				stillIterating <- FALSE
				outMessage <- 'error: epi iteration got worse'
			}else{
				if(solEpi$value > optValue - 0.01){ stillIterating <- FALSE
				}else{ m1 <- epiFun1(epiPar); m2 <- epiFun2(epiPar); m3 <- epiFun3(epiPar)}

				optValue <- solEpi$value					
			}
		}
	}
	list(optEpi = epiPar, optTests = testsPar, optValue = optValue, outMessage = outMessage, parIter = parIter, LLiter = LLiter)		
}

getPhSize2 <- function(parEpi){
	p <- unlist(getEpiPar(parEpi)$parTransm)		
	out <- p['pUU2']
	if(numTransmGroups > 1){
		pUV2 <- getPhsig(p['pUU2'],p['sigUV'],p['kU'])
		pVV2 <- getPhsig(p['pVU2'],p['sigVV'],p['kV'])
	
		out <- c(p['pUU2'], p['pVU2'],
	  		   pUV2 = pUV2, pVV2 = pVV2)
	}
	if(numTransmGroups > 2){
		pUW2 <- getPhsig(p['pUU2'],p['sigUW'],p['kU'])
		pVW2 <- getPhsig(p['pVU2'],p['sigVW'],p['kV'])
		pWV2 <- getPhsig(p['pWU2'],p['sigWV'],p['kW'])
		pWW2 <- getPhsig(p['pWU2'],p['sigWW'],p['kW'])
		
		out <- c(p['pUU2'], p['pVU2'], p['pWU2'],
	  		   pUV2 = pUV2, pVV2 = pVV2, pWV2 = pWV2,
	  		   pUW2 = pUW2, pVW2 = pVW2, pWW2 = pWW2)
	}
	out
}

getInterval <- function(parName){
	nm <- substring(parName,1,1)
	if(nchar(parName) == 2 & nm == 'p'){ interval <- c(0,1)
	}else if(nm == 'p'){ interval <- c(0.001,0.999)
	}else if(nm == 'k'){ interval <- c(0,Inf)
	}else if(nm == 's'){ interval <- c(0,1000)}
	interval
}

getCIepi <- function(optEpi, optTests, optValue, epiFun1, epiFun2, epiFun3, testsFun){
	CIepi <- matrix(0,length(optEpi),2)
	solPhSize2 <- getPhSize2(optEpi)
	for(i in seq_along(optEpi)){
		parName <- names(optEpi)[i]

		ciFun <- function(p){
			par <- optEpi
			par[i] <- p
			if(parName %in% c('pUU2','pVU2','pWU2')){
				frmGrp <- substr(parName,2,2)
				if(numTransmGroups > 1) par[paste0('sig',frmGrp,'V')] <- getSigXY(p, solPhSize2[paste0('p',frmGrp,'V2')], k)
				if(numTransmGroups > 2) par[paste0('sig',frmGrp,'W')] <- getSigXY(p, solPhSize2[paste0('p',frmGrp,'W2')], k)
			}
			LLratioStat <- -2*(optValue + testsFun(epiFun1(par), epiFun2(par), epiFun3(par), optTests))
			LLratioStat - qchisq(0.95,1)
		}
		interval <- getInterval(parName)
		lowCI <- ifelse(ciFun(interval[1]) < 0, interval[1], uniroot(ciFun, c(interval[1],optEpi[i]),tol=1e-10)$root)
		highCI <- ifelse(ciFun(interval[2]) < 0, interval[2], uniroot(ciFun, c(optEpi[i],min(interval[2],1e5)),tol=1e-10)$root)
		CIepi[i,] <- c(lowCI, highCI)
	}
	rownames(CIepi) <- names(sol0$optEpi)
	CIepi
}

getCItests <- function(solParEpi, solParTests, solVal, getMfun1, getMfun2, getMfun3, getLLtestsFun){

	CI <- matrix(0,length(solParTests),2)
	colnames(CI) <- c('low','high')
	m1 <- getMfun1(solParEpi)
	m2 <- getMfun2(solParEpi)
	m3 <- getMfun3(solParEpi)

	for(i in seq_along(solParTests)){
		ciFun <- function(p){
			par <- solParTests
			par[i] <- p
			LLratioStat <- -2*(solVal + getLLtestsFun(m1, m2, m3, par))
			LLratioStat - qchisq(0.95,1)
		}
		nm <- names(solParTests)[i]
		lowCI <- uniroot(ciFun, c(0.01,solParTests[i]),tol=1e-10)$root
		highCI <- ifelse(ciFun(1) < 0, 1, uniroot(ciFun, c(solParTests[i],1),tol=1e-10)$root)

		CI[i,] <- c(lowCI, highCI)
	}
	CI
}

getOutputTable <- function(optEpi, optTests, CIepi, CItests){
	solPhSize2 <- getPhSize2(optEpi)
	solSigSize <- {}
	if(transmissionSizeModel == 'full') solSigSize <- optEpi[c('sig34','sig5plus')]
	mle <- c(optEpi[seq(npc)], solPhSize2, solSigSize, optTests)

	CIphSize2 <- cbind(getPhSize2(CIepi[,1]), getPhSize2(CIepi[,2]))
	CIsigSize <- {}
	if(transmissionSizeModel == 'full') CIsigSize <- CIepi[c('sig34','sig5plus'),]
	tableCI <- rbind(CIepi[seq(npc),], CIphSize2, CIsigSize, CItests)

	cbind(mle, tableCI)
}

getPval <- function(LL, optLL, df) pchisq(2*(optLL-LL), df=df, lower.tail=FALSE)
getAIC <- function(LL, numParam) 2*(numParam - LL)


