getHouseholdDataframe <- function(dataFile, ageCutoffs){

	dRaw <- read.table(dataFile, header=TRUE)

	isV <- (dRaw$numVaccineDoses > 0)

	EuroimmunResult <- dRaw$EuroimmunResult
	EuroimmunResult[isV & dRaw$EuroimmunResult == 'P'] <- 'X'

	SiemensResult <- dRaw$SiemensResult
	SiemensResult[isV & dRaw$SiemensResult == 'P'] <- 'X'

	ageEx <- c(a = 3, b = 10, c = 16, d = 22, e = 35, f = 55, g = 73, h = 85)
 
	age <- ageEx[dRaw$ageCategory]

	numAgeGroups <- length(ageCutoffs) + 1

	ac <- c(0,ageCutoffs,Inf)
	
	ageGrp <- rep('x',nrow(dRaw))
	for(a in 1:numAgeGroups) ageGrp[age >= ac[a] & age < ac[a+1]] <- letters[a]

	res <- data.frame(id=dRaw$houseID, 
				results=paste0(dRaw$priorTest, dRaw$AbbottResult, EuroimmunResult, SiemensResult, ageGrp))

	ures <- unique(res$results)
	resAgg <- aggregate(list(res$results),list(id=res$id),function(x) tabulate(factor(x,levels=ures)))
	resList <- resAgg[,2]
	rlist <- lapply(resList,function(x) c(x,rep(0,length(ures)-length(x))))
	rmatrix <- do.call(rbind,rlist)
	rdf <- data.frame(rmatrix)
	names(rdf) <- ures

	rdf
}

getHouseholdDataframeAgeOnly <- function(dataFile, ageCutoffs){

	dRaw <- read.table(dataFile, header=TRUE)

	ageEx <- c(a = 3, b = 10, c = 16, d = 22, e = 35, f = 55, g = 73, h = 85)
 
	age <- ageEx[dRaw$ageCategory]

	numAgeGroups <- length(ageCutoffs) + 1

	ac <- c(0,ageCutoffs,Inf)
	
	ageGrp <- rep('x',nrow(dRaw))
	for(a in 1:numAgeGroups) ageGrp[age >= ac[a] & age < ac[a+1]] <- letters[a]

	res <- data.frame(id=dRaw$houseID, results=ageGrp)

	ures <- sort(unique(res$results))
	resAgg <- aggregate(list(res$results),list(id=res$id),function(x) tabulate(factor(x,levels=ures)))
	resList <- resAgg[,2]
	rlist <- lapply(resList,function(x) c(x,rep(0,length(ures)-length(x))))
	rmatrix <- do.call(rbind,rlist)
	rdf <- data.frame(rmatrix)
	names(rdf) <- ures

	rdf
}

getHouseholdComboFrequency <- function(df,sizeMinMax){
	NdfRaw <- df[rowSums(df) >= sizeMinMax[1] & rowSums(df) <= sizeMinMax[2],]
	Ndf <- NdfRaw[,colSums(NdfRaw)>0]

	cf <- aggregate(list(freq=rep(1,nrow(Ndf))),as.list(Ndf), sum)	
 	cf
}

getCalculationElements <- function(cf, numGrp){

	cff <- cf[,setdiff(names(cf),'freq')]

	result <- names(cff)
	nt <- nchar(result[1])

	grp <- letters[1:numGrp]

	grpCol <- matrix(FALSE,numGrp,length(result))
	grpCount <- matrix(0,nrow(cff),numGrp); colnames(grpCount) <- grp
	for(r in 1:numGrp){
		grpCol[r,] <- substr(result,nt,nt) == grp[r]
		grpCount[,r] <- rowSums(cff[,grpCol[r,],drop=FALSE])
	}

	#forward / backward products:
	nc <- ncol(cff)
	fw <- matrix(0,nrow(cff),nc)
	bk <- fw
	fw[,1] <- cff[,1]+1
	bk[,1] <- cff[,nc]+1
	for(i in 2:nc){
		fw[,i] <- (cff[,i]+1)*fw[,i-1]
		bk[,i] <- (cff[,nc-i+1]+1)*bk[,i-1] 
	}
	numTerms <- fw[,nc]

	cfRep <- apply(cff,2,function(x) rep(x,numTerms))

	vseq <- Vectorize(seq.default, c('from', 'to'))
	vrep <- Vectorize(rep.int)

	iInf <- matrix(0,sum(numTerms),nc)
	iInf[,1] <- rep(unlist(vseq(0,cff[,1])),rep(bk[,nc-1],fw[,1]))
	for(i in 2:(nc-1)) iInf[,i] <- rep(unlist(vrep(vseq(0,cff[,i]),fw[,i-1])),rep(bk[,nc-i], fw[,i]))
	iInf[,nc] <- unlist(vrep(vseq(0,cff[,nc]),fw[,nc-1]))

	colnames(iInf) <- names(cff)
	iUninf <- cfRep - iInf

	participantCol <- substr(result,1,nt-1) != paste(rep('X',nt-1),collapse='')

	testCombo <- gsub('P','O',gsub('N','O',result))
	uniqueTestCombo <- unique(testCombo)
	participantTotCol <- substr(uniqueTestCombo,1,nt-1) != paste(rep('X',nt-1),collapse='')

	grpColCombo <- matrix(FALSE, numGrp, length(uniqueTestCombo))
	for(r in 1:numGrp) grpColCombo[r,] <- substr(uniqueTestCombo,nt,nt) == grp[r]

	iInfTot <- matrix(0,nrow(iInf),length(uniqueTestCombo))
	iUninfTot <- iInfTot

	for(i in seq_along(uniqueTestCombo)){
		tcCol <- testCombo == uniqueTestCombo[i]
		iInfTot[,i] <- rowSums(matrix(iInf[,tcCol],nrow=nrow(iInf)))
		iUninfTot[,i] <- rowSums(matrix(iUninf[,tcCol],nrow=nrow(iInf)))
	}
	colnames(iInfTot) <- uniqueTestCombo
	colnames(iUninfTot) <- uniqueTestCombo

	lGrp <- matrix(0,nrow(iInf),numGrp)
	for(r in 1:numGrp) lGrp[,r] <- rowSums(iInf[,grpCol[r,],drop=FALSE])

	niL <- apply(lGrp+1,1,prod)
	infGrp <- matrix(0,sum(niL),numGrp)

	lProd <- function(col) apply(cbind(lGrp[,col]+1,1),1,prod)

	infGrp[,1] <- rep(sequence(lGrp[,1]+1,0), rep(lProd(2:numGrp), lGrp[,1]+1))
	if(numGrp > 2)
		for(r in 2:(numGrp-1)) 
			infGrp[,r] <- rep(sequence(rep(lGrp[,r]+1,lProd(1:(r-1))),0), rep(lProd((r+1):numGrp), lProd(1:r)))
	infGrp[,numGrp] <- sequence(rep(lGrp[,numGrp]+1,lProd(1:(numGrp-1))),0)

	lGrpRep <- apply(lGrp, 2, function(x) rep(x,niL))

	nGrp <- matrix(0,nrow(iInf),numGrp)
	for(r in 1:numGrp) nGrp[,r] <- rep(rowSums(cff[,grpCol[r,],drop=FALSE]),numTerms)

	nGrpRep <- apply(nGrp, 2, function(x) rep(x,niL))
	uninfGrp <- nGrpRep - infGrp

	ldbCoef <- rep(0,nrow(nGrpRep))

	for(r in 1:numGrp)
		ldbCoef <- ldbCoef + lchoose(nGrpRep[,r],infGrp[,r])

	mTerms <- rep(1:nrow(lGrp),niL)

	truePos <- trueNeg <- falsePos <- falseNeg <- array(0,dim=c(nrow(iInf),nt-1,numGrp))

	for(i in 1:(nt-1)){
		iResult <- substr(result,i,i)
		Pcol <- (iResult == 'P')
		Ncol <- (iResult == 'N')

		for(r in 1:numGrp){
			PgrpCol <- Pcol & grpCol[r,]
			NgrpCol <- Ncol & grpCol[r,]
			truePos[,i,r] <- rowSums(matrix(iInf[,PgrpCol],nrow=nrow(iInf)))
			falseNeg[,i,r] <- rowSums(matrix(iInf[,NgrpCol],nrow=nrow(iInf))) 
			falsePos[,i,r] <- rowSums(matrix(cfRep[,PgrpCol],nrow=nrow(iInf))) - truePos[,i,r]
			trueNeg[,i,r] <- rowSums(matrix(cfRep[,NgrpCol],nrow=nrow(iInf))) - falseNeg[,i,r]
		}
	}

	iInfParticipant <- iInf[,participantCol]
	iUninfParticipant <- iUninf[,participantCol]
	iInfTotParticipant <- iInfTot[,participantTotCol]
	iUninfTotParticipant <- iUninfTot[,participantTotCol]

	infCoef <- apply(factorial(iInfTotParticipant),1,prod) / apply(factorial(iInfParticipant),1,prod)
	uninfCoef <- apply(factorial(iUninfTotParticipant),1,prod) / apply(factorial(iUninfParticipant),1,prod)  

	hcoef <- rep(1,nrow(iInfTot))

	for(r in 1:numGrp){
		iITP <- as.matrix(iInfTotParticipant[,substr(colnames(iInfTotParticipant),nt,nt)==grp[r]])
		iUTP <- as.matrix(iUninfTotParticipant[,substr(colnames(iUninfTotParticipant),nt,nt)==grp[r]])
		iTP <- iITP + iUTP
	
		lastInf <- 0
		lastUninf <- 0
		for(i in 1:ncol(iTP)){
			hcoef <- hcoef * dhyper(iITP[,i], lGrp[,r]-lastInf, nGrp[,r]-lGrp[,r]-lastUninf, iTP[,i])
			lastInf <- lastInf + iITP[,i]
			lastUninf <- lastUninf + iUTP[,i]
		}
	}

	coef <- infCoef * uninfCoef * hcoef 

	facProd <- apply(factorial(cff),1,prod)
	resultTot <- matrix(0,nrow(cff),length(uniqueTestCombo))

	for(i in seq_along(uniqueTestCombo)){
		tcCol <- testCombo == uniqueTestCombo[i]
		resultTot[,i] <- rowSums(as.matrix(cff[,tcCol]))
	}
	facTotProd <- apply(factorial(resultTot),1,prod)
	adjustCoef <- facProd / facTotProd

	list(lGrp = lGrp, infGrp = infGrp, nGrp = nGrp, niL = niL, grpCount = grpCount,
	     nGrpRep = nGrpRep, lGrpRep = lGrpRep, ldbCoef = ldbCoef, uninfGrp = uninfGrp, 
	     mTerms = mTerms, truePos = truePos, falseNeg = falseNeg, trueNeg = trueNeg, falsePos = falsePos,
	     coef = coef, numTerms = numTerms, adjustCoef = adjustCoef)
}

getCalculationElementsAgeOnly <- function(cf, numGrp){
	cff <- cf[,setdiff(names(cf),'freq')]

	result <- names(cff)
	nt <- nchar(result[1])

	grp <- letters[1:numGrp]

	grpCol <- matrix(FALSE,numGrp,length(result))
	grpCount <- matrix(0,nrow(cff),numGrp); colnames(grpCount) <- grp
	for(r in 1:numGrp){
		grpCol[r,] <- substr(result,nt,nt) == grp[r]
		grpCount[,r] <- rowSums(cff[,grpCol[r,],drop=FALSE])
	}

	#forward / backward products:
	nc <- ncol(cff)
	fw <- matrix(0,nrow(cff),nc)
	bk <- fw
	fw[,1] <- cff[,1]+1
	bk[,1] <- cff[,nc]+1
	for(i in 2:nc){
		fw[,i] <- (cff[,i]+1)*fw[,i-1]
		bk[,i] <- (cff[,nc-i+1]+1)*bk[,i-1] 
	}
	numTerms <- fw[,nc]

	cfRep <- apply(cff,2,function(x) rep(x,numTerms))

	vseq <- Vectorize(seq.default, c('from', 'to'))
	vrep <- Vectorize(rep.int)

	iInf <- matrix(0,sum(numTerms),nc)
	iInf[,1] <- rep(unlist(vseq(0,cff[,1])),rep(bk[,nc-1],fw[,1]))
	for(i in 2:(nc-1)) iInf[,i] <- rep(unlist(vrep(vseq(0,cff[,i]),fw[,i-1])),rep(bk[,nc-i], fw[,i]))
	iInf[,nc] <- unlist(vrep(vseq(0,cff[,nc]),fw[,nc-1]))

	colnames(iInf) <- names(cff)

	lGrp <- nGrp <- matrix(0,nrow(iInf),numGrp)
	for(r in 1:numGrp){
		lGrp[,r] <- rowSums(iInf[,grpCol[r,],drop=FALSE])
		nGrp[,r] <- rep(rowSums(cff[,grpCol[r,],drop=FALSE]),numTerms)
	}

	niL <- apply(lGrp+1,1,prod)
	infGrp <- matrix(0,sum(niL),numGrp)

	lProd <- function(col) apply(cbind(lGrp[,col]+1,1),1,prod)

	infGrp[,1] <- rep(sequence(lGrp[,1]+1,0), rep(lProd(2:numGrp), lGrp[,1]+1))
	if(numGrp > 2)
		for(r in 2:(numGrp-1)) 
			infGrp[,r] <- rep(sequence(rep(lGrp[,r]+1,lProd(1:(r-1))),0), rep(lProd((r+1):numGrp), lProd(1:r)))
	infGrp[,numGrp] <- sequence(rep(lGrp[,numGrp]+1,lProd(1:(numGrp-1))),0)

	lGrpRep <- apply(lGrp, 2, function(x) rep(x,niL))
	nGrpRep <- apply(nGrp, 2, function(x) rep(x,niL))
	uninfGrp <- nGrpRep - infGrp

	ldbCoef <- rep(0,nrow(nGrpRep))
	for(r in 1:numGrp) ldbCoef <- ldbCoef + lchoose(nGrpRep[,r],infGrp[,r])

	mTerms <- rep(1:nrow(lGrp),niL)

	list(lGrp = lGrp, infGrp = infGrp, nGrp = nGrp, niL = niL, grpCount = grpCount,
	     nGrpRep = nGrpRep, lGrpRep = lGrpRep, ldbCoef = ldbCoef, uninfGrp = uninfGrp, 
	     mTerms = mTerms, numTerms = numTerms, iInf = iInf, cfRep = cfRep)
}

