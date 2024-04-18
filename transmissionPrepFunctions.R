getTransmissionPrepUVW <- function(ageCutoffs, ageCutoffsTransm, d){

	numGrp <- length(ageCutoffs) + 1
	grp <- letters[1:numGrp]

	numWgroups <- sum(ageCutoffs <= ageCutoffsTransm[1])
	wGroups <- letters[1:numWgroups] 
	vGroups <- letters[(numWgroups+1):sum(ageCutoffs <= ageCutoffsTransm[2])]

	lu <- lv <- lw <- nu <- nv <- nw <- rep(0,nrow(d$lGrp))
	infU <- infV <- infW <- rep(0,nrow(d$infGrp))

	for(r in 1:numGrp){
		if(grp[r] %in% vGroups){
			lv <- lv + d$lGrp[,r] 
			nv <- nv + d$nGrp[,r]
			infV <- infV + d$infGrp[,r]
		}else if(grp[r] %in% wGroups){
			lw <- lw + d$lGrp[,r] 
			nw <- nw + d$nGrp[,r]
			infW <- infW + d$infGrp[,r]
		}else{
			lu <- lu + d$lGrp[,r]
			nu <- nu + d$nGrp[,r]
			infU <- infU + d$infGrp[,r]
		}
	}

	luRep <- rep(lu,d$niL)
	lvRep <- rep(lv,d$niL)
	lwRep <- rep(lw,d$niL)
	nuRep <- rep(nu,d$niL)
	nvRep <- rep(nv,d$niL)
	nwRep <- rep(nw,d$niL)

	uninfU <- nuRep - infU
	uninfV <- nvRep - infV
	uninfW <- nwRep - infW

	numU <- numV <- numW <- 0

	for(g in 1:numGrp){
		if(grp[g] %in% vGroups){
			numV <- numV + d$grpCount[,g]
		}else if(grp[g] %in% wGroups){
			numW <- numW + d$grpCount[,g]
		}else{
			numU <- numU + d$grpCount[,g]
		}
	}
	Mmax <- max(numU+numV+numW)

	MuMax <- max(numU)
	Mu <- 0:MuMax
	Mv <- rep(0,MuMax+1)
	Mw <- matrix(0,MuMax+1,max(numV)+1)
	ind <- {}
	indX <- {}
	for(xu in 0:MuMax){
		Mv[xu+1] <- max(numV[numU >= xu])
		for(zu in 0:(MuMax-xu)){
			xvMax <- max(numV[numU >= xu+zu])
			for(xv in 0:xvMax){
				Mw[xu+1,xv+1] <- max(numW[numU >= xu & numV >= xv])
				for(zv in 0:(xvMax-xv)){
					xwMax <- max(numW[numU >= xu+zu & numV >= xv+zv])
					for(xw in 0:xwMax){
						if(zu==0 & zv==0) indX <- rbind(indX,c(xu,xv,xw))
						for(zw in 0:(xwMax-xw)){
							yu <- rep(0:zu,each=(zv+1)*(zw+1))
							yv <- rep(rep(0:zv,each=zw+1),zu+1)
							yw <- rep(0:zw,(zu+1)*(zv+1))
							ind <- rbind(ind,cbind(xu,xv,xw,yu,yv,yw,zu,zv,zw))
						}
					}
				}
			}
		}
	}

	indX <- ind[,1]+ind[,2]+ind[,3]
	indY <- ind[,4]+ind[,5]+ind[,6]

	MvMax <- max(Mv)
	MwMax <- max(Mw)
	xui <- cumsum(c(0,as.vector(table(ind[,1]))))
	zui <- matrix(0,MuMax+1,MuMax+2)
	xvi <- array(0,c(MuMax+1,MuMax+1,MvMax+2))
	zvi <- array(0,c(MuMax+1,MvMax+1,MvMax+2))
	xwi <- array(0,c(MuMax+1,MvMax+1,MwMax+2))

	for(i in 0:MuMax){
		zui[i+1,1:(MuMax+2-i)] <- cumsum(c(0,as.vector(table(ind[ind[,1]==i,7]))))
		for(j in 0:(MuMax-i)){
			xvi[i+1,j+1,1:(Mv[i+j+1]+2)] <- cumsum(c(0,as.vector(table(ind[ind[,1]==i & ind[,7]==j,2]))))	
		}
		for(j in 0:Mv[i+1]){
			zvi[i+1,j+1,1:(Mv[i+1]-j+2)] <- cumsum(c(0,choose(Mw[i+1,(j+1):(Mv[i+1]+1)]+3,3)*(1:(Mv[i+1]-j+1))))
			xwi[i+1,j+1,1:(Mw[i+1,j+1]+2)] <- cumsum(c(0,choose((Mw[i+1,j+1]+2):2,2)))
		}
	}
	zwi <- c(0,cumsum(1:(MwMax+1)))

	getIndXYZ <- function(xu,xv,xw,yu,yv,yw,zu,zv,zw){
		1 + xui[xu+1] + zui[xu+1 + (MuMax+1)*zu] + xvi[1 + xu + (MuMax+1)*zu + (MuMax+1)^2*xv] + 
		(zu+1)*zvi[1 + (xu+zu) + (MuMax+1)*xv + (MuMax+1)*(MvMax+1)*zv] +
		(zu+1)*(zv+1)*xwi[1 + (xu+zu) + (MuMax+1)*(xv+zv) + (MuMax+1)*(MvMax+1)*xw] + (zu+1)*(zv+1)*zwi[zw+1] + (zw+1)*(zv+1)*yu + (zw+1)*yv + yw
	}

	ptInds <- getIndXYZ(infU,infV,infW,luRep-infU,lvRep-infV,lwRep-infW,uninfU,uninfV,uninfW)

	ptCoef <- apply(choose(d$nGrpRep-d$infGrp, d$lGrpRep-d$infGrp), 1, prod) / 
		    choose(nuRep-infU, luRep-infU) / choose(nvRep-infV, lvRep-infV) / choose(nwRep-infW, lwRep-infW) 

	pInd00Z <- which(rowSums(ind[,1:6])==0)
	pIndX00 <- which(rowSums(ind[,4:9])==0)

	phInit <- rep(0,nrow(ind))
	phInit[pInd00Z] <- 1
	phInit[pIndX00] <- 1

	get_yu <- function(xu,xv,xw) ind[ind[,1]==xu & ind[,2]==xv & ind[,3]==xw, 4]
	get_yv <- function(xu,xv,xw) ind[ind[,1]==xu & ind[,2]==xv & ind[,3]==xw, 5]
	get_yw <- function(xu,xv,xw) ind[ind[,1]==xu & ind[,2]==xv & ind[,3]==xw, 6]
	get_zu <- function(xu,xv,xw) ind[ind[,1]==xu & ind[,2]==xv & ind[,3]==xw, 7]
	get_zv <- function(xu,xv,xw) ind[ind[,1]==xu & ind[,2]==xv & ind[,3]==xw, 8]  
	get_zw <- function(xu,xv,xw) ind[ind[,1]==xu & ind[,2]==xv & ind[,3]==xw, 9]

	yu100 <- get_yu(1,0,0)
	yu010 <- get_yu(0,1,0)
	yu001 <- get_yu(0,0,1)
	yv100 <- get_yv(1,0,0)
	yv010 <- get_yv(0,1,0)
	yv001 <- get_yv(0,0,1)
	yw100 <- get_yw(1,0,0)
	yw010 <- get_yw(0,1,0)
	yw001 <- get_yw(0,0,1)
	zu100 <- get_zu(1,0,0)
	zu010 <- get_zu(0,1,0)
	zu001 <- get_zu(0,0,1)
	zv100 <- get_zv(1,0,0)
	zv010 <- get_zv(0,1,0)
	zv001 <- get_zv(0,0,1)
	zw100 <- get_zw(1,0,0)
	zw010 <- get_zw(0,1,0)
	zw001 <- get_zw(0,0,1)

	pInd100YZ <- getIndXYZ(1,0,0,yu100,yv100,yw100,zu100,zv100,zw100)
	pInd010YZ <- getIndXYZ(0,1,0,yu010,yv010,yw010,zu010,zv010,zw010)
	pInd001YZ <- getIndXYZ(0,0,1,yu001,yv001,yw001,zu001,zv001,zw001)

	xCount <- tabulate(indX)
	xVal <- 2:(Mmax-1)
	xCum <- cumsum(c(0,xCount[xVal]))
	
	xu <- xv <- xw <- yu <- yv <- yw <- zu <- zv <- zw <- ni <- phInd <- rep(0,xCum[length(xCum)])
	
	xuPos <- xvPos <- xwPos <- rep(FALSE,length(xu))

	numU <- numV <- numW <- numUs <- numVs <- numWs <- rep(0,length(xVal))

	jStack <- matrix(0,length(xVal),2)

	for(i in seq_along(xVal)){
		x <- xVal[i]
		xRow <- (indX == x)
		j <- (xCum[i]+1):xCum[i+1]
		jStack[i,] <- c(xCum[i]+1,xCum[i+1])
		xu[j] <- ind[xRow,1]
		xv[j] <- ind[xRow,2]
		xw[j] <- ind[xRow,3]
		yu[j] <- ind[xRow,4]
		yv[j] <- ind[xRow,5]
		yw[j] <- ind[xRow,6]
		zu[j] <- ind[xRow,7]
		zv[j] <- ind[xRow,8]
		zw[j] <- ind[xRow,9]

		xuPos[j] <- (xu[j] > 0)
		xvPos[j] <- (xv[j] > 0)
		xwPos[j] <- (xw[j] > 0)

		numU[i] <- sum(xuPos[j])
		numV[i] <- sum(xvPos[j])
		numW[i] <- sum(xwPos[j])

		phInd[j] <- getIndXYZ(xu[j],xv[j],xw[j],yu[j],yv[j],yw[j],zu[j],zv[j],zw[j])

		ni[j] <- (yu[j]+1)*(yv[j]+1)*(yw[j]+1)
		xus <- rep(xu[j],ni[j])
		xvs <- rep(xv[j],ni[j])
		xws <- rep(xw[j],ni[j])

		numUs[i] <- sum(xus > 0)
		numVs[i] <- sum(xvs > 0)
		numWs[i] <- sum(xws > 0)
	}

	uCum <- cumsum(c(0,numU))
	vCum <- cumsum(c(0,numV))
	wCum <- cumsum(c(0,numW))

	usCum <- cumsum(c(0,numUs))
	vsCum <- cumsum(c(0,numVs))
	wsCum <- cumsum(c(0,numWs))

	indU1 <- indU2 <- tapplyIndexU <- rep(0,usCum[length(usCum)])
	indV1 <- indV2 <- tapplyIndexV <- rep(0,vsCum[length(vsCum)])
	indW1 <- indW2 <- tapplyIndexW <- rep(0,wsCum[length(wsCum)])

	xuPosFrac <- rep(0,uCum[length(uCum)])
	xvPosFrac <- rep(0,vCum[length(vCum)])
	xwPosFrac <- rep(0,wCum[length(wCum)])

	for(i in seq_along(xVal)){
		j <- (xCum[i]+1):xCum[i+1]	

		if(uCum[i+1] > uCum[i]) xuPosFrac[(uCum[i]+1):uCum[i+1]] <- xu[j][xu[j]>0]/xVal[i]
		if(vCum[i+1] > vCum[i]) xvPosFrac[(vCum[i]+1):vCum[i+1]] <- xv[j][xv[j]>0]/xVal[i]
		if(wCum[i+1] > wCum[i]) xwPosFrac[(wCum[i]+1):wCum[i+1]] <- xw[j][xw[j]>0]/xVal[i]

		xus <- rep(xu[j],ni[j])
		xvs <- rep(xv[j],ni[j])
		xws <- rep(xw[j],ni[j])
		yus <- rep(yu[j],ni[j])
		yvs <- rep(yv[j],ni[j])
		yws <- rep(yw[j],ni[j])
		zus <- rep(zu[j],ni[j])
		zvs <- rep(zv[j],ni[j])
		zws <- rep(zw[j],ni[j])

		ius <- rep(sequence(yu[j]+1,0),rep((yv[j]+1)*(yw[j]+1),yu[j]+1))
		ivs <- rep(sequence(rep(yv[j]+1,yu[j]+1),0),rep(yw[j]+1,(yu[j]+1)*(yv[j]+1)))
		iws <- sequence(rep(yw[j]+1,(yu[j]+1)*(yv[j]+1)),0) 

		ymiu <- yus-ius
		ymiv <- yvs-ivs
		ymiw <- yws-iws
		zmiu <- zus-ius
		zmiv <- zvs-ivs
		zmiw <- zws-iws

		if(usCum[i+1] > usCum[i]){
			uj <- (usCum[i]+1):usCum[i+1]
			gu <- (xus > 0)
			indU1[uj] <- getIndXYZ(xus[gu]-1,xvs[gu],xws[gu],ius[gu],ivs[gu],iws[gu],zus[gu],zvs[gu],zws[gu])
			indU2[uj] <- getIndXYZ(1,0,0,ymiu[gu],ymiv[gu],ymiw[gu],zmiu[gu],zmiv[gu],zmiw[gu])
			tapplyIndexU[uj] <- rep(1:sum(xu[j]>0),ni[j][xu[j]>0])
		}

		if(vsCum[i+1] > vsCum[i]){
			vj <- (vsCum[i]+1):vsCum[i+1]
			gv <- (xvs > 0)
			indV1[vj] <- getIndXYZ(xus[gv],xvs[gv]-1,xws[gv],ius[gv],ivs[gv],iws[gv],zus[gv],zvs[gv],zws[gv])
 			indV2[vj] <- getIndXYZ(0,1,0,ymiu[gv],ymiv[gv],ymiw[gv],zmiu[gv],zmiv[gv],zmiw[gv])
			tapplyIndexV[vj] <- rep(1:sum(xv[j]>0),ni[j][xv[j]>0])
		}

		if(wsCum[i+1] > wsCum[i]){
			wj <- (wsCum[i]+1):wsCum[i+1]
			gw <- (xws > 0)
			indW1[wj] <- getIndXYZ(xus[gw],xvs[gw],xws[gw]-1,ius[gw],ivs[gw],iws[gw],zus[gw],zvs[gw],zws[gw])
			indW2[wj] <- getIndXYZ(0,0,1,ymiu[gw],ymiv[gw],ymiw[gw],zmiu[gw],zmiv[gw],zmiw[gw])
			tapplyIndexW[wj] <- rep(1:sum(xw[j]>0),ni[j][xw[j]>0])
		}
	}

	niVals <- ni

	#The code below generates values to use in the "getT" function:

	yCount <- tabulate(indY)
	yVal <- 1:(Mmax-1)
	yCum <- cumsum(c(0,yCount[yVal]))
	
	ptInd <- rep(0,yCum[length(yCum)])

	niCount <- rep(0,length(yVal))

	for(i in seq_along(yVal)){
		y <- yVal[i]
		yRow <- (indY == y)

		xu <- ind[yRow,1]
		xv <- ind[yRow,2]
		xw <- ind[yRow,3]
		yu <- ind[yRow,4]
		yv <- ind[yRow,5]
		yw <- ind[yRow,6]
		zu <- ind[yRow,7]
		zv <- ind[yRow,8]
		zw <- ind[yRow,9]

		j <- (yCum[i]+1):yCum[i+1]
	
		ptInd[j] <- getIndXYZ(xu,xv,xw,yu,yv,yw,zu,zv,zw)

		niCount[i] <- sum((yu+1)*(yv+1)*(yw+1))
	}

	niCum <- cumsum(c(0,niCount))

	hInd <- tInd <- rep(0,niCum[length(niCum)])
	niValsT <- rep(0,yCum[length(yCum)])

	for(i in seq_along(yVal)){
	
		y <- yVal[i]
		yRow <- (indY == y)

		xu <- ind[yRow,1]
		xv <- ind[yRow,2]
		xw <- ind[yRow,3]
		yu <- ind[yRow,4]
		yv <- ind[yRow,5]
		yw <- ind[yRow,6]
		zu <- ind[yRow,7]
		zv <- ind[yRow,8]
		zw <- ind[yRow,9]

		iu <- rep(sequence(yu+1,0),rep((yv+1)*(yw+1),yu+1))
		iv <- rep(sequence(rep(yv+1,yu+1),0),rep(yw+1,(yv+1)*(yu+1)))
		iw <- sequence(rep(yw+1,(yu+1)*(yv+1)),0)

		ni <- (yu+1)*(yv+1)*(yw+1)

		xus <- rep(xu,ni)
		xvs <- rep(xv,ni)
		xws <- rep(xw,ni)
		yus <- rep(yu,ni)
		yvs <- rep(yv,ni)
		yws <- rep(yw,ni)
		zus <- rep(zu,ni)
		zvs <- rep(zv,ni)
		zws <- rep(zw,ni)

		j <- (niCum[i]+1):niCum[i+1]

		hInd[j] <- getIndXYZ(xus, xvs, xws, iu, iv, iw, zus, zvs, zws)
		tInd[j] <- getIndXYZ(iu, iv, iw, yus-iu, yvs-iv, yws-iw, zus-iu, zvs-iv, zws-iw)

		niValsT[(yCum[i]+1):yCum[i+1]] <- ni
	}

	list(phInit = phInit,
	     pInd100YZ = pInd100YZ,
	     yu100 = yu100, zu100 = zu100, yv100 = yv100, zv100 = zv100, yw100 = yw100, zw100 = zw100,
	     pInd010YZ = pInd010YZ,
	     yu010 = yu010, zu010 = zu010, yv010 = yv010, zv010 = zv010, yw010 = yw010, zw010 = zw010,
	     pInd001YZ = pInd001YZ,
	     yu001 = yu001, zu001 = zu001, yv001 = yv001, zv001 = zv001, yw001 = yw001, zw001 = zw001,
	     xVal = xVal,
	     jStack = jStack,
	     niVals = niVals,
	     usCum = usCum, vsCum = vsCum, wsCum = wsCum,
	     uCum = uCum, vCum = vCum, wCum = wCum,
	     xuPos = xuPos, xvPos = xvPos, xwPos = xwPos,
	     xuPosFrac = xuPosFrac, xvPosFrac = xvPosFrac, xwPosFrac = xwPosFrac,
	     indU1 = indU1, indU2 = indU2, indV1 = indV1, indV2 = indV2, indW1 = indW1, indW2 = indW2,
	     phInd = phInd,
	     yVal = yVal, 
	     yCum = yCum, niCum = niCum,
	     hInd = hInd, tInd = tInd,	
	     yCount = yCount,
	     niValsT = niValsT,
	     ptInd = ptInd,
	     ptInds = ptInds, ptCoef = ptCoef)
}

getTransmissionPrepUV <- function(ageCutoffs, ageCutoffsTransm, d){

	numGrp <- length(ageCutoffs) + 1
	grp <- letters[1:numGrp]

	numVgroups <- sum(ageCutoffs <= ageCutoffsTransm)
	vGroups <- letters[1:numVgroups] 

	lu <- lv <- nu <- nv <- rep(0,nrow(d$lGrp))

	infU <- infV <- rep(0,nrow(d$infGrp))

	for(r in 1:numGrp){
		if(grp[r] %in% vGroups){
			lv <- lv + d$lGrp[,r] 
			nv <- nv + d$nGrp[,r]
			infV <- infV + d$infGrp[,r]
		}else{
			lu <- lu + d$lGrp[,r]
			nu <- nu + d$nGrp[,r]
			infU <- infU + d$infGrp[,r]
		}
	}

	luRep <- rep(lu,d$niL)
	lvRep <- rep(lv,d$niL)
	nuRep <- rep(nu,d$niL)
	nvRep <- rep(nv,d$niL)

	uninfU <- nuRep - infU
	uninfV <- nvRep - infV

	numU <- numV <- 0

	for(g in 1:numGrp){
		if(grp[g] %in% vGroups){
			numV <- numV + d$grpCount[,g]
		}else{
			numU <- numU + d$grpCount[,g]
		}
	}

	Mmax <- max(numU+numV)

	MuMax <- max(numU)
	Mu <- 0:MuMax
	Mv <- rep(0,MuMax+1)

	ind <- {}
	for(xu in 0:MuMax){
		Mv[xu+1] <- max(numV[numU >= xu])
		for(zu in 0:(MuMax-xu)){
			xvMax <- max(numV[numU >= xu+zu])
			for(xv in 0:xvMax){
				for(zv in 0:(xvMax-xv)){
					yu <- rep(0:zu,each=zv+1)
					yv <- rep(0:zv,zu+1)
					ind <- rbind(ind,cbind(xu,xv,yu,yv,zu,zv))
				}
			}
		}
	}

	indX <- ind[,1]+ind[,2]
	indY <- ind[,3]+ind[,4]

	MvMax <- max(Mv)
	xui <- cumsum(c(0,as.vector(table(ind[,1]))))
	zui <- matrix(0,MuMax+1,MuMax+2)
	xvi <- matrix(0,MuMax+1,MvMax+2)

	for(i in 0:MuMax){
		zui[i+1,1:(MuMax+2-i)] <- cumsum(c(0,as.vector(table(ind[ind[,1]==i,5]))))
		xvi[i+1,1:(Mv[i+1]+2)] <- cumsum(c(0,choose((Mv[i+1]+2):2,2)))
	}
	zvi <- c(0,cumsum(1:(MvMax+1)))

	getIndXYZ <- function(xu,xv,yu,yv,zu,zv)
		1 + xui[xu+1] + zui[xu+1 + (MuMax+1)*zu] + (zu+1)*xvi[xu+zu+1 + (MuMax+1)*xv] + (zu+1)*zvi[zv+1] + (zv+1)*yu + yv

	ptInds <- getIndXYZ(infU,infV,luRep-infU,lvRep-infV,uninfU,uninfV)

	ptCoef <- apply(choose(d$nGrpRep-d$infGrp, d$lGrpRep-d$infGrp), 1, prod) / 
		    choose(nuRep-infU, luRep-infU) / choose(nvRep-infV, lvRep-infV) 

	pInd00Z <- which(rowSums(ind[,1:4])==0)
	pIndX00 <- which(rowSums(ind[,3:6])==0)

	phInit <- rep(0,nrow(ind))
	phInit[pInd00Z] <- 1
	phInit[pIndX00] <- 1

	get_yu <- function(xu,xv) ind[ind[,1]==xu & ind[,2]==xv, 3]
	get_yv <- function(xu,xv) ind[ind[,1]==xu & ind[,2]==xv, 4]
	get_zu <- function(xu,xv) ind[ind[,1]==xu & ind[,2]==xv, 5]
	get_zv <- function(xu,xv) ind[ind[,1]==xu & ind[,2]==xv, 6]  

	yu10 <- get_yu(1,0)
	yu01 <- get_yu(0,1)
	yv10 <- get_yv(1,0)
	yv01 <- get_yv(0,1)
	zu10 <- get_zu(1,0)
	zu01 <- get_zu(0,1)
	zv10 <- get_zv(1,0)
	zv01 <- get_zv(0,1)

	pInd10YZ <- getIndXYZ(1,0,yu10,yv10,zu10,zv10)
	pInd01YZ <- getIndXYZ(0,1,yu01,yv01,zu01,zv01)

	xCount <- tabulate(indX)
	xVal <- 2:(Mmax-1)
	xCum <- cumsum(c(0,xCount[xVal]))

	xu <- xv <- yu <- yv <- zu <- zv <- ni <- phInd <- rep(0,rev(xCum)[1])

	xuPos <- xvPos <- rep(FALSE,length(xu))

	numU <- numV <- numUs <- numVs <- rep(0,length(xVal))

	jStack <- matrix(0,length(xVal),2)

	for(i in seq_along(xVal)){
		x <- xVal[i]
		xRow <- (indX == x)
		j <- (xCum[i]+1):xCum[i+1]
		jStack[i,] <- c(xCum[i]+1,xCum[i+1])
		xu[j] <- ind[xRow,1]
		xv[j] <- ind[xRow,2]
		yu[j] <- ind[xRow,3]
		yv[j] <- ind[xRow,4]
		zu[j] <- ind[xRow,5]
		zv[j] <- ind[xRow,6]
	
		xuPos[j] <- (xu[j] > 0)
		xvPos[j] <- (xv[j] > 0)

		numU[i] <- sum(xuPos[j])
		numV[i] <- sum(xvPos[j])

		phInd[j] <- getIndXYZ(xu[j],xv[j],yu[j],yv[j],zu[j],zv[j])

		ni[j] <- (yu[j]+1)*(yv[j]+1)
		xus <- rep(xu[j],ni[j])
		xvs <- rep(xv[j],ni[j])

		numUs[i] <- sum(xus > 0)
		numVs[i] <- sum(xvs > 0)
	}

	uCum <- cumsum(c(0,numU))
	vCum <- cumsum(c(0,numV))

	usCum <- cumsum(c(0,numUs))
	vsCum <- cumsum(c(0,numVs))

	indU1 <- indU2 <- tapplyIndexU <- rep(0,rev(usCum)[1])
	indV1 <- indV2 <- tapplyIndexV <- rep(0,rev(vsCum)[1])

	xuPosFrac <- rep(0,rev(uCum)[1])
	xvPosFrac <- rep(0,rev(vCum)[1])

	for(i in seq_along(xVal)){
		j <- (xCum[i]+1):xCum[i+1]	

		xuPosFrac[(uCum[i]+1):uCum[i+1]] <- xu[j][xu[j]>0]/xVal[i]
		xvPosFrac[(vCum[i]+1):vCum[i+1]] <- xv[j][xv[j]>0]/xVal[i]

		xus <- rep(xu[j],ni[j])
		xvs <- rep(xv[j],ni[j])
		yus <- rep(yu[j],ni[j])
		yvs <- rep(yv[j],ni[j])
		zus <- rep(zu[j],ni[j])
		zvs <- rep(zv[j],ni[j])

		ius <- rep(sequence(yu[j]+1,0),rep(yv[j]+1,yu[j]+1))
		ivs <- sequence(rep(yv[j]+1,yu[j]+1),0) 

		ymiu <- yus-ius
		ymiv <- yvs-ivs
		zmiu <- zus-ius
		zmiv <- zvs-ivs

		uj <- (usCum[i]+1):usCum[i+1]
		vj <- (vsCum[i]+1):vsCum[i+1]

		gu <- (xus > 0)
		indU1[uj] <- getIndXYZ(xus[gu]-1,xvs[gu],ius[gu],ivs[gu],zus[gu],zvs[gu])
		indU2[uj] <- getIndXYZ(1,0,ymiu[gu],ymiv[gu],zmiu[gu],zmiv[gu])
		tapplyIndexU[uj] <- rep(1:sum(xu[j]>0),ni[j][xu[j]>0])

		gv <- (xvs > 0)
		indV1[vj] <- getIndXYZ(xus[gv],xvs[gv]-1,ius[gv],ivs[gv],zus[gv],zvs[gv])
 		indV2[vj] <- getIndXYZ(0,1,ymiu[gv],ymiv[gv],zmiu[gv],zmiv[gv])
		tapplyIndexV[vj] <- rep(1:sum(xv[j]>0),ni[j][xv[j]>0])
	}

	niVals <- ni

	#The code below generates values to use in the "getT" function:

	yCount <- tabulate(indY)
	yVal <- 1:(Mmax-1)
	yCum <- cumsum(c(0,yCount[yVal]))

	ptInd <- rep(0,rev(yCum)[1])

	niCount <- rep(0,length(yVal))

	for(i in seq_along(yVal)){
		y <- yVal[i]
		yRow <- (indY == y)

		xu <- ind[yRow,1]
		xv <- ind[yRow,2]
		yu <- ind[yRow,3]
		yv <- ind[yRow,4]
		zu <- ind[yRow,5]
		zv <- ind[yRow,6]

		j <- (yCum[i]+1):yCum[i+1]
	
		ptInd[j] <- getIndXYZ(xu,xv,yu,yv,zu,zv)

		niCount[i] <- sum((yu+1)*(yv+1))
	}

	niCum <- cumsum(c(0,niCount))

	hInd <- tInd <- rep(0,rev(niCum)[1])
	niValsT <- rep(0,rev(yCum)[1])

	for(i in seq_along(yVal)){
	
		y <- yVal[i]
		yRow <- (indY == y)

		xu <- ind[yRow,1]
		xv <- ind[yRow,2]
		yu <- ind[yRow,3]
		yv <- ind[yRow,4]
		zu <- ind[yRow,5]
		zv <- ind[yRow,6]

		iu <- rep(sequence(yu+1,0),rep(yv+1,yu+1))
		iv <- sequence(rep(yv+1,yu+1),0)

		ni <- (yu+1)*(yv+1)

		xus <- rep(xu,ni)
		xvs <- rep(xv,ni)
		yus <- rep(yu,ni)
		yvs <- rep(yv,ni)
		zus <- rep(zu,ni)
		zvs <- rep(zv,ni)

		j <- (niCum[i]+1):niCum[i+1]

		hInd[j] <- getIndXYZ(xus, xvs, iu, iv, zus, zvs)
		tInd[j] <- getIndXYZ(iu, iv, yus-iu, yvs-iv, zus-iu, zvs-iv)

		niValsT[(yCum[i]+1):yCum[i+1]] <- ni
	}

	list(phInit = phInit,
		pInd10YZ = pInd10YZ,
		yu10 = yu10, zu10 = zu10, yv10 = yv10, zv10 = zv10,
		pInd01YZ = pInd01YZ,
		yu01 = yu01, zu01 = zu01, yv01 = yv01, zv01 = zv01,
		xVal = xVal,
		jStack = jStack,
		niVals = niVals,
		usCum = usCum, vsCum = vsCum,
		uCum = uCum, vCum = vCum,
		xuPos = xuPos, xvPos = xvPos,
		xuPosFrac = xuPosFrac, xvPosFrac = xvPosFrac,
		indU1 = indU1, indU2 = indU2, indV1 = indV1, indV2 = indV2,
		phInd = phInd,
		yVal = yVal, 
		yCum = yCum, niCum = niCum,
		hInd = hInd, tInd = tInd,	
		yCount = yCount,
		niValsT = niValsT,
		ptInd = ptInd,
		ptInds = ptInds, ptCoef = ptCoef)
}

getTransmissionPrepU <- function(ageCutoffs, ageCutoffsTransm, d){

	numGrp <- length(ageCutoffs) + 1
	grp <- letters[1:numGrp]

	lu <- nu <- rep(0,nrow(d$lGrp))
	
	infU <- rep(0,nrow(d$infGrp))

	for(r in 1:numGrp){
		lu <- lu + d$lGrp[,r]
		nu <- nu + d$nGrp[,r]
		infU <- infU + d$infGrp[,r]
	}

	luRep <- rep(lu,d$niL)
	nuRep <- rep(nu,d$niL)

	uninfU <- nuRep - infU

	numU <- 0
	for(g in 1:numGrp){
		numU <- numU + d$grpCount[,g]
	}

	Mmax <- max(numU)
	Mu <- 0:Mmax

	ind <- {}
	for(xu in 0:Mmax){
		for(zu in 0:(Mmax-xu)){
			yu <- 0:zu
			ind <- rbind(ind,cbind(xu,yu,zu))
		}
	}

	indX <- ind[,1]
	indY <- ind[,2]

	xui <- cumsum(c(0,as.vector(table(ind[,1]))))
	zui <- matrix(0,Mmax+1,Mmax+2)

	for(i in 0:Mmax){
		zui[i+1,1:(Mmax+2-i)] <- cumsum(c(0,as.vector(table(ind[ind[,1]==i,3]))))
	}

	getIndXYZ <- function(xu,yu,zu)
		1 + xui[xu+1] + zui[xu+1 + (Mmax+1)*zu] + yu

	ptInds <- getIndXYZ(infU,luRep-infU,uninfU)

	ptCoef <- apply(choose(d$nGrpRep-d$infGrp, d$lGrpRep-d$infGrp), 1, prod) / 
		    choose(nuRep-infU, luRep-infU) 

	pInd00Z <- which(rowSums(ind[,1:2])==0)
	pIndX00 <- which(rowSums(ind[,2:3])==0)

	phInit <- rep(0,nrow(ind))
	phInit[pInd00Z] <- 1
	phInit[pIndX00] <- 1

	get_yu <- function(xu) ind[ind[,1]==xu, 2]
	get_zu <- function(xu) ind[ind[,1]==xu, 3]  

	yu1 <- get_yu(1)
	zu1 <- get_zu(1)

	pInd1YZ <- getIndXYZ(1,yu1,zu1)

	xCount <- tabulate(indX)
	xVal <- 2:(Mmax-1)
	xCum <- cumsum(c(0,xCount[xVal]))

	xu <- yu <- zu <- ni <- phInd <-  rep(0,rev(xCum)[1])

	xuPos <- rep(FALSE,length(xu))

	numU <- numUs <- rep(0,length(xVal))

	jStack <- matrix(0,length(xVal),2)

	for(i in seq_along(xVal)){
		x <- xVal[i]
		xRow <- (indX == x)
		j <- (xCum[i]+1):xCum[i+1]
		jStack[i,] <- c(xCum[i]+1,xCum[i+1])
		xu[j] <- ind[xRow,1]
		yu[j] <- ind[xRow,2]
		zu[j] <- ind[xRow,3]
	
		xuPos[j] <- (xu[j] > 0)

		numU[i] <- sum(xuPos[j])

		phInd[j] <- getIndXYZ(xu[j],yu[j],zu[j])

		ni[j] <- (yu[j]+1)
		xus <- rep(xu[j],ni[j])

		numUs[i] <- sum(xus > 0)
	}

	uCum <- cumsum(c(0,numU))
	usCum <- cumsum(c(0,numUs))

	indU1 <- indU2 <- tapplyIndexU <- rep(0,rev(usCum)[1])

	xuPosFrac <- rep(0,rev(uCum)[1])

	for(i in seq_along(xVal)){
		j <- (xCum[i]+1):xCum[i+1]	

		xuPosFrac[(uCum[i]+1):uCum[i+1]] <- xu[j][xu[j]>0]/xVal[i]

		xus <- rep(xu[j],ni[j])
		yus <- rep(yu[j],ni[j])
		zus <- rep(zu[j],ni[j])

		ius <- sequence(yu[j]+1,0)

		ymiu <- yus-ius
		zmiu <- zus-ius

		uj <- (usCum[i]+1):usCum[i+1]

		gu <- (xus > 0)
		indU1[uj] <- getIndXYZ(xus[gu]-1,ius[gu],zus[gu])
		indU2[uj] <- getIndXYZ(1,ymiu[gu],zmiu[gu])
		tapplyIndexU[uj] <- rep(1:sum(xu[j]>0),ni[j][xu[j]>0])
	}

	niVals <- ni

	#The code below generates values to use in the "getT" function:

	yCount <- tabulate(indY)
	yVal <- 1:(Mmax-1)
	yCum <- cumsum(c(0,yCount[yVal]))

	ptInd <- rep(0,rev(yCum)[1])

	niCount <- rep(0,length(yVal))

	for(i in seq_along(yVal)){
		y <- yVal[i]
		yRow <- (indY == y)

		xu <- ind[yRow,1]
		yu <- ind[yRow,2]
		zu <- ind[yRow,3]

		j <- (yCum[i]+1):yCum[i+1]
	
		ptInd[j] <- getIndXYZ(xu,yu,zu)

		niCount[i] <- sum(yu+1)
	}

	niCum <- cumsum(c(0,niCount))

	hInd <- tInd <- rep(0,rev(niCum)[1])
	niValsT <- rep(0,rev(yCum)[1])

	for(i in seq_along(yVal)){
	
		y <- yVal[i]
		yRow <- (indY == y)

		xu <- ind[yRow,1]
		yu <- ind[yRow,2]
		zu <- ind[yRow,3]

		iu <- sequence(yu+1,0)

		ni <- (yu+1)

		xus <- rep(xu,ni)
		yus <- rep(yu,ni)
		zus <- rep(zu,ni)

		j <- (niCum[i]+1):niCum[i+1]

		hInd[j] <- getIndXYZ(xus, iu, zus)
		tInd[j] <- getIndXYZ(iu, yus-iu, zus-iu)

		niValsT[(yCum[i]+1):yCum[i+1]] <- ni
	}

	list(phInit = phInit,
		pInd1YZ = pInd1YZ,
		yu1 = yu1, zu1 = zu1,
		xVal = xVal,
		jStack = jStack,
		niVals = niVals,
		usCum = usCum,
		uCum = uCum,
		xuPos = xuPos,
		xuPosFrac = xuPosFrac,
		indU1 = indU1, indU2 = indU2,
		phInd = phInd,
		yVal = yVal, 
		yCum = yCum, niCum = niCum,
		hInd = hInd, tInd = tInd,	
		yCount = yCount,
		niValsT = niValsT,
		ptInd = ptInd,
		ptInds = ptInds, ptCoef = ptCoef)
}

