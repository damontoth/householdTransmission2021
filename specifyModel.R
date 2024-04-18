modelSettings <- read.table('resultTableSettings.txt', header=TRUE, row.names='resultName')

importationModel <- modelSettings[model,'importationModel']
ascertainmentModel <- modelSettings[model,'ascertainmentModel']
transmissionAgeModel <- modelSettings[model,'transmissionAgeModel']
transmissionSizeModel <- modelSettings[model,'transmissionSizeModel']
transmissibilityVarianceModel <- modelSettings[model,'transmissibilityVarianceModel']
transmissionAgeCutoff1 <- modelSettings[model,'transmissionAgeCutoff1']
transmissionAgeCutoff2 <- modelSettings[model,'transmissionAgeCutoff2']

if(importationModel == 'full'){
	ageCutoffsPcom <- c(6,13,19,25,45,65,80)
	parImportInit <- c(pA = 0.032, pB = 0.020, pC = 0.151, pD = 0.126,
			  	 pE = 0.067, pF = 0.067, pG = 0.052, pH = 1e-6)
}else{
	ageCutoffsPcom <- c(13,25,80)
	parImportInit <- c(pComA = 0.02, pComB = 0.14, pComC = 0.07, pComD = 1e-7)
}
npc <- length(parImportInit)

if(ascertainmentModel == 'full'){
	ageCutoffsPhi <- c(6,13,19,25,45,65,80)
	phiSinit <- c(phiSa = 0.16, phiSb = 0.36, phiSc = 0.53, phiSd = 0.62,
			  phiSe = 0.59, phiSf = 0.73, phiSg = 0.73, phiSh = 1-1e-6)
}else{
	ageCutoffsPhi <- c(13,45)
	phiSinit <- c(phiSa = 0.29, phiSb = 0.58, phiSc = 0.73)
}
nphi <- length(phiSinit)

testsPar <- c(phiSinit, phiAA = 0.72, phiAE = 0.993, phiAS = 0.87,
      	  piS = 0.992, piAA = 0.997, piAE = 0.981, piAS = 0.9985)

if(transmissibilityVarianceModel == 'exponential'){ k <- 1
}else if(transmissibilityVarianceModel == 'noVariance'){ k <- Inf
}else k <- NA

parTransmFix <- {}

if(!is.na(transmissionAgeCutoff1) & !is.na(transmissionAgeCutoff2)){
	ageCutoffsTransm <- c(transmissionAgeCutoff1, transmissionAgeCutoff2)
	
	if(transmissionAgeModel == 'full'){
		parTransmInit <- c(pUU2 = 0.68, pVU2 = 0.45, pWU2 = 0.42,
					 sigUV = 0.52, sigVV = 0.96, sigWV = 2.35,
					 sigUW = 0.25, sigVW = 1.09, sigWW = 1.23)

		getEpiPar <- function(x) list(parImport = x[parImportInds],
							parTransm = c(x[parTransmInds], parTransmFix))

	}else if(transmissionAgeModel == 'equalContact'){

		parTransmInit <- c(pUU2 = 0.58, pVU2 = 0.45, pWU2 = 0.53, sigUV = 1.05, sigUW = 0.9)
	
		getEpiPar <- function(x) list(parImport = x[parImportInds], 
						 	parTransm = c(x[parTransmInds], parTransmFix,
									  sigVV = as.numeric(x['sigUV']), sigWV = as.numeric(x['sigUV']),
									  sigVW = as.numeric(x['sigUW']), sigWW = as.numeric(x['sigUW'])))
					
	}else if(transmissionAgeModel == 'reciprocalTransmission'){

		parTransmInit <- c(pUU2 = 0.68, pVU2 = 0.49, pWU2 = 0.36, 
				    	 sigVV = 0.84, sigVW = 1.17, sigWW = 1.14)

		getEpiPar <- function(x) list(parImport = x[parImportInds], 
			     				parTransm = c(x[parTransmInds], parTransmFix,
						 			  sigUV = getSigXY(x['pUU2'], x['pVU2'], k),
				      				  sigUW = getSigXY(x['pUU2'], x['pWU2'], k),
						 			  sigWV = as.numeric(x['sigVW']) * getSigXY(x['pWU2'], x['pVU2'], k)))

	}else if(transmissionAgeModel == 'parsimonious'){
	
		parTransmInit <- c(pUU2 = 0.66, pVU2 = 0.47, pWU2 = 0.35)
		#parTransmInit <- c(pUU2 = 0.57, pVU2 = 0.30, pWU2 = 0.16)

		getEpiPar <- function(x) list(parImport = x[parImportInds],
							parTransm = c(x[parTransmInds], parTransmFix, sigVV = 1, sigVW = 1, sigWW = 1,
									  sigUV = getSigXY(x['pUU2'], x['pVU2'], k),
				      				  sigUW = getSigXY(x['pUU2'], x['pWU2'], k),
									  sigWV = getSigXY(x['pWU2'], x['pVU2'], k)))

	}else if(transmissionAgeModel == 'equalHSARfromU'){
		parTransmInit <- c(pUU2 = 0.68, pVU2 = 0.45, pWU2 = 0.42,
				     sigVV = 1, sigWV = 1, sigVW = 1, sigWW = 1)

		getEpiPar <- function(x) list(parImport = x[parImportInds], 
			     				parTransm = c(x[parTransmInds], parTransmFix, sigUV = 1, sigUW = 1))

	}else if(transmissionAgeModel == 'equalHSARfromV'){
		parTransmInit <- c(pUU2 = 0.68, pVU2 = 0.45, pWU2 = 0.42,
				     sigUV = 1, sigWV = 1, sigUW = 1, sigWW = 1)

		getEpiPar <- function(x) list(parImport = x[parImportInds], 
			     				parTransm = c(x[parTransmInds], parTransmFix, sigVV = 1, sigVW = 1))

	}else if(transmissionAgeModel == 'equalHSARfromW'){
		parTransmInit <- c(pUU2 = 0.68, pVU2 = 0.45, pWU2 = 0.42,
				     sigUV = 1, sigVV = 1, sigUW = 1, sigVW = 1)

		getEpiPar <- function(x) list(parImport = x[parImportInds], 
			     				parTransm = c(x[parTransmInds], parTransmFix, sigWV = 1, sigWW = 1))

	}else if(transmissionAgeModel == 'equalHSARtoU'){
		parTransmInit <- c(pUU2 = 0.59, sigUV = 0.9, sigVV = 0.8, sigWV = 1.2,
				     sigUW = 0.7, sigVW = 0.9, sigWW = 0.9)

		getEpiPar <- function(x) list(parImport = x[parImportInds], 
			     				parTransm = c(x[parTransmInds], parTransmFix,
						 			  pVU2 = as.numeric(x['pUU2']), pWU2 = as.numeric(x['pUU2'])))		

	}else if(transmissionAgeModel == 'equalHSARtoV'){
		parTransmInit <- c(pUU2 = 0.65, pVU2 = 0.36, pWU2 = 0.34,
				     sigUV = 0.5, sigUW = 0.2, sigVW = 0.8, sigWW = 1.05)

		getEpiPar <- function(x) list(parImport = x[parImportInds], 
			     				parTransm = c(x[parTransmInds], parTransmFix,
						 			  sigVV = as.numeric(x['sigUV']) * getSigXY(x['pVU2'], x['pUU2'], k),
						 			  sigWV = as.numeric(x['sigUV']) * getSigXY(x['pWU2'], x['pUU2'], k)))		

	}else if(transmissionAgeModel == 'equalHSARtoW'){
		parTransmInit <- c(pUU2 = 0.65, pVU2 = 0.36, pWU2 = 0.42,
				    sigUV = 0.6, sigVV = 1.2, sigWV = 1.8, sigUW = 0.4)

		getEpiPar <- function(x) list(parImport = x[parImportInds], 
			     				parTransm = c(x[parTransmInds], parTransmFix,
						 			  sigVW = as.numeric(x['sigUW']) * getSigXY(x['pVU2'], x['pUU2'], k),
						 			  sigWW = as.numeric(x['sigUW']) * getSigXY(x['pWU2'], x['pUU2'], k)))
	
	}else if(transmissionAgeModel == 'equalHSARtoAndFromV'){
		parTransmInit <- c(pUU2 = 0.68, pVU2 = 0.45, pWU2 = 0.42,
				     sigUW = 1, sigWW = 1)

		getEpiPar <- function(x) list(parImport = x[parImportInds], 
			     				parTransm = c(x[parTransmInds], parTransmFix, sigVV = 1, sigVW = 1,
						 			  sigUV = getSigXY(x['pUU2'], x['pVU2'], k),
						 			  sigWV = getSigXY(x['pWU2'], x['pVU2'], k)))
	}

}else if(!is.na(transmissionAgeCutoff1)){
	ageCutoffsTransm <- transmissionAgeCutoff1

	if(transmissionAgeModel == 'full'){
		parTransmInit <- c(pUU2 = 0.670, pVU2 = 0.403, sigUV = 0.328, sigVV = 1.36)

		getEpiPar <- function(x) list(parImport = x[parImportInds],
							parTransm = c(x[parTransmInds], parTransmFix))

	}else if(transmissionAgeModel == 'reciprocalTransmission'){
		parTransmInit <- c(pUU2 = 0.67, pVU2 = 0.403, sigVV = 1)

		getEpiPar <- function(x) list(parImport = x[parImportInds],
							parTransm = c(x[parTransmInds], parTransmFix, 
									  sigUV = getSigXY(x['pUU2'], x['pVU2'], k)))

	}else if(transmissionAgeModel == 'equalHSARexceptUtoU'){
		parTransmInit <- c(pUU2 = 0.67, pVU2 = 0.44)

		getEpiPar <- function(x) list(parImport = x[parImportInds],
							parTransm = c(x[parTransmInds], parTransmFix, sigVV = 1,
									  sigUV = getSigXY(x['pUU2'], x['pVU2'], k)))
	}

}else{
	ageCutoffsTransm = {}
	parTransmInit <- c(pUU2 = 0.54)
	getEpiPar <- function(x) list(parImport = x[parImportInds],
						parTransm = c(x[parTransmInds], parTransmFix))
}

if(is.na(k)){
	parTransmInit <- c(parTransmInit, kU = 1, kV = 1, kW = 1)
}else{
	parTransmFix <- c(parTransmFix, kU = k, kV = k, kW = k)
}

if(transmissionSizeModel == 'full'){
	#parTransmInit <- c(parTransmInit, sig34 = 0.43, sig5plus = 0.22)
	parTransmInit <- c(parTransmInit, sig34 = 0.51, sig5plus = 0.29)
}else{
	parTransmFix <- c(parTransmFix, sig34 = 1, sig5plus = 1)
}

epiPar <- c(parImportInit, parTransmInit)
nph <- length(parTransmInit)

numTransmGroups <- length(ageCutoffsTransm) + 1

ageCutoffs <- sort(unique(c(ageCutoffsPcom,ageCutoffsPhi,ageCutoffsTransm)))
getParInd <- function(x,ac) max(seq_along(ac)[c(0,ageCutoffs)[x] >= ac])
getParInd <- Vectorize(getParInd,'x')
parImportInds <- getParInd(seq_along(c(0,ageCutoffs)), c(0,ageCutoffsPcom))
parTransmInds <- npc + seq(nph)

phiSInds <- getParInd(seq_along(c(0,ageCutoffs)), c(0,ageCutoffsPhi))
phiInds <- c(phiSInds,rep(sequence(3,max(phiSInds)+1),each=length(phiSInds)))
piInds <- rep(sequence(4,max(phiSInds)+4),each=length(phiSInds))

houseSize1 <- c(1,2)
houseSize2 <- c(3,4)
houseSize3 <- c(5,Inf)