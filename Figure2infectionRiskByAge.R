rm(list=ls())

model <- 'TableFparsimonious'
source('specifyModel.R')
dataFile <- 'householdData2021.txt'

source('dataPreparationFunctions.R')
source('transmissionPrepFunctions.R')
source('MLEfunctions.R')

getTransmissionPrep <- getTransmissionPrepUVW
getH <- getH_UVW

df <- getHouseholdDataframeAgeOnly(dataFile, ageCutoffs)

cf1 <- getHouseholdComboFrequency(df, houseSize1)
cf2 <- getHouseholdComboFrequency(df, houseSize2)
cf3 <- getHouseholdComboFrequency(df, houseSize3)

numGrp <- length(ageCutoffs)+1

d1 <- getCalculationElementsAgeOnly(cf1, numGrp)
d2 <- getCalculationElementsAgeOnly(cf2, numGrp)
d3 <- getCalculationElementsAgeOnly(cf3, numGrp)

p1 <- getTransmissionPrep(ageCutoffs, ageCutoffsTransm, d1)
p2 <- getTransmissionPrep(ageCutoffs, ageCutoffsTransm, d2)
p3 <- getTransmissionPrep(ageCutoffs, ageCutoffsTransm, d3)

pCom <- c(pComA = 0.021,
	    pComB = 0.021,
 	    pComC = 0.137,
	    pComD = 0.137,
 	    pComE = 0.066,
	    pComF = 0.066,
	    pComG = 0.066,
	    pComH = 0)

epiPar <- c(pCom,
		pUU2 = 0.66,
		pVU2 = 0.47,
		pWU2 = 0.35,
		sig34 = 0.51,
		sig5plus = 0.29)

M1 <- epiFun1(epiPar)
M2 <- epiFun2(epiPar)
M3 <- epiFun3(epiPar)

tst1 <- c(M1)*ifelse(d1$iInf == 0, 0, d1$iInf/d1$cfRep)
ntRep1 <- rep(1:length(d1$numTerms),d1$numTerms)
infProbs1 <- matrix(0,nrow(cf1),numGrp)
infProb1 <- rep(0,numGrp)
for(j in 1:numGrp){
	infProbs1[,j] <- tapply(tst1[,j],ntRep1,sum)
	infProb1[j] <- sum(infProbs1[,j]*cf1$freq*cf1[,j])/sum(cf1$freq*cf1[,j])
}

tst2 <- c(M2)*ifelse(d2$iInf == 0, 0, d2$iInf/d2$cfRep)
ntRep2 <- rep(1:length(d2$numTerms),d2$numTerms)
infProbs2 <- matrix(0,nrow(cf2),numGrp)
infProb2 <- rep(0,numGrp)
for(j in 1:numGrp){
	infProbs2[,j] <- tapply(tst2[,j],ntRep2,sum)
	infProb2[j] <- sum(infProbs2[,j]*cf2$freq*cf2[,j])/sum(cf2$freq*cf2[,j])
}

tst3 <- c(M3)*ifelse(d3$iInf == 0, 0, d3$iInf/d3$cfRep)
ntRep3 <- rep(1:length(d3$numTerms),d3$numTerms)
infProbs3 <- matrix(0,nrow(cf3),numGrp)
infProb3 <- rep(0,numGrp)
for(j in 1:numGrp){
	infProbs3[,j] <- tapply(tst3[,j],ntRep3,sum)
	infProb3[j] <- sum(infProbs3[,j]*cf3$freq*cf3[,j])/sum(cf3$freq*cf3[,j])
}

infProbs <- rbind(infProbs1, infProbs2, infProbs3)
infProb <- rep(0,numGrp)
cfAge <- rbind(cf1, cf2, cf3)
for(j in 1:numGrp) infProb[j] <- sum(infProbs[,j]*cfAge$freq*cfAge[,j])/sum(cfAge$freq*cfAge[,j])

out <- cbind(infTot = infProb, infCom = pCom, infHouse = infProb-pCom,
		 housePortion = (infProb-pCom)/infProb)
rownames(out) <- c('0-5','6-12','13-18','19-24','25-44','45-64','65-79','80+')

barplot(t(out[,2:3]), col=c("black","grey"), border="white", space=0.04, font.axis=2, ylim=c(0,0.25),
	  xlab='age group',ylab='infection probability')


