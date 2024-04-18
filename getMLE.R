rm(list=ls())

# Available models:
# "Table2"  "Table3equalContact"  "Table3reciprocalTransmission"   "Table4"       
# "TableDgamma"  "TableDnoVariance"                    
# "TableErank2"  "TableErank3"  "TableErank4"  "TableErank5"  "TableErank6"                         
# "TableErank7"  "TableErank8"  "TableErank9"  "TableErank10"                        
# "TableFconsolidateVandU"   "TableFconsolidateWandV"   "TableFconsolidateWandVreciprocal"    
# "TableFequalHSARfromU"   "TableFequalHSARfromV"   "TableFequalHSARfromW"                
# "TableFequalHSARtoAndFromV"   "TableFequalHSARtoU"   "TableFequalHSARtoV"                  
# "TableFequalHSARtoW"   "TableFparsimonious"                  
# "TableGequalHSARallAges"   "TableGequalHSARallAgesAllSizes"      
# "TableGequalHSARexceptUtoU"   "TableGnoSizeEffect"                  
# "TableHfullTransmissionConsolidateAge"

model <- 'Table4' #choose one from above

source('specifyModel.R')

dataFile <- 'householdData2021.txt'

source('dataPreparationFunctions.R')
source('transmissionPrepFunctions.R')
source('MLEfunctions.R')

if(numTransmGroups == 3){
	getTransmissionPrep <- getTransmissionPrepUVW
	getH <- getH_UVW
}else if(numTransmGroups == 2){
	getTransmissionPrep <- getTransmissionPrepUV
	getH <- getH_UV
}else{
	getTransmissionPrep <- getTransmissionPrepU
	getH <- getH_U
}

df <- getHouseholdDataframe(dataFile, ageCutoffs)

cf1 <- getHouseholdComboFrequency(df, houseSize1)
cf2 <- getHouseholdComboFrequency(df, houseSize2)
cf3 <- getHouseholdComboFrequency(df, houseSize3)

numGrp <- length(ageCutoffs)+1

d1 <- getCalculationElements(cf1, numGrp)
d2 <- getCalculationElements(cf2, numGrp)
d3 <- getCalculationElements(cf3, numGrp)

p1 <- getTransmissionPrep(ageCutoffs, ageCutoffsTransm, d1)
p2 <- getTransmissionPrep(ageCutoffs, ageCutoffsTransm, d2)
p3 <- getTransmissionPrep(ageCutoffs, ageCutoffsTransm, d3)

initTest <- testsFun(epiFun1(epiPar), epiFun2(epiPar), epiFun3(epiPar), testsPar)

sol0 <- getSolIter(epiPar, epiFun1, epiFun2, epiFun3, testsPar, testsFun)
CIepi <- getCIepi(sol0$optEpi, sol0$optTests, sol0$optValue, epiFun1, epiFun2, epiFun3, testsFun)
CItests <- getCItests(sol0$optEpi, sol0$optTests, sol0$optValue, epiFun1, epiFun2, epiFun3, testsFun)

outTable <- getOutputTable(sol0$optEpi, sol0$optTests, CIepi, CItests)

print(outTable)
print(paste('Log likelihood =', round(sol0$optValue*100)/100))