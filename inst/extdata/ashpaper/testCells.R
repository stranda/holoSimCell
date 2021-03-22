library(holoSimCell)

enm_model <- 10

load(file=paste0(system.file(package="holoSimCell"),"/extdata/landscapes/",enmScenarios$enms$rasterStackName[enm_model],".rda"))  #loads a landscape created from the rasterStack previously

plotCellsOnSuit(enmScenarios$refugeCellNum[[enm_model]],landscape)
