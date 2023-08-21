library(tidyverse)
library(randomForest)
library(terra)
library(pROC)

#load training data ###LINK TO FIGSHARE###
trainingData <- read.delim('training_data.txt')
presenceData <- trainingData[trainingData$occurrence=='present',]
backgroundData <- trainingData[trainingData$occurrence=='absent',]
worldStack <- rast('stack.tif')

fTable <- data.frame(table(presenceData$species))

#prepare training dataset for a target species
run_model <- function(target_species) {

#pull environmental data for all target species presences
targetPresences <- presenceData[presenceData$species == target_species,]
targetPresences$occurrence <- "present"
targetPresences <- targetPresences %>% dplyr::select(-species)

#attach background absence data
backgroundAbsences <- backgroundData %>% dplyr::select(-species)
backgroundAbsences$occurrence <- "absent"
trainingData <- rbind(targetPresences, backgroundAbsences)

#format training data
trainingData$occurrence <- factor(trainingData$occurrence)
trainingData$KG <- factor(trainingData$KG)
trainingData$plantClass <- factor(trainingData$plantClass)
trainingData$ecoflor <- factor(trainingData$ecoflor)
trainingData$geoClass <- factor(trainingData$geoClass)
trainingData$glac <- factor(trainingData$glac)
trainingData$soilClass <- factor(trainingData$soilClass)
trainingData$biom <- factor(trainingData$biom)
trainingData$forestType <- factor(trainingData$forestType)

#train model
num_presence <- length(trainingData$occurrence[trainingData$occurrence=='present'])
trained_model <- randomForest(occurrence ~ ., data=trainingData,  sampsize = c("present"=num_presence, "absent"=num_presence), ntree = 100)

#generate summary statistics
confusion <- data.frame(t(trained_model$confusion[,3]))
confusion$species <- target_species

importance <- data.frame(t(trained_model$importance))
importance$species <- target_species

ROC <- roc(trainingData$occurrence, trained_model$votes[,2])

ROCdf <- data.frame(species=rep(target_species, length(ROC$sensitivities)), specificities=(1-ROC$specificities), sensitivities=ROC$sensitivities)
auc <- data.frame(species=target_species, AUC=ROC$auc)

#run trained model on global raster stack
prediction <- predict(worldStack, trained_model)
prediction = prediction - 1

filename <- paste(gsub(" ","_", target_species), "tif", sep ='.')
writeRaster(prediction, filename, gdal=c("COMPRESS=DEFLATE"))

write.table(confusion, 'error.txt', sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
write.table(importance, 'importance.txt', sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
write.table(ROCdf, 'ROC.txt', sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
write.table(auc, 'AUC.txt', sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)

}

species2do <- fTable[fTable$Freq>4,]$Var1

for (sp in species2do) {
run_model(sp)
}