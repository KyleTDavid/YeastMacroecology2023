library(tidyverse)
library(readxl)
library(CoordinateCleaner)

#Filter and Assemble Coordinates

#GlobalFungi occurrences (https://doi.org/10.1038/s41597-020-0567-7)
#taxonomic key
taxa <- read.delim("Saccharomycetes_OTUTAB/Saccharomycetes_SHs_tax.txt", header=FALSE)
#sample collection metadata
samples <- read.delim("Saccharomycetes_OTUTAB/Saccharomycetes_REL4_TABLE_TAB_SH.txt.metadata.txt", sep='\t', na.strings = "NA_")
#format occurrences table
gf <- read.delim("Saccharomycetes_OTUTAB/Saccharomycetes_REL4_TABLE_TAB_SH.txt.otutab.txt", sep='\t', header=TRUE)
gf <- pivot_longer(gf, 2:ncol(gf))
gf <- gf[gf$value > 0,][1:2]
#attach metadata
gf <- merge(gf, samples, by.x = 'name', by.y = 'id')
gf$species <- taxa[match(gf[['SH']], taxa[['V1']]), 'V9']
#remove sp.s
gf <- gf[!grepl("sp.", gf$species),]
#select relevant columns
gf <- gf %>% select(species, longitude, latitude)
#add source
gf$source <- "GlobalFungi v4"

#GBIF occurrences (https://doi.org/10.15468/dl.n4fkqs)
#input is GBIF coordinates file
gbif <- read.csv("0216703-220831081235567.csv", sep='\t', header=TRUE, comment.char='', quote="")
#remove rows without species name
gbif = gbif[which(!is.na(gbif$species)),] 
gbif = gbif[which(gbif$species!=""),] 
#convert latitude, longitude and coordinateUncertainty to numeric values
gbif[,c(2:4)] = apply(gbif[,c(2:4)], 2, function(x) as.numeric(as.character(x)))
#filter entries if reported coordinate uncertainty is greater than 1km
gbif <- gbif[is.na(gbif$coordinateUncertaintyInMeters) | gbif$coordinateUncertaintyInMeters < 1000,]
#filter out unnecessary columns 
gbif <- gbif %>% select(species, decimalLongitude, decimalLatitude)
gbif <- rename(gbif, longitude=decimalLongitude, latitude=decimalLatitude)
gbif$source <- "GBIF.org (14 December 2022) GBIF Occurrence Download https://doi.org/10.15468/dl.n4fkqs"

#samples from Peris et al. (https://doi.org/10.6084/m9.figshare.17185874)
peris <- read_excel('TableS1.xlsx', skip = 1)
#only wild collected samples
peris <- peris[peris$`Wild/Anthropic`=='Wild',]
peris <- peris %>% select(Species, Longitude, Latitude)
peris <- na.omit(peris)
peris$source <- "Peris, David, et al. \"Macroevolutionary diversity of traits and genomes in the model yeast genus Saccharomyces.\" bioRxiv (2022)."
peris <- rename(peris, species=Species, latitude=Latitude, longitude=Longitude)

#samples from Spurley et al. (https://doi.org/10.6084/m9.figshare.c.5599392)
spurley <- read_excel('Spurley_S1.xlsx', na = 'NA')
spurley <- spurley[spurley$Subphylum == 'Saccharomycotina',]
spurley <- spurley %>% select(OTU, Lattitude, Longitude) %>% rename(species=OTU, latitude=Lattitude, longitude=Longitude)
spurley = na.omit(spurley)
spurley$source <- "Spurley, William J., et al. \"Substrate, temperature, and geographical patterns among nearly 2000 natural yeast isolates.\" Yeast 39.1-2 (2022): 55-68. APA"
                      
#combine and fix synonyms (Dataset S4)
speciesCoord <- rbind(gf,gbif,peris,spurley)
synonyms <- read.delim("synonyms.txt", na.strings = "")
synonyms <- synonyms[!is.na(synonyms$Species),]
key <- synonyms$Species[match(unlist(speciesCoord$species), synonyms$Synonym)]
speciesCoord$species <- if_else(!is.na(key), key, speciesCoord$species)
                      
#remove sp.s
speciesCoord <- speciesCoord[!grepl("sp.", speciesCoord$species),]
speciesCoord <- unique(speciesCoord)

#don't care about industrial hybrids
speciesCoord <- speciesCoord[!(speciesCoord$species %in% c('Saccharomyces bayanus', 'Saccharomyces pastorianus', 'Saccharomyces uvarum')),]

#resolution of coordinates must be at least 0.01
length_lat = nchar(sub("^[^.]*", "", speciesCoord$latitude))-1
length_lon = nchar(sub("^[^.]*", "", speciesCoord$longitude))-1
speciesCoord = speciesCoord[which(length_lat>=2 & length_lon>=2),] 

#clean coordinates
flags <- clean_coordinates(x = speciesCoord,
                    lon = "longitude",
                    lat = "latitude",
                    species = "species",
                    tests = c("capitals", "centroids", "equal", "gbif", "institutions", "zeros")
)        
                      
speciesCoord <- speciesCoord[flags$.summary,]
                      
#write to file ### LINK TO FIGSHARE ###
write.table(speciesCoord, 'species_coordinates.txt', quote=FALSE, row.names = FALSE, sep='\t')
