library(tidyverse)
library(MASS)
library(terra)
library(readxl)
library(geiger)
library(caper)

#read in ecoregion data (Dataset S7)
df <- read.delim('eco_data.txt')

#run regressions and generate summary statistics for every variable (Dataset S2)
sumstat <- data.frame(var=character(), p=double(), m=double())
for (i in 3:98) {
model <- glm.nb(df$richness ~ df[,i])
scaled_model <- lm(df$richness ~ scale(df[,i]))
p <- summary(model)$coefficients[2,4]
m <- summary(scaled_model)$coefficients[2,1]
r2 <- summary(scaled_model)$adj.r.squared
sumstat <- rbind(sumstat, data.frame(var=names(df)[i], p=p, m=m))
    }
sumstat$fdr <- p.adjust(sumstat$p, method = "fdr")

#filter insignificant variables
idx <- match(c('ecoregion', sumstat[sumstat$fdr <= 0.05,]$var), names(df))
df_sig <- df[,c(idx, 99)] 

#construct principal components
pca_out <- data.frame(ecoregion = df_sig$ecoregion, index = 1:length(df_sig$ecoregion))
pca_var <- data.frame(matrix(ncol=2,nrow=0))
colnames(pca_var) <- colnames(c('variable', 'variance_explained'))

formulas <- c(~humidMax+humidMean+humidMin, ~phDeep+phShallow, 
              ~h2o10Shallow+h2o10Deep+h2o33Deep+h2o33Shallow, ~windMax+windMean+windMean, ~clayDeep+clayShallow, ~siltDeep+siltShallow+sandDeep+sandShallow,
              ~h2o15kDeep+h2o15kShallow, ~sine+cosine, ~coarseEarthShallow+coarseEarthDeep)
              
              
names <- c('humid.pca', 'ph.pca', 'h2o.pca', 'wind.pca', 'clay.pca', 'sand.pca', 'h2o15k.pca', 'angle.pca', 'coarseEarth.pca')

for (i in 1:9) {
    pca <- prcomp(formulas[[i]], df_sig[,2:71], scale=FALSE)
    pca_out <- merge(pca_out, pca$x[,1], by.x='index', by.y='row.names', all=TRUE)
    names(pca_out)[length(names(pca_out))] <- names[i]
    
    pca_var <- rbind(pca_var, data.frame(name=names[i], variance_explained=pca$sdev[1]^2/sum(pca$sdev^2)))
    }
    
scaled_formulas <- c(~prec+precMax+moist+moistMax, ~moistMin+precMin, ~grow+NPP+soilco2, ~tempRange+tempStd+sunRange+iso+evapMin+sunMin+tempMin, ~rough+rugged+slope, ~Cdeep+Cshallow+Ndeep+Nshallow+OrgoCdeep+OrgoCshallow+fineEarthDeep+fineEarthShallow,
             ~sun+snow, ~evapMax+tempMax, ~forestBiomass+biomassAbove+biomassBelow)
              
scaled_names <- c('wet.pca', 'wetMin.pca', 'productivity.pca', 'tempRange.pca', 'topo.pca', 'soilRichness.pca', 'temp.pca', 'tempMax.pca', 'biomass.pca')

for (i in 1:9) {
    pca <- prcomp(scaled_formulas[[i]], df_sig[,2:82], scale=TRUE)
    pca_out <- merge(pca_out, pca$x[,1], by.x='index', by.y='row.names', all=TRUE)
    names(pca_out)[length(names(pca_out))] <- scaled_names[i]
    
    pca_var <- rbind(pca_var, data.frame(name=scaled_names[i], variance_explained=pca$sdev[1]^2/sum(pca$sdev^2)))
    }

eco_pca <- df_sig %>% dplyr::select(area, cultivated, dry, evapRange, EVIdiv, forest, forestIntegrity, geoClassDiv,
                             humidRange, impact, moistRange, plantRich, precVar, sunMax, tempDayRange, temperate, EVIvar,
                             tempStdLGM, windRange, precVarLGM, precMinLGM)

df_pca <- cbind(pca_out %>% dplyr::select(-index), eco_pca)
df_pca <- df_pca %>% mutate_at(.vars = c("wet.pca", "clay.pca", "sand.pca", "h2o15k.pca", "productivity.pca", "humid.pca", "wetMin.pca", "soilRichness.pca", "temp.pca", "biomass.pca"), function(x) {return(-x)})
df_pca$richness <- df_sig$richness

#run regressions and generate summary statistics for significant variables and principal components (Dataset S3)
sumstat_pca <- data.frame(var=character(), p=double(), coef=double())

for (i in 2:40) {
    model <- glm.nb(df_pca$richness ~ df_pca[,i])
    scaled_model <- lm(df_pca$richness ~ scale(df_pca[,i]))

    p <- summary(model)$coefficients[2,4]
    coef <- summary(scaled_model)$coefficients[2,1]
    sumstat_pca <- rbind(sumstat_pca, data.frame(var=names(df_pca)[i], p=p, coef=coef))
    }
sumstat_pca$fdr <- p.adjust(sumstat_pca$p, method = "fdr")

#relative importance analysis for highest performing variables
L <- c('forest', 'productivity.pca', 'h2o15k.pca', 'dry', 'ph.pca', 'wet.pca', 'geoClassDiv', 'biomass.pca', 'plantRich', 'EVIdiv', 'wetMin.pca', 'temperate', 'topo.pca', 'forestIntegrity', 'sunMax', 'humid.pca')
res <- Map(combn, list(L), seq_along(L), simplify = FALSE)
res <- unlist(res, recursive = FALSE)

for (R in res) {
aic <- AIC(glm.nb(`df_pca$richness` ~ ., data=cbind(df_pca$richness, data.frame(df_pca[,unlist(R)]))))
write(paste(c(paste(unlist(R), collapse='+')), aic, sep='\t'), "tmpaic.txt", append=TRUE)
}

aic <- read.delim('aic.txt', header=FALSE, col.names = c('var', 'AIC'))

aic$delta_AIC <- aic$AIC - min(aic$AIC)
aic$relative_likelihood <- exp(-0.5 * aic$delta_AIC)
aic$Akaike_weight <- aic$relative_likelihood / sum(aic$relative_likelihood)

#make sure forest variable is renamed "forest1" to avoid overlap with "forestIntegrity"
outdf <- data.frame(var=character(), relimpo=numeric())
L <- c('forest1', 'productivity.pca', 'h2o15k.pca', 'dry', 'ph.pca', 'wet.pca', 'geoClassDiv', 'biomass.pca', 'plantRich', 'EVIdiv', 'wetMin.pca', 'temperate', 'topo.pca', 'forestIntegrity', 'sunMax', 'humid.pca')
for (var in L) {
    outdf <- rbind(outdf, data.frame(var=var, relimpo=sum(aic[grepl(var, aic$var),]$Akaike_weight)))
    }

#generalist specialist classifications from https://doi.org/10.1101/2023.06.19.545611
breadth <- read.delim('Generalist_Specialist_Table_names_10_22.txt')

#range size, absolute latitude, and overlapping species richness for each species (Dataset S10)
occupancy <- read.delim('species_lat_rich.txt')

#load tree from https://doi.org/10.1101/2023.06.19.545611 and reconcile taxonomies 
tree <- read.tree('1154yeasts_1403OGs_ml_timetree_newnames.tree')
taxonomy <- read_excel('Draft_TableS1_20230201.xlsx')

breadth$species <- taxonomy$NEW_Tip_ID[match(unlist(breadth$assembly_fullID_updated), taxonomy$assembly_fullID_updated)]
occupancy$species <- taxonomy$NEW_Tip_ID[match(unlist(occupancy$species), taxonomy$Species)]
occupancy <- na.omit(occupancy)

yeast_compare <- comparative.data(phy = tree, data = occupancy, names.col = species, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

#PGLS for latitude and range size
latitude.pgls<-pgls(range_size ~ absolute_latitude, data = yeast_compare, lambda='ML')
summary(latitude.pgls)

#PGLS for co-occurring diversity and range size
diversity.pgls<-pgls(range_size ~ richness, data = yeast_compare, lambda='ML')
summary(diversity.pgls)

breadth$Carb_Class <- ifelse(breadth$Carb_Class=='Specialist', 'Specialist', 'Other')
df <- merge(breadth, occupancy)

#phylogenetic ANOVA for specialist class and range size
anova.df <- df %>% dplyr::select(c(species, Carb_Class, occupancy)) %>% na.omit()
dat <- anova.df[3]
rownames(dat) <- unlist(anova.df$species)
group <- factor(anova.df$Carb_Class)
names(group) <- unlist(anova.df$species)
phy.manova <- aov.phylo(dat~group, tree, nsim = 1000)

summary(phy.manova)
