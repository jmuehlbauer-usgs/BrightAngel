##### Bright Angel trophic cascade analysis #####
## Last updated 3 October 2019 by J.D. Muehlbauer


##### Set up workspace ##### 

## Load/install requisite packages
source('https://github.com/jmuehlbauer-usgs/R-packages/blob/master/packload.r?raw=TRUE')
packload(c('devtools', 'lubridate', 'plots', 'bugR', 'glmmTMB'))

## Set working directory (depends on whether Jeff or Megan)
if(Sys.info()[6]=='jmuehlbauer'){
	setwd('C:/Users/jmuehlbauer/Documents/Projects/BrightAngel')
} else{
	setwd('C:/Users/mdaubert/Documents/Bright Angel/R')
}

## Create folder for putting figures in
dir.create('Figures', showWarnings = FALSE)


##### Read in data #####

## Read and write drift and benthic data from GCMRC
	## Note: Only need to run this function once. Then only if the data need updating for some reason
BAreadwrite <- function(){
	## Get all raw data from database
	dbdat <- 'https://raw.githubusercontent.com/jmuehlbauer-usgs/Database/master/'
	files <- c('DriftSample', 'DriftSpecimen', 'tbl_BenthicSample', 
		'tbl_BenthicSpecimen', 'SpeciesList')
	gcmrc <- lapply(paste0(dbdat, files, '.csv'), read.csv)
		names(gcmrc) <- c('dsamp', 'dspec', 'bsamp', 'bspec', 'sppl')
	## Format dates and times for use
	gcmrc1 <- gcmrc
		gcmrc1$dsamp$Date <- as.Date(gcmrc$dsamp$Date, format = '%m/%d/%Y')
		gcmrc1$bsamp$Date <- as.Date(gcmrc$bsamp$SampleDate, format = '%m/%d/%Y')
		gcmrc1$bsamp$ProcessDate <- as.Date(gcmrc$bsamp$DateProcessed, format = '%m/%d/%Y')
	## Add barcodes to benthic specimens
	gcmrc1$bspec$BarcodeID <- gcmrc1$bsamp[match(gcmrc1$bspec$SampleID, gcmrc1$bsamp$SampleID), 'BarcodeID']
	## Subset to just the samples of interest (Bright Angel in 2016-January 2017)
	gcmrc2 <- gcmrc1
		gcmrc2$dsamp <- gcmrc1$dsamp[gcmrc1$dsamp$Reach == 'BrightAngel' & (year(gcmrc1$dsamp$Date) == 
			2016 | (year(gcmrc1$dsamp$Date) == 2017 & month(gcmrc1$dsamp$Date) == 1)),]
			gcmrc2$dsamp <- droplevels(gcmrc2$dsamp)
		gcmrc2$bsamp <- gcmrc1$bsamp[gcmrc1$bsamp$River == 'Bright Angel' & (year(gcmrc1$bsamp$Date) == 
			2016 | (year(gcmrc1$bsamp$Date) == 2017 & month(gcmrc1$bsamp$Date) == 1)),]
			gcmrc2$bsamp <- droplevels(gcmrc2$bsamp)
	## Format benthic sample data to parallel drift, rename or drop columns
	gcmrc2$bsamp$Region <- 'GrandCanyon'
	gcmrc2$bsamp$Reach <- 'BrightAngel'
	gcmrc2$bsamp$SampleNumber <- gcmrc2$bsamp$DatasheetSampleNo
	gcmrc2$bsamp$Depth <- gcmrc2$bsamp$SampleDepth
	gcmrc3 <- gcmrc2
	gcmrc3$bsamp <- gcmrc2$bsamp[,c('BarcodeID', 'TripID', 'Region', 'Reach', 'Date', 'SampleNumber', 
		'RiverMile', 'Depth', 'GearID', 'SampleArea', 'EntererSample', 'Processor', 
		'ProcessDate', 'ProcessTime', 'EntererSpecimen', 'Notes')]
		gcmrc3$bsamp$TripID <- ifelse(month(gcmrc3$bsamp$Date) == 6, 'BA20160608',
			ifelse(month(gcmrc3$bsamp$Date) == 11, 'BA20161108', 
			ifelse(month(gcmrc3$bsamp$Date) == 1, 'BA20170118', 'BA20160831')))
		gcmrc3$bsamp$TripID <- as.factor(gcmrc3$bsamp$TripID)
	## Format benthic specimen data to include barcodes and total counts, rename or drop columns
	gcmrc3$bspec$CountC <- rowSums(gcmrc3$bspec[,6:37], na.rm = TRUE)
	gcmrc3$bspec$CountF <- rowSums(gcmrc3$bspec[,40:56], na.rm = TRUE)
	gcmrc3$bspec$CountTotal <- rowSums(gcmrc3$bspec[,c(6:37, 40:56)], na.rm = TRUE)
	colnames(gcmrc3$bspec)[which(colnames(gcmrc3$bspec)=='Cpt5')] <- 'C0'
	colnames(gcmrc3$bspec)[which(colnames(gcmrc3$bspec)=='Fpt5')] <- 'F0'
	ccols <- paste0('C', 0:30)
	fcols <- paste0('F', 0:15)
	gcmrc4 <- gcmrc3
	gcmrc4$bspec <- gcmrc3$bspec[, c('BarcodeID', 'SpeciesID', ccols, 'CExtra', fcols, 'FExtra', 
		'CountC', 'CountF', 'CountTotal', 'Notes')]
	## Convert NA counts to 0s
	gcmrc4$bspec[is.na(gcmrc4$bspec)] <- 0
		gcmrc4$bspec <- droplevels(gcmrc4$bspec)	
	## Merge benthic sample and specimen data together. Same for drift
		d1 <- merge(gcmrc4$dsamp, gcmrc4$dspec, by = 'BarcodeID')
		b1 <- merge(gcmrc4$bsamp, gcmrc4$bspec, by = 'BarcodeID')	
	## Get densities and concentrations
	d1$Concentration <- d1$CountTotal / d1$Volume
	b1$Density <- b1$CountTotal / b1$SampleArea
	## Write data
	write.csv(d1, 'Data/DriftData.csv', row.names = FALSE)
	write.csv(b1, 'Data/BenthicData.csv', row.names = FALSE)
	write.csv(gcmrc4$sppl, 'Data/SpeciesList.csv', row.names = FALSE)

}
#BAreadwrite()

## Get data from GitHub
gitdat <- 'https://raw.githubusercontent.com/jmuehlbauer-usgs/BrightAngel/master/Data/'
gitfiles <- c('DriftData', 'BenthicData', 'WhitingData')
dat <- lapply(paste0(gitdat, gitfiles, '.csv'), read.csv, colClasses = c('Date' = 'character'))
		names(dat) <- c('Drift', 'Benthic', 'Whiting')
spp <- read.csv('https://raw.githubusercontent.com/jmuehlbauer-usgs/Database/master/SpeciesList.csv')


##### Clean up data #####

## Convert dates usable format
dat <- lapply(dat, function(x) {
	if(length(grep('/', x[,'Date'][1])) > 0){x[,'Date'] <- as.Date(x[,'Date'], format = '%m/%d/%Y')
	} else{x[,'Date'] <- as.Date(x[,'Date'])}
	return(x)
	})

## Give Whiting Data a barcode name for consistency
dat$Whiting$BarcodeID <- paste(dat$Whiting$Date, dat$Whiting$SampleID)

## Add functional feeding groups
dat <- lapply(dat, transform, FFG = spp[match(SpeciesID, spp$SpeciesID), 'FFG'])

## Add season
dat <- lapply(dat, transform, Season = factor(ifelse(month(Date) == 11, 'November', 
	ifelse(month(Date) == 1, 'January', ifelse(month(Date) == 6, 'June', 'September'))),
	levels = c('November', 'January', 'June', 'September')))

## Add sample area, raw counts, biomass, and mean size to Whiting data
	## Note: Whiting used a 0.086 m^2 Hess, but combined 2 samples per each of his 6 "samples"
dat$Whiting$Area <- 0.086 * 2
dat$Whiting$CountTotal <- round(dat$Whiting$Density * dat$Whiting$Area)
dat$Whiting$Biomass <- round(dat$Whiting$Biomass * dat$Whiting$Area, 2)
Wmatch <- match(dat$Whiting$SpeciesID, spp$SpeciesID)
dat$Whiting$Size <- round(exp(log(dat$Whiting$Biomass / dat$Whiting$CountTotal / 
	spp[Wmatch, 'RegressionA']) / spp[Wmatch, 'RegressionB']), 2)

## Add biomass and mean size for GCMRC data
for(i in 1 : 2){
	t1 <- dat[[i]]
	sizes <- if('C30' %in% colnames(t1)){1:30} else{1:20}
	reps <- c(0.5, 1:15, 0.5, sizes)
	mycols <- c(paste0('F', 0:15), 'C0', paste0('C', sizes))
	t2 <- t1[, mycols]
	lsize <- matrix(reps, ncol = length(reps), nrow = dim(t2)[1], byrow = TRUE)
	ABs <- spp[match(t1$SpeciesID, spp$SpeciesID), c('RegressionA', 'RegressionB')]
	t3 <- round(t2 * (lsize^ABs$RegressionB) * ABs$RegressionA, 5)
	dat[[i]]$Biomass <- round((rowSums(t3[, 1: 16]) * 
		ifelse(t1$CountF > 0, t1$CountF / (t1$CountF - t1$FExtra), 0)) + 
		(rowSums(t3[, 17 : (17 + length(sizes + 1))]) * 
		ifelse(t1$CountC > 0, t1$CountC / (t1$CountC - t1$CExtra), 0)), 2)
	dat[[i]]$Size <- round(rowSums(t2 * lsize) / (t1$CountTotal - (t1$CExtra + t1$FExtra)), 2)
}

## Rename Benthic Area column for consistency
colnames(dat$Benthic)[which(colnames(dat$Benthic) == 'SampleArea')] <- 'Area'

## Keep only columns of interest
dat0 <- lapply(dat, function(x){
	mycols <- c('BarcodeID', 'Date', 'Season', 
		ifelse('RiverMile' %in% colnames(x), 'RiverMile', 'BarcodeID'), 
		ifelse('Density' %in% colnames(x), 'Area', 'Volume'), 
		'SpeciesID', 'FFG', 'CountTotal', 'Size', 'Biomass')
	x[, mycols]
	})
dat0 <- lapply(dat0, setNames, nm = c('BarcodeID', 'Date', 'Season', 'RiverMile', 'Unit', 'SpeciesID', 
	'FFG', 'Count', 'Size', 'Biomass'))
	dat0$Whiting$RiverMile <- 1000
	
	
##### Clean up taxa #####

## Get taxa list of present taxa
taxa <- rbind(spp[spp$SpeciesID %in% dat0$Drift$SpeciesID,c('SpeciesID', 'Description')], 
	spp[spp$SpeciesID %in% dat0$Benthic$SpeciesID, c('SpeciesID', 'Description')], 
	spp[spp$SpeciesID %in% dat0$Whiting$SpeciesID,c('SpeciesID', 'Description')])
	taxa <- taxa[match(unique(taxa$SpeciesID), taxa$SpeciesID),]
	taxa <- taxa[order(taxa$Description),]

## Convert all taxa in different life stages to same Species ID (e.g., CHIL, CHIP, CHIA all become CHIA)
dat1 <- lapply(dat0, transform, SpeciesID = as.character(SpeciesID))
origt <- c('MCYA', 'CERA', 'CERP', 'CHIA', 'CHIP', 'CULP', 'EMPA', 'WIEA', 'SIMA', 'SIMP', 'BASP', 'BAET', 
	'LEPA', 'CAPA', 'TRIA', 'TRIP', 'HYSP', 'HYDA')
newt <- c('MCYL', 'CERL', 'CERL', 'CHIL', 'CHIL', 'CULL', 'EMPL', 'WIEL', 'SIML', 'SIML', 'BAEL', 'BAEL', 
	'LEPL', 'CAPL', 'TRIL', 'TRIL', 'HYDE', 'HYDL')
dat1 <- lapply(dat1, transform, SpeciesID = ifelse(SpeciesID %in% origt, 
	newt[match(SpeciesID, origt)], SpeciesID))

## Build table to look for congenerics to combine (only Whiting and Benthic)
	## Note: Due to possible errors/discrepancies between our data and Whiting's that aren't real.
unq <- unique(unlist(lapply(dat1, function(x) unique(x[,'SpeciesID']))))
cnt <- lapply(dat1, function(x) tapply(x[, 'Count'], x[, 'SpeciesID'], sum))
spptab <- data.frame(SpeciesID = unq, Drift = NA, Benthic = NA, Whiting = NA)
	spptab$Description <- spp[match(spptab$SpeciesID, spp$SpeciesID), 'Description']
for(i in 1:3){
	spptab[, i + 1][match(names(cnt[[i]]), spptab$SpeciesID)] <- cnt[[i]]
}
spptabW <- spptab[order(-spptab$Benthic, -spptab$Benthic, spptab$SpeciesID),]
spptabB <- spptab[order(-spptab$Whiting, -spptab$Benthic, spptab$SpeciesID),]
spptab1 <- spptab[order(spptab$Description),]
	## Note: Whiting saw Chloroperlidae, which we didn't, but there is no obvious substitution.
	## Cutting Nematoda, Lymnaeidae, Veliidae, Veneroida, which Whiting doesn't have 
		## (probably didn't count them).
	## Whiting probably didn't separate Empidid or Tabanid taxa. He also counted Ostracods (we didn't).

## Combine congenerics
origt1 <- c('ELML', 'ELOA', 'PROB', 'HEMR', 'WIEL', 'SILV', 'TABS', 'DICL', 'DRAL', 'DAML', 'COEN', 
	'CAPL', 'TRIL', 'HYOS', 'HYDL', 'HYLA', 'POLY', 'RHCL')
newt1 <- c('MCYL', 'MCYL', 'CERL', 'EMPL', 'EMPL', 'TABL', 'TABL', 'TIPL', 'LIBE', 'ARGI', 'ARGI', 
	'CAPN', 'HYDE', 'HYDE', 'LETR', 'LETR', 'POLL', 'RHYL')
dat2 <- lapply(dat1, transform, SpeciesID = ifelse(SpeciesID %in% origt1, 
	newt1[match(SpeciesID, origt1)], SpeciesID))

## Cut oddballs, delete 0 count rows
dat2 <- lapply(dat2, function(x) x[!(x[, 'SpeciesID'] %in% c('NEMA', 'LYMN', 'RHAG', 'CLAM', 'OSTR')),])
dat2 <- lapply(dat2, function(x) x[x[, 'Count'] != 0,])

## Combine rows of same taxa from different life stages
combinefx <- function(oldlist, groupby){
	tlist <- list()
	for(i in 1:length(oldlist)){
		t1 <- oldlist[[i]]
			t1$BarGrp <- paste(t1$BarcodeID, t1[,groupby])
		tCnt <- aggregate(t1$Count, by = list(t1$BarGrp), sum, na.rm = TRUE)
		tSize <- aggregate(t1$Count * t1$Size, by = list(t1$BarGrp), sum, na.rm = TRUE)
			tSize[, 2] <- round(tSize[, 2] / tCnt[, 2], 2)
		tBio <- aggregate(t1$Biomass, by = list(t1$BarGrp), sum, na.rm = TRUE)
		t2 <- t1[match(unique(t1$BarGrp), t1$BarGrp),]
			t2$Count <- tCnt[match(t2$BarGrp, tCnt[, 1]), 2]
			t2$Size <- tSize[match(t2$BarGrp, tSize[, 1]), 2]
			t2$Biomass <- tBio[match(t2$BarGrp, tBio[, 1]), 2]
		tlist[[i]] <- subset(t2, select = -BarGrp)
	}
	names(tlist) <- names(oldlist)
	return(tlist)
}
dat3 <- combinefx(dat2, 'SpeciesID')
	dat3 <- lapply(dat3, transform, FFG = spp[match(SpeciesID, spp$SpeciesID), 'FFG'])

## Limit to only aquatic taxa
dat4 <- lapply(dat3, function(x){x[x[, 'SpeciesID'] %in% spp[spp$Habitat == 'Aquatic', 'SpeciesID'],]})

## Reassign FFGs in case any got messed up by the SpeciesID combining
	## Combine shredders into Collector-Gatherers--
		## None found for our benthics, so a problem modeling otherwise.
dat4 <- lapply(dat4, function(x){
	x[,'FFG'] <- spp[match(x[, 'SpeciesID'], spp$SpeciesID), 'FFG']
	levels(x[, 'FFG'])[which(levels(x[, 'FFG']) == 'Shredder')] <- 'CollectorGatherer'
	return(x)
	})

## Convert SpeciesIDs back to factors
dat4 <- lapply(dat4, transform, SpeciesID = as.factor(SpeciesID))

## Sort by Date, BarcodeID, SpeciesID, drop old factor levels
dat4 <- lapply(dat4, function(x) x[order(x[, 'Date'], x[, 'BarcodeID'], x[, 'SpeciesID']),])
	dat4 <- lapply(dat4, droplevels)

## Keep only 1000 m benthic site
dat5 <- dat4
dat5$Benthic <- dat4$Benthic[dat4$Benthic$RiverMile == 1000,]
	dat5$Benthic <- droplevels(dat5$Benthic)
dat5 <- lapply(dat5, subset, select = -RiverMile)

## Add consistent factor levels
levspp <- sort(unique(unlist(lapply(dat5, function(x) levels(x[, 'SpeciesID'])))))
levFFG <- unique(unlist(lapply(dat5, function(x) levels(x[, 'FFG']))))
	levFFG <- ifelse(length(levFFG) == 6, c('ScraperGrazer', 'CollectorFilterer', 
		'CollectorGatherer', 'Generalist', 'Predator'), levFFG)
dat5 <- lapply(dat5, function(x){
	levspp1 <- levspp[!(levspp %in% levels(x[, 'SpeciesID']))]
	levFFG1 <- levFFG[!(levFFG %in% levels(x[, 'FFG']))]
	levels(x[, 'SpeciesID']) <- c(levels(x[, 'SpeciesID']), levspp[!(levspp %in% levels(x[, 'SpeciesID']))])
	levels(x[, 'FFG']) <- c(levels(x[, 'FFG']), levFFG[!(levFFG %in% levels(x[, 'FFG']))])
	x[, 'SpeciesID'] <- factor(x[, 'SpeciesID'], levels = levspp)
	x[, 'FFG'] <- factor(x[, 'FFG'], levels = c('ScraperGrazer', 'CollectorFilterer',
		'CollectorGatherer', 'Generalist', 'Predator'))
	return(x)
	})


##### Group specimens by Sample #####

## Create new list, with samples summed.
samp <- combinefx(dat5, 'BarcodeID')
	samp <- lapply(samp, subset, select = -SpeciesID)

## Combine benthic and Whiting into single dataframe for analysis
samp1 <- rbind(samp$Whiting, samp$Benthic)
	samp1$Study <- factor(ifelse(year(samp1$Date) < 2015, 'Pre', 'Post'), levels = c('Pre', 'Post'))
	samp1 <- subset(samp1, select = -FFG)
	
## Build some models for count data
sampmod <- list()
	sampmod[[1]] <- sampmod[[2]] <- sampmod[[3]] <- list()
	names(sampmod) <- c('Count', 'Size', 'Biomass')
	
## Find best data distribution (family) to use
sampmod[[1]]$lm1 <- glmmTMB(Count ~ 1 + offset(Unit), data = samp1)
sampmod[[1]]$po1 <- glmmTMB(Count ~ 1 + offset(log(Unit)), family = 'poisson', data = samp1)
sampmod[[1]]$nb1 <- glmmTMB(Count ~ 1 + offset(log(Unit)), family = 'nbinom2', data = samp1)

## Add a random effect for sample or season
sampmod[[1]]$nb2 <- glmmTMB(Count ~ 1 + (1 | BarcodeID) + offset(log(Unit)), 
	family = 'nbinom2', data = samp1)
sampmod[[1]]$nb3 <- glmmTMB(Count ~ 1 + (1 | Season) + offset(log(Unit)), 
	family = 'nbinom2', data = samp1)
sampmod[[1]]$nb4 <- glmmTMB(Count ~ 1 + (1 | BarcodeID) + (1 | Season) + offset(log(Unit)), 
	family = 'nbinom2', data = samp1)
	
## Add pre-post fixed effect
sampmod[[1]]$nb5 <- suppressWarnings(glmmTMB(Count ~ 0 + Study + (1 | BarcodeID) + offset(log(Unit)), 
	family = 'nbinom2', data = samp1))
sampmod[[1]]$nb6 <- glmmTMB(Count ~ 0 + Study + (1 | Season) + offset(log(Unit)), 
	family = 'nbinom2', data = samp1)
sampmod[[1]]$nb7 <- suppressWarnings(glmmTMB(Count ~ 0 + Study + (1 | BarcodeID) + (1 | Season) + 
	offset(log(Unit)), family = 'nbinom2', data = samp1))
	## Note: nb7 doesn't fit well and isn't a good model, hence the warning suppression.

## Compare models using AIC
AICtable <- function(modlist){
	lapply(modlist, function(x) data.frame(
		AIC = sapply(x, function(y) round(AIC(y), 2)),
		FixedEffect = sapply(x, function(y){
			att1 <- attributes(y$modelInfo$reTrms$cond$terms$fixed)$term.labels
			att2 <- ifelse(length(att1) > 0, att1, 1)
			return(att2)
			}),
		RandomEffect = sapply(x, function(y) {
			pas1 <- paste(y$modelInfo$grpVar, collapse = ', ')
			pas2 <- ifelse(pas1 == '', 'None', pas1)
			return(pas2)
			})
		))
	}
sampmodAIC <- AICtable(sampmod)
	## Negative binomial distribution is best bet (over linear/Gaussian and Poisson).
	## A Season random effect is co-equivalent with no random effect.
		## Having this random effect and makes logical sense though given uneven sample effort by season.
	## Including a pre-post term improves the Count (density) model.

## Build models for size and biomass data
	## Keeping random effect as-is from count model.
sampmod[[2]]$lm6 <- glmmTMB(Size ~ 0 + Study + (1 | Season),
	family = 'gaussian', data = samp1)
sampmod[[2]]$ln6 <- glmmTMB(log(Size) ~ 0 + Study + (1 | Season),
	family = 'gaussian', data = samp1)
sampmod[[2]]$ga6 <- glmmTMB(Size ~ 0 + Study + (1 | Season),
	family = Gamma(link = 'log'), data = samp1)
sampmod[[3]]$lm6 <- glmmTMB(Biomass ~ 0 + Study + (1 | Season) + offset(Unit), 
	family = 'gaussian', data = samp1)
sampmod[[3]]$ln6 <- glmmTMB(log(Biomass) ~ 0 + Study + (1 | Season) + offset(log(Unit)), 
	family = 'gaussian', data = samp1)
sampmod[[3]]$ga6 <- glmmTMB(Biomass ~ 0 + Study + (1 | Season) + offset(log(Unit)), 
	family = Gamma(link = 'log'), data = samp1)
sampmodAIC <- AICtable(sampmod)
	for(i in 1 : length(sampmodAIC)){
		sampmodAIC[[i]][,c('ModelPre', 'ModelPost')] <- t(sapply(sampmod[[i]], function(x){
			if(grepl('Study', x$call[2])){
				if(x$modelInfo$family$link == 'log' | grepl('log', x$call[2]) == TRUE){
					round(exp(fixef(x)[[1]]), 2)
				} else{round(fixef(x)[[1]], 2)}
			} else{rep(NA, 2)}
		}))
	}
sampmodAIC$SampMeans <- data.frame(Metric = c('Count', 'Size', 'Biomass'),
	Pre = rep(NA, length(sampmodAIC)), Post = rep(NA, length(sampmodAIC)))
	for(i in 1 : dim(sampmodAIC$SampMeans)[1]){
		if(sampmodAIC$SampMeans$Metric[i] == 'Size'){
			sampmodAIC$SampMeans[i, 2:3] <- round(tapply(samp1[,
				as.character(sampmodAIC$SampMeans$Metric[i])], 
				samp1$Study, mean, na.rm = TRUE), 2)
		} else{
			sampmodAIC$SampMeans[i, 2:3] <- round(tapply(samp1[,
				as.character(sampmodAIC$SampMeans$Metric[i])] / samp1$Unit, 
				samp1$Study, mean, na.rm = TRUE), 2)
		}
	}
	## Looking at AICs and model estimates vs. sample means, use nb6 for Count, lm6 for Size, ga6 for Biom.
		## "6" indicates a Season random effect, nb is negative binomial, lm is gaussian, ga is gamma.
		
## Get best model for Count, Size, and Biomass
bestmod <- function(mods){
	lapply(mods, function(x) {
	x[names(which(sapply(x, AIC) == min(sapply(x, AIC), na.rm = TRUE)))][[1]]
	})
}
	## Note: function works later, but not here (due to lognormal).
bestmod0 <- with(sampmod, list(Count$nb6, Size$lm6, Biomass$ga6))
	names(bestmod0) <- sampmodAIC$SampMeans$Metric


##### Get predicted values and plot change #####

## Get parameters to predict
predparm0 <- data.frame(Study = levels(samp1$Study), Season = 'Generic', Unit = 1)
	
## Get model-predicted values and confidence intervals
preds0l <- lapply(bestmod0, predict, newdata = predparm0, type = 'link', se.fit = TRUE, 
	allow.new.levels = TRUE)
preds0r <- lapply(bestmod0, predict, newdata = predparm0, type = 'response', se.fit = TRUE, 
	allow.new.levels = TRUE)
fits0 <- list()
for(i in 1 : length(bestmod0)){
	fits0[[i]] <- data.frame(Study = predparm0$Study)
	fits0[[i]]$LinkFit <- preds0l[[i]]$fit
	fits0[[i]]$LinkSE <- preds0l[[i]]$se.fit
	fits0[[i]]$Fit <- preds0r[[i]]$fit
	fits0[[i]]$SELower <- preds0r[[i]]$fit - preds0r[[i]]$se.fit
	fits0[[i]]$SEUpper <- preds0r[[i]]$fit + preds0r[[i]]$se.fit
	if(bestmod0[[i]][6]$modelInfo$family$link == 'log'){
		fits0[[i]]$CILower95 <- exp(fits0[[i]]$LinkFit - (qnorm(0.975) * fits0[[i]]$LinkSE))
		fits0[[i]]$CIUpper95 <- exp(fits0[[i]]$LinkFit + (qnorm(0.975) * fits0[[i]]$LinkSE))
		fits0[[i]]$CILower50 <- exp(fits0[[i]]$LinkFit - (qnorm(0.75) * fits0[[i]]$LinkSE))
		fits0[[i]]$CIUpper50 <- exp(fits0[[i]]$LinkFit + (qnorm(0.75) * fits0[[i]]$LinkSE))			
	} else{
		fits0[[i]]$CILower95 <- fits0[[i]]$LinkFit - (qnorm(0.975) * fits0[[i]]$LinkSE)
		fits0[[i]]$CIUpper95 <- fits0[[i]]$LinkFit + (qnorm(0.975) * fits0[[i]]$LinkSE)
		fits0[[i]]$CILower50 <- fits0[[i]]$LinkFit - (qnorm(0.75) * fits0[[i]]$LinkSE)
		fits0[[i]]$CIUpper50 <- fits0[[i]]$LinkFit + (qnorm(0.75) * fits0[[i]]$LinkSE)		
	}
}	
	fits0 <- lapply(fits0, function(x){cbind(Study = x[, 1], round(x[, -1], 2))})
	names(fits0) <- names(bestmod0)

## Plot panel graph of overall change in Density, Size, and Biomass
plotpred <- function(){
par(mfrow = c(3, 1), mar = c(2, 5, 0.1, 0.1), cex = 1)	
for(i in 1:3){
		mydat <- fits0[[i]]
		plot(c(0.5, 2.5), c(min(mydat$CILower95), max(mydat$CIUpper95)), xlab = '', ylab = '', 
			axes = FALSE, type = 'n')
		if(i != 3){axis(1, at = 1:2, labels = FALSE)
		} else{axis(1, at = 1:2, labels = c('2011', '2016'))
		}
		axis(2, las = 2)
		if(i == 1){mtext(side = 2, bquote('Invertebrate density (#*'*m^-2*')'), line = 3.5)}
		if(i == 2){mtext(side = 2, 'Mean invertebrate size (mm)', line = 3.5)}
		if(i == 3){mtext(side = 2, bquote('Invertebrate biomass (mg*'*m^-2*')'), line = 3.5)}
		box(bty = 'l')
		with(mydat, arrows(x0 = 1:2, y0 = CILower95, y1 = CIUpper95, code = 3, angle = 90, length = 0.05))
		with(mydat, arrows(x0 = 1:2, y0 = CILower50, y1 = CIUpper50, length = 0, col = 'grey50', lwd = 4))
		points(mydat$Fit, pch = c(17, 16), col = c(2, 4), cex = 1.5)
		legend('topright', LETTERS[i], bty = 'n')
	}
}
plotTypes(plotpred, 'ModelDensitySizeBiomassOverall', 'Figures', width = 4)


##### Group specimens by FFG #####

## Create new list, with samples summed by FFG
ffg <- combinefx(dat5, 'FFG')
	ffg <- lapply(ffg, subset, select = -SpeciesID)

## Get relative counts and biomasses
ffg <- lapply(ffg, function(x){
	cnttot <- tapply(x[, 'Count'], x[, 'BarcodeID'], sum)
	biotot <- tapply(x[, 'Biomass'], x[, 'BarcodeID'], sum)
	x[, 'RelCount'] <- round(x[, 'Count'] / cnttot[match(x[, 'BarcodeID'], names(cnttot))], 4)
	x[, 'RelBiomass'] <- round(x[, 'Biomass'] / biotot[match(x[, 'BarcodeID'], names(biotot))], 4)
	return(x)
	})

## Combine benthic and Whiting into single dataframe for analysis
ffg1 <- rbind(ffg$Whiting, ffg$Benthic)
	ffg1$Study <- factor(ifelse(year(ffg1$Date) < 2015, 'Pre', 'Post'), levels = c('Pre', 'Post'))
	ffg1 <- ffg1[, c('BarcodeID', 'Date', 'Study', 'Season', 'Unit', 'FFG', 'Count', 'Size', 'Biomass', 
		'RelCount', 'RelBiomass')]

	
##### Build model for FFG change #####

## Build lists
ffgmod <- list()
	ffgmod[[1]] <- ffgmod[[2]] <- ffgmod[[3]] <- list()
	names(ffgmod) <- names(sampmod)

## Start model selection using basic best model structure from sampmod (above)
ffgmod[[1]]$nb1 <- glmmTMB(Count ~ 0 + Study + (1 | Season) + offset(log(Unit)), 
	family = 'nbinom2', data = ffg1)
ffgmod[[1]]$nb2 <- glmmTMB(Count ~ 0 + Study * FFG + (1 | Season) + offset(log(Unit)), 
	family = 'nbinom2', data = ffg1)
ffgmod[[2]]$lm1 <- glmmTMB(Size ~ 0 + Study + (1 | Season), family = 'gaussian', data = ffg1)
ffgmod[[2]]$lm2 <- glmmTMB(Size ~ 0 + Study * FFG + (1 | Season), family = 'gaussian', data = ffg1)
ffgmod[[3]]$ga1 <- glmmTMB(Biomass ~ 0 + Study + (1 | Season) + offset(log(Unit)), 
	family = Gamma(link = 'log'), data = ffg1)
ffgmod[[3]]$ga2 <- glmmTMB(Biomass ~ 0 + Study * FFG + (1 | Season) + offset(log(Unit)), 
	family = Gamma(link = 'log'), data = ffg1)

## Get best model for Count, Size, and Biomass
ffgmodAIC <- AICtable(ffgmod)
	## Adding FFG is a good idea for all three model types
bestmod1 <- bestmod(ffgmod)


##### Get predicted values and plot FFG change #####

## Get parameters to predict
predparm1 <- data.frame(FFG = rep(levels(ffg1$FFG), 2), Season = 'Generic', Unit = 1,
	Study = rep(levels(samp1$Study), c(length(levels(ffg1$FFG)), length(levels(ffg1$FFG)))))

## Get model-predicted values and confidence intervals
preds1l <- lapply(bestmod1, predict, newdata = predparm1, type = 'link', se.fit = TRUE, 
	allow.new.levels = TRUE)
preds1r <- lapply(bestmod1, predict, newdata = predparm1, type = 'response', se.fit = TRUE, 
	allow.new.levels = TRUE)
fits1 <- list()	
for(i in 1 : length(bestmod1)){
	fits1[[i]] <- data.frame(Study = predparm1$Study, FFG = predparm1$FFG)
	fits1[[i]]$LinkFit <- preds1l[[i]]$fit
	fits1[[i]]$LinkSE <- preds1l[[i]]$se.fit
	fits1[[i]]$Fit <- preds1r[[i]]$fit
	fits1[[i]]$SELower <- preds1r[[i]]$fit - preds1r[[i]]$se.fit
	fits1[[i]]$SEUpper <- preds1r[[i]]$fit + preds1r[[i]]$se.fit
	if(bestmod1[[i]][6]$modelInfo$family$link == 'log'){
		fits1[[i]]$CILower95 <- exp(fits1[[i]]$LinkFit - (qnorm(0.975) * fits1[[i]]$LinkSE))
		fits1[[i]]$CIUpper95 <- exp(fits1[[i]]$LinkFit + (qnorm(0.975) * fits1[[i]]$LinkSE))
		fits1[[i]]$CILower50 <- exp(fits1[[i]]$LinkFit - (qnorm(0.75) * fits1[[i]]$LinkSE))
		fits1[[i]]$CIUpper50 <- exp(fits1[[i]]$LinkFit + (qnorm(0.75) * fits1[[i]]$LinkSE))			
	} else{
		fits1[[i]]$CILower95 <- fits1[[i]]$LinkFit - (qnorm(0.975) * fits1[[i]]$LinkSE)
		fits1[[i]]$CIUpper95 <- fits1[[i]]$LinkFit + (qnorm(0.975) * fits1[[i]]$LinkSE)
		fits1[[i]]$CILower50 <- fits1[[i]]$LinkFit - (qnorm(0.75) * fits1[[i]]$LinkSE)
		fits1[[i]]$CIUpper50 <- fits1[[i]]$LinkFit + (qnorm(0.75) * fits1[[i]]$LinkSE)		
	}
}	
	fits1 <- lapply(fits1, function(x){cbind(Study = x[, 1], FFG = x[, 2], round(x[, c(-1, -2)], 2))})
	names(fits1) <- names(bestmod1)

## Plot combined graph of Whiting's and our data
plotmod <- function(){
	par(mfrow = c(length(fits1), 1), mar = c(1.4, 4.5, 0.1, 0.7), oma = c(1.5, 0, 0, 0), cex = 1)
	d1 <- 0.09
	for(i in 1: length(fits1)){
		mydat <- fits1[[i]]
		plot(c(1, 5), c(min(mydat$CILower95), max(mydat$CIUpper95)), xlab = '', ylab = '', axes = FALSE, 
			type = 'n')
		if(i == 1){
			axis(1, at = 1:5, padj = 1, mgp = c(3, 0.2, 0), labels = FALSE)
			mtext(side = 2, bquote('Invertebrate density (#*'*m^-2*')'), line = 3.2)
			legend('topright', legend = c(2011, 2016), pt.cex = 1.5, pch = c(17, 16), col = c(2, 4), 
				bty = 'n')
		}
		if(i == 2){
			axis(1, at = 1:5, padj = 1, mgp = c(3, 0.2, 0), labels = FALSE)
			mtext(side = 2, 'Mean invertebrate length (mm)', line = 3.2)
		}
		if (i == 3){
			axis(1, at = 1:5, padj = 1, mgp = c(3, 0.2, 0), labels = c('Scraper/\nGrazer', 
				'Collector\nFilterer', 'Collector\nGatherer', 'Generalist', 'Predator'))
			mtext(side = 2, bquote('Invertebrate biomass (mg*'*m^-2*')'), line = 3.2)
		}
		axis(2, las = 2)
		box(bty = 'l')
		with(mydat[mydat$Study == 'Pre',], arrows(x0 = (1:5) - d1, y0 = CILower95, y1 = CIUpper95, 
			code = 3, angle = 90, length = 0.05))
		with(mydat[mydat$Study == 'Post',], arrows(x0 = (1:5) + d1, y0 = CILower95, y1 = CIUpper95, 
			code = 3, angle = 90, length = 0.05))
		with(mydat[mydat$Study == 'Pre',], arrows(x0 = (1:5) - d1, y0 = CILower50, y1 = CIUpper50, 
			length = 0, col = 'grey50', lwd = 4))
		with(mydat[mydat$Study == 'Post',], arrows(x0 = (1:5) + d1, y0 = CILower50, y1 = CIUpper50, 
			length = 0, col = 'grey50', lwd = 4))
		points((1:5) - d1, mydat[mydat$Study == 'Pre', 'Fit'], pch = 17, col = 2, cex = 1.5)
		points((1:5) + d1, mydat[mydat$Study == 'Post', 'Fit'], pch = 16, col = 4, cex = 1.5)
		legend('topleft', legend = LETTERS[i], bty = 'n')
	}
}
plotTypes(plotmod, 'ModelDensitySizeBiomassFFG', 'Figures')


##### Combine overall and FFG plots into a single panel graph #####

## Use compilation of plotpred and plotmod
plotall <- function(){
	layout(matrix(1:6, 3, 2), widths = c(1, 3))
	par(mar = c(2, 2.8, 0.2, 1), oma = c(1, 2, 0, 0), cex = 1)
	for(i in 1:length(fits0)){
		mydat <- fits0[[i]]
		plot(c(0.5, 2.5), c(min(mydat$CILower95), max(mydat$CIUpper95)), xlab = '', ylab = '', 
			axes = FALSE, type = 'n')
		if(i != 3){axis(1, at = 1:2, labels = FALSE)
		} else{axis(1, at = 1:2, labels = c('2011', '2016'))
		}
		axis(2, las = 2)
		if(i == 1){mtext(side = 2, bquote('Invertebrate density (#*'*m^-2*')'), line = 3.5)}
		if(i == 2){mtext(side = 2, 'Mean invertebrate size (mm)', line = 3.5)}
		if(i == 3){mtext(side = 2, bquote('Invertebrate biomass (mg*'*m^-2*')'), line = 3.5)}
		box(bty = 'l')
		with(mydat, arrows(x0 = 1:2, y0 = CILower95, y1 = CIUpper95, code = 3, angle = 90, length = 0.05))
		with(mydat, arrows(x0 = 1:2, y0 = CILower50, y1 = CIUpper50, length = 0, col = 'grey50', lwd = 4))
		points(mydat$Fit, pch = c(17, 16), col = c(2, 4), cex = 1.5)
		legend('top', LETTERS[i], bty = 'n')
	}
	d1 <- 0.09
	for(i in 1:length(fits1)){
		mydat <- fits1[[i]]
		plot(c(1, 5), c(min(mydat$CILower95), max(mydat$CIUpper95)), xlab = '', ylab = '', axes = FALSE, 
			type = 'n')
		if(i == 1){
			axis(1, at = 1:5, padj = 1, mgp = c(3, 0.2, 0), labels = FALSE)
			legend('topright', legend = c(2011, 2016), pt.cex = 1.5, pch = c(17, 16), col = c(2, 4), 
				bty = 'n')
		}
		if(i == 2){
			axis(1, at = 1:5, padj = 1, mgp = c(3, 0.2, 0), labels = FALSE)
		}
		if (i == 3){
			axis(1, at = 1:5, padj = 1, mgp = c(3, 0.2, 0), labels = c('Scraper/\nGrazer', 
				'Collector\nFilterer', 'Collector\nGatherer', 'Generalist', 'Predator'))
		}
		axis(2, las = 2)
		box(bty = 'l')
		with(mydat[mydat$Study == 'Pre',], arrows(x0 = (1:5) - d1, y0 = CILower95, y1 = CIUpper95, 
			code = 3, angle = 90, length = 0.05))
		with(mydat[mydat$Study == 'Post',], arrows(x0 = (1:5) + d1, y0 = CILower95, y1 = CIUpper95, 
			code = 3, angle = 90, length = 0.05))
		with(mydat[mydat$Study == 'Pre',], arrows(x0 = (1:5) - d1, y0 = CILower50, y1 = CIUpper50, 
			length = 0, col = 'grey50', lwd = 4))
		with(mydat[mydat$Study == 'Post',], arrows(x0 = (1:5) + d1, y0 = CILower50, y1 = CIUpper50, 
			length = 0, col = 'grey50', lwd = 4))
		points((1:5) - d1, mydat[mydat$Study == 'Pre', 'Fit'], pch = 17, col = 2, cex = 1.5)
		points((1:5) + d1, mydat[mydat$Study == 'Post', 'Fit'], pch = 16, col = 4, cex = 1.5)
		if(par()$mfrow[2] == 2){legend('topleft', legend = LETTERS[i + 3], bty = 'n')
		} else{legend('topleft', legend = LETTERS[i + 3], bty = 'n')}
	}
}
plotTypes(plotall, 'ModelDensitySizeBiomassAll', 'Figures', width = 8, height = 8)


##### Ordination analysis #####

## Put Whiting and our benthic data in a single dataframe
ord1 <- rbind(dat5$Whiting, dat5$Benthic)
	ord1$Density <- round(ord1$Count / ord1$Unit)
	ord1$Study <- factor(ifelse(year(ord1$Date) < 2015, 'Pre', 'Post'), levels = c('Pre', 'Post'))

## Get data in an ordination-friendly matrix
ord2 <- matrix(nrow = length(unique(ord1$BarcodeID)), ncol = length(unique(ord1$SpeciesID)))
	rows <- rownames(ord2) <- sort(unique(ord1$BarcodeID))
	cols <- colnames(ord2) <- sort(unique(ord1$SpeciesID))
	ord2[is.na(ord2)] <- 0
rownum <-match(ord1$BarcodeID, rows)
colnum <- match(ord1$SpeciesID, cols)
for(i in 1 : dim(ord1)[1]){
	ord2[rownum[i], colnum[i]] <- ord1$Density[i]
}

## Check cvs
cv2 <- cv(ord2)
	## 86% by row and 162% by columns. Not horrendous.

## Remove rare species
ord3 <- delRare(ord2)
	delspp <- colnames(ord2)[colnames(ord2) %in% colnames(ord3) == FALSE]
	## Removed 1 taxon (LIBE).
cv3 <- cv(ord3)
	## Dropped cv on columns to 158%.

## Relativize by species max
ord4 <- rel(ord3)
cv4 <- cv(ord4)
	## Dropped cv on columns to 50%, and on rows to 79%.
	## Probably OK for analysis from this point.

## Run an NMS with stepdown
#NMS(ord4, maxruns = 10000)
	## Scree plot is typical; suggests 2-3 dimensions would be best.
	## Going to work with 2D for ease.
	## Note that this only needs to be run once, hence why it's greyed out now.
pts2D <- read.csv('NMS Output/NMSPoints2D.csv', row.names = 1)
spp2D <- read.csv('NMS Output/NMSSpecies2D.csv', row.names = 1)

## Check axes R2 of 2D solution
r2 <- axisR2(ord4, pts2D)
	## 86% of total variation is explained by the 2 axis (73% on Axis 1). Pretty good.

## Build a matrix of grouping variables
unqrow <- match(unique(rownum), rownum)
env1 <- ord1[unqrow, c('Study', 'Season')]
	rownames(env1) <- ord1[unqrow, 'BarcodeID']

## See how groupings load on the ordination
envpts <- envfit(pts2D, env1)
	## 30% for pre-post and 12% for season. Neither is great.

## Plot groupings on ordination
plotord <- function(){
	pchs <- ifelse(as.numeric(env1$Season) == 1, 9,
		ifelse(as.numeric(env1$Season) == 2, 10,
		ifelse(as.numeric(env1$Season) == 3, 12, 13)))
	par(mfrow = c(1, 1), mar = c(4, 4, 0.1, 0.1), cex = 1)
	plot(c(-2, 2), c(-2, 2), type = 'n', xlab = 'Axis 1', ylab = 'Axis 2', axes = FALSE)
	axis(1)
	axis(2, las = 2)
	box(bty = 'l')
	points(pts2D, col = ifelse(env1$Study == 0, 2, 4), pch = pchs)
	ordiellipse(pts2D, env1$Study, kind = 'sd', conf = 0.95, col = c(2, 4))
	text(spp2D, rownames(spp2D), col = 1, cex = 0.6)
	legend('topright', legend = c('2011', '2016', levels(env1$Season), '95% CI'), 
		col = c(2, 4, rep(1, 4), 8), pch = c(rep(15, 2), sort(unique(pchs)), 1), 
		pt.cex = c(1.5, 1.5, rep(1, 4), 1.5), bty = 'n')
}
plotTypes(plotord, 'Ordination', 'Figures', height = 6.5)
	## See the loss of predators
	### STOPPED HERE.
	## NOTE: Consider groups spp2D by FFG instead and showing that instead of sample points 
		## (bigger deal, and less cluttered).

## See how taxa load on the ordination
envspp <- envfit(pts2D, ord4)
envspp1 <- data.frame(round(cbind(envspp$vectors[[1]], envspp$vectors[[2]], envspp$vectors[[4]]), 4))
	colnames(envspp1) <- c('Axis1', 'Axis2', 'R2', 'p')
	envspp1 <- envspp1[order(-envspp1$R2),]

## Limit to only taxa with R2 >= 25%
envspp2 <- envspp1[envspp1$R2 >= 0.25,]
	envspp2$FFG <- spp[match(rownames(envspp2), spp$SpeciesID), 'FFG']
	## TABL, PETL, PLAN, MEGL, OCHR, HELI are most strongly loading on the ordination.
	## No real theme to these in terms of FFG, although most are loading in the 
		## direction of Whiting's samples (MEGL the opposite).


##### Look at differences in sample richness and diversity #####

## Get data matrices for Whiting and our Benthic data
ordlist <- list(Whiting = ord2[substr(rownames(ord2), 1, 1) != 'B', ], 
	Benthic = ord2[substr(rownames(ord2), 1, 1) == 'B', ])
	ordlist <- lapply(ordlist, function(x) x[, which(colSums(x) > 0)])

## Get overall richness and diversity
ordlist0 <- lapply(ordlist, function(x){
	rownames(x) <- rep(1, dim(x)[1])
	return(x)
	})
rich0 <- data.frame(t(sapply(ordlist0, function(x) data.frame(
	AllRich = dim(x)[2], AllShann = round(diversity(colSums(x)), 2), 
	AllSimp = round(diversity(colSums(x), 'simpson'), 2)))))
	
## Get sample-level richness and diversity
rich1 <- lapply(ordlist, function(x) data.frame(
	Richness = rowSums(x != 0), Shannon = round(diversity(x), 4), 
	Simpson = round(diversity(x, 'simpson'), 4)))
rich0$MeanRich <- sapply(rich1, function(x) round(mean(x$Richness), 2))
rich0$MeanShann <- sapply(rich1, function(x) round(mean(x$Shannon), 2))
rich0$MeanSimp <- sapply(rich1, function(x) round(mean(x$Simpson), 2))

## Model richness and diversity, using normal distribution
rich2 <- rbind(rich1$Whiting, rich1$Benthic)
rich2$Season <- samp1[match(rownames(rich2), samp1$BarcodeID), 'Season']
rich2$Study <- samp1[match(rownames(rich2), samp1$BarcodeID), 'Study']
richmod <- list()
	richmod$Richness$r1 <- glmmTMB(Richness ~ 1, data = rich2)
	richmod$Richness$r2 <- glmmTMB(Richness ~ 1 + (1 | Season), data = rich2)
	richmod$Richness$r3 <- glmmTMB(Richness ~ 0 + Study + (1 | Season), data = rich2)
	richmod$Shannon$r1 <- glmmTMB(Shannon ~ 1, data = rich2)
	richmod$Shannon$r2 <- glmmTMB(Shannon ~ 1 + (1 | Season), data = rich2)
	richmod$Shannon$r3 <- glmmTMB(Shannon ~ 0 + Study + (1 | Season), data = rich2)
	richmod$Simpson$r1 <- glmmTMB(Simpson ~ 1, data = rich2)
	richmod$Simpson$r2 <- glmmTMB(Simpson ~ 1 + (1 | Season), data = rich2)
	richmod$Simpson$r3 <- glmmTMB(Simpson ~ 0 + Study + (1 | Season), data = rich2)	
	names(richmod) <- colnames(rich2)[1 : 3]
richmodAIC <- AICtable(richmod)
	## Adding a season random effect doesn't really affect models; adding pre-post really helps.

## Get best model for richness and diversity
bestmod2 <- bestmod(richmod)

## Get parameters to predict
predparm2 <- data.frame(Study = levels(rich2$Study), Season = 'Generic')
	
## Get model-predicted values and confidence intervals
preds2 <- lapply(bestmod2, predict, newdata = predparm2, type = 'response', se.fit = TRUE, 
	allow.new.levels = TRUE)
cis2v95 <- lapply(bestmod2, confint)
cis2v50 <- lapply(bestmod2, confint, level = 0.5)
fits2 <- list()
for(i in 1 : length(bestmod2)){
	fits2[[i]] <- data.frame(Study = levels(samp1$Study))
	fits2[[i]]$Fit <- preds2[[i]]$fit
	fits2[[i]]$SELower <- preds2[[i]]$fit - preds2[[i]]$se.fit
	fits2[[i]]$SEUpper <- preds2[[i]]$fit + preds2[[i]]$se.fit
	fits2[[i]]$CILower95 <- cis2v95[[i]][1:2, 1]
	fits2[[i]]$CIUpper95 <- cis2v95[[i]][1:2, 2]
	fits2[[i]]$CILower50 <- cis2v50[[i]][1:2, 1]
	fits2[[i]]$CIUpper50 <- cis2v50[[i]][1:2, 2]
}
	fits2 <- lapply(fits2, function(x){cbind(Study = x[, 1], round(x[, -1], 2))})
	names(fits2) <- names(bestmod2)
	## Decline in all metrics from pre to post.

## Plot panel graph of overall change in Richness and Shannon's and Simpson's diversity
plotpred2 <- function(){
	ylabs <- c('Richness (taxa per sample)', "Shannon's diversity (H)", "Simpson's diversity (D)")
	par(mfrow = c(3, 1), mar = c(2, 5, 0.1, 0.1), cex = 1)	
	for(i in 1:3){
		mydat <- fits2[[i]]
		plot(c(0.5, 2.5), c(min(mydat$CILower95), max(mydat$CIUpper95)), xlab = '', ylab = '', axes = FALSE, 
			type = 'n')
		if(i != 3){axis(1, at = 1:2, labels = FALSE)
		} else {axis(1, at = 1:2, labels = c('2011', '2016'))}
		axis(2, las = 2)
		mtext(side = 2, ylabs[i], line = 3.2)
		box(bty = 'l')
		with(mydat, arrows(x0 = 1:2, y0 = CILower95, y1 = CIUpper95, code = 3, angle = 90, length = 0.05))
		with(mydat, arrows(x0 = 1:2, y0 = CILower50, y1 = CIUpper50, length = 0, col = 'grey50', lwd = 4))
		points(mydat$Fit, pch = c(17, 16), col = c(2, 4), cex = 1.5)
		legend('topright', legend = LETTERS[i], bty= 'n')
	}
}
plotTypes(plotpred2, 'ModelRichnessDiversity', 'Figures', width = 4)
	
##### Plot species accumulation curves #####

## Get species accumulation for Whiting and benthic data (1000 m only)
spac0 <- suppressWarnings(lapply(ordlist, specaccum))

## Compare species accumulation from all benthic data (not just 1000 m)
	## As a synoptic comparison
ord1all <- rbind(dat4$Whiting, dat4$Benthic)
ord2all <- matrix(nrow = length(unique(ord1all$BarcodeID)), ncol = length(unique(ord1all$SpeciesID)))
	rows <- rownames(ord2all) <- sort(unique(ord1all$BarcodeID))
	cols <- colnames(ord2all) <- sort(unique(ord1all$SpeciesID))
	ord2all[is.na(ord2all)] <- 0
rownumall <-match(ord1all$BarcodeID, rows)
colnumall <- match(ord1all$SpeciesID, cols)
for(i in 1 : dim(ord1all)[1]){
	ord2all[rownumall[i], colnumall[i]] <- ord1all$Count[i]
}
ordlist1 <- list(Whiting = ord2all[substr(rownames(ord2all), 1, 1) != 'B', ], 
	Benthic = ord2all[substr(rownames(ord2all), 1, 1) == 'B', ])
	ordlist1 <- lapply(ordlist1, function(x) x[, which(colSums(x) > 0)])
spac1 <- suppressWarnings(lapply(ordlist1, specaccum))

plotaccum <- function(){
	par(mar = c(4, 4, 0.4, 0.6), cex = 1, fig = c(0, 1, 0, 1), xaxs = 'i', yaxs = 'i')
	plot(spac0[[1]], col = 2, lwd = 1, ci.type = 'poly', ci.col = 'lightpink', ci.lty = 0, 
		xlab = 'Number of samples', ylab = 'Cumulative taxa', bty = 'n', axes = FALSE, xlim = c(0, 25),
		ylim = c(0, 30))
	plot(spac0[[2]], add = TRUE, col = 4, lwd = 1, ci.type = 'poly', ci.col = 'lightblue', ci.lty = 0)
	axis(1, at = seq(0, 25, 5))
	axis(2, las = 2)
	box(bty = 'l')
	legend(17, 21, legend = c('2011', '2016', '', '2016 synoptic', '(Creek-wide)'), 
		col = c(2, 4, NA, 4, NA), lwd = c(1, 1, 0, 1, 0), lty = c(1, 1, 0, 2, 0), bty = 'n')
	par(fig = c(0.57, 0.97, 0.1, 0.5), xpd = TRUE, new = TRUE)
	plot(spac0[[1]], col = 2, lwd = 1, ci.type = 'poly', ci.col = 'lightpink', ci.lty = 0, 
		xlab = '', ylab = '', bty = 'n', axes = FALSE, xlim = c(0, 80), ylim = c(0, 30))
	plot(spac1[[2]], add = TRUE, col = 4, lwd = 1, lty = 2, ci.type = 'poly', ci.col = 'lightblue', 
		ci.lty = 0)
	axis(1, padj = -1, at = seq(0, 80, 20), tcl = -0.35)
	axis(1, padj = -1, at = seq(10, 70, 20), labels = FALSE, tick = TRUE, tcl = -0.2)
	axis(2, las = 2, hadj = 0.5, at = seq(0, 30, 10), tcl = -0.35)
	axis(2, at = seq(5, 25, 10), labels = FALSE, tick = TRUE, tcl = -0.2)	
	mtext(side = 1, '# samples', line = 1.4)
	mtext(side = 2, '# taxa', line = 1.7)
	box(bty = 'l')
}
plotTypes(plotaccum, 'SpeciesAccumulation', 'Figures', height = 6.5)
	## Note: Warnings here have to do with plotting and are overridden. Can be ignored.

### NEXT STEPS: 
	## Could be interesting to see if predators have gotten bigger, for instance 
		## (are large Corydalus and Odonates invulnerable to dace predation?)