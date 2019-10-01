##### Bright Angel trophic cascade analysis #####
## Last updated 1 October 2019 by J.D. Muehlbauer


##### Set up workspace ##### 

## Load/install requisite packages
source('https://github.com/jmuehlbauer-usgs/R-packages/blob/master/packload.r?raw=TRUE')
packload(c('devtools', 'lubridate', 'plots', 'bugR', 'MASS', 'glmmTMB'))

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
	mycols <- c('BarcodeID', 'Date', 'Season', ifelse('Density' %in% colnames(x), 'Area', 'Volume'),
		'SpeciesID', 'FFG', 'CountTotal', 'Size', 'Biomass')
	x[, mycols]
	})
dat0 <- lapply(dat0, setNames, nm = c('BarcodeID', 'Date', 'Season', 'Unit', 'SpeciesID', 'FFG', 
	'Count', 'Size', 'Biomass'))

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
		tCnt <- aggregate(t1$Count, by = list(t1$BarGrp), sum)
		tSize <- aggregate(t1$Count * t1$Size, by = list(t1$BarGrp), sum)
			tSize[, 2] <- round(tSize[, 2] / tCnt[, 2], 2)
		tBio <- aggregate(t1$Biomass, by = list(t1$BarGrp), sum)
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
dat4 <- lapply(dat4, function(x){
	x[,'FFG'] <- spp[match(x[, 'SpeciesID'], spp$SpeciesID), 'FFG']
	return(x)
	})

# Convert SpeciesIDs back to factors
dat4 <- lapply(dat4, transform, SpeciesID = as.factor(SpeciesID))

## Sort by Date, BarcodeID, SpeciesID, drop old factor levels
dat4 <- lapply(dat4, function(x) x[order(x[, 'Date'], x[, 'BarcodeID'], x[, 'SpeciesID']),])
	dat4 <- lapply(dat4, droplevels)

## Add consistent factor levels
levspp <- sort(unique(unlist(lapply(dat4, function(x) levels(x[, 'SpeciesID'])))))
levFFG <- unique(unlist(lapply(dat4, function(x) levels(x[, 'FFG']))))
	levFFG <- ifelse(length(levFFG) == 6, c('Shredder', 'ScraperGrazer', 'CollectorFilterer', 
		'CollectorGatherer', 'Generalist', 'Predator'), levFFG)
dat4 <- lapply(dat4, function(x){
	levspp1 <- levspp[!(levspp %in% levels(x[, 'SpeciesID']))]
	levFFG1 <- levFFG[!(levFFG %in% levels(x[, 'FFG']))]
	levels(x[, 'SpeciesID']) <- c(levels(x[, 'SpeciesID']), levspp[!(levspp %in% levels(x[, 'SpeciesID']))])
	levels(x[, 'FFG']) <- c(levels(x[, 'FFG']), levFFG[!(levFFG %in% levels(x[, 'FFG']))])
	x[, 'SpeciesID'] <- factor(x[, 'SpeciesID'], levels = levspp)
	x[, 'FFG'] <- factor(x[, 'FFG'], levels = c('Shredder', 'ScraperGrazer', 'CollectorFilterer',
		'CollectorGatherer', 'Generalist', 'Predator'))
	return(x)
	})


##### Group specimens by Sample #####

## Create new list, with samples summed.
samp <- combinefx(dat4, 'BarcodeID')
	samp <- lapply(samp, subset, select = -SpeciesID)

## Combine benthic and Whiting into single dataframe for analysis
samp1 <- rbind(samp$Whiting, samp$Benthic)
	samp1$Study <- factor(ifelse(year(samp1$Date) < 2015, 'Pre', 'Post'), levels = c('Pre', 'Post'))

## Build some models for count data
sampmod <- list()
	sampmod[[1]] <- sampmod[[2]] <- sampmod[[3]] <- list()
	names(sampmod) <- c('Count', 'Size', 'Biomass')
	
## Find best data distribution (family) to use
sampmod[[1]]$lm1 <- glm(Count ~ 1 + offset(Unit), data = samp1)
sampmod[[1]]$po1 <- glm(Count ~ 1 + offset(log(Unit)), family = 'poisson', data = samp1)
sampmod[[1]]$nb1 <- glm.nb(Count ~ 1 + offset(log(Unit)), data = samp1)

## Add a random effect for sample or season
sampmod[[1]]$nb2 <- suppressWarnings(glmmTMB(Count ~ 1 + (1 | BarcodeID) + offset(log(Unit)), 
	family = 'nbinom2', data = samp1))
sampmod[[1]]$nb3 <- suppressWarnings(glmmTMB(Count ~ 1 + (1 | Season) + offset(log(Unit)), 
	family = 'nbinom2', data = samp1))
sampmod[[1]]$nb4 <- suppressWarnings(glmmTMB(Count ~ 1 + (1 | BarcodeID) + (1 | Season) + 
	offset(log(Unit)), family = 'nbinom2', data = samp1))
	
## Add pre-post fixed effect
sampmod[[1]]$nb5 <- suppressWarnings(glmmTMB(Count ~ Study + (1 | BarcodeID) + offset(log(Unit)), 
	family = 'nbinom2', data = samp1))
sampmod[[1]]$nb6 <- suppressWarnings(glmmTMB(Count ~ Study + (1 | Season) + offset(log(Unit)), 
	family = 'nbinom2', data = samp1))
sampmod[[1]]$nb7 <- suppressWarnings(glmmTMB(Count ~ Study + (1 | BarcodeID) + (1 | Season) + 
	offset(log(Unit)), family = 'nbinom2', data = samp1))	

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
	## A Season random effect fits best and makes logical sense.
	## Including a pre-post term improves the Count (density) model.

## Build models for size and biomass data
	## Using normal distribution for sizes, based on histogram.
	## Using Gamma distribution for biomass (continuous, right-tailed data).
	## Keeping random effect as-is from count model.
sampmod[[2]]$lm3 <- suppressWarnings(glmmTMB(Size ~ 1 + (1 | Season) + offset(log(Unit)), 
	family = 'gaussian', data = samp1))
sampmod[[2]]$lm6 <- suppressWarnings(glmmTMB(Size ~ Study + (1 | Season) + offset(log(Unit)), 
	family = 'gaussian', data = samp1))	
sampmod[[3]]$ga3 <- suppressWarnings(glmmTMB(Biomass ~ 1 + (1 | Season) + offset(log(Unit)), 
	family = Gamma(link = 'log'), data = samp1))
sampmod[[3]]$ga6 <- suppressWarnings(glmmTMB(Biomass ~ Study + (1 | Season) + offset(log(Unit)), 
	family = Gamma(link = 'log'), data = samp1))	
sampmodAIC <- AICtable(sampmod)
	## Adding a pre-post fixed effect improves Biomass and Size models as well.

## Get best model for Count, Size, and Biomass
bestmod <- function(mods){
	lapply(mods, function(x) {
	x[names(which(sapply(x, AIC) == min(sapply(x, AIC), na.rm = TRUE)))][[1]]
	})
}
bestmod0 <- bestmod(sampmod)

## Get parameters to predict
predparm0 <- data.frame(Study = levels(samp1$Study), Season = 'Generic', Unit = 1)
	
## Get model-predicted values and confidence intervals
preds0l <- lapply(bestmod0, predict, newdata = predparm0, type = 'link', se.fit = TRUE, 
	allow.new.levels = TRUE)
preds0r <- lapply(bestmod0, predict, newdata = predparm0, type = 'response', se.fit = TRUE, 
	allow.new.levels = TRUE)
fits0 <- list()
for(i in 1 : length(bestmod0)){
	fits0[[i]] <- data.frame(Study = levels(samp1$Study), Unit = 1)
	fits0[[i]]$LinkFit <- preds0l[[i]]$fit
	fits0[[i]]$LinkSE <- preds0l[[i]]$se.fit
	fits0[[i]]$Fit <- round(preds0r[[i]]$fit, 2)
	fits0[[i]]$SELower <- round(preds0r[[i]]$fit - preds0r[[i]]$se.fit, 2)
	fits0[[i]]$SEUpper <- round(preds0r[[i]]$fit + preds0r[[i]]$se.fit, 2)
	if(bestmod0[[i]][6]$modelInfo$family$link == 'log'){
		fits0[[i]]$CILower <- round(exp(fits0[[i]]$LinkFit - (qnorm(0.975) * fits0[[i]]$LinkSE)), 2)
		fits0[[i]]$CIUpper <- round(exp(fits0[[i]]$LinkFit + (qnorm(0.975) * fits0[[i]]$LinkSE)), 2)	
	} else{
		fits0[[i]]$CILower <- round(fits0[[i]]$LinkFit - (qnorm(0.975) * fits0[[i]]$LinkSE), 2)
		fits0[[i]]$CIUpper <- round(fits0[[i]]$LinkFit + (qnorm(0.975) * fits0[[i]]$LinkSE), 2)	
	}
}
	names(fits0) <- names(bestmod0)

## Plot panel graph of overall change in Density, Size, and Biomass
plotpred <- function(){
	par(mfrow = c(3, 1), mar = c(3, 5, 0.2, 1), cex = 1)
	for(i in 1:3){
		mydat <- fits0[[i]]
		plot(c(0.5, 2.5), c(min(mydat$SELower), max(mydat$SEUpper)), xlab = '', ylab = '', axes = FALSE, 
			type = 'n')
		if(i != 3){axis(1, at = 1:2, labels = FALSE)
		} else{axis(1, at = 1:2, labels = c('2011', '2016'))
		}
		axis(2, las = 2)
		if(i == 1){mtext(side = 2, bquote('Invertebrate density ('*m^-2*')'), line = 3.5)}
		if(i == 2){mtext(side = 2, 'Mean invertebrate size (mm)', line = 3.5)}
		if(i == 3){mtext(side = 2, bquote('Invertebrate biomass (mg*'*m^-2*')'), line = 3.5)}
		box(bty = 'l')
		points(mydat$Fit, pch = 16, cex = 1.5)
		with(mydat, arrows(x0 = 1:2, y0 = SELower, y1 = SEUpper, code = 3, angle = 90, length = 0.05))
	}
}
plotTypes(plotpred, 'ModelOverall', 'Figures', width = 4)
	

##### Group specimens by FFG #####

## Create new list, with samples summed by FFG
ffg <- combinefx(dat4, 'FFG')
	ffg <- lapply(ffg, subset, select = -SpeciesID)

## Get relative counts and biomasses
ffg <- lapply(ffg, function(x){
	cnttot <- tapply(x[, 'Count'], x[, 'BarcodeID'], sum)
	biotot <- tapply(x[, 'Biomass'], x[, 'BarcodeID'], sum)
	x[, 'RelCount'] <- round(x[, 'Count'] / cnttot[match(x[, 'BarcodeID'], names(cnttot))], 4)
	x[, 'RelBiomass'] <- round(x[, 'Biomass'] / biotot[match(x[, 'BarcodeID'], names(biotot))], 4)
	return(x)
	})

## Get mean, standard error, and relative mean by FFG, by season
ffgcnt <- lapply(ffg, function(x) round(tapply(x[, 'Count'], list(x[, 'FFG'], x[, 'Season']), mean)))
ffgbio <- lapply(ffg, function(x) round(tapply(x[, 'Biomass'], list(x[, 'FFG'], x[, 'Season']), mean), 2))
ffgsize <- lapply(ffg, function(x) round(tapply(x[, 'Size'], list(x[, 'FFG'], x[, 'Season']), mean), 2))
#ffgsem <- lapply(ffg, function(x) round(tapply(x[, 'Count'], list(x[, 'FFG'], x[, 'Season']), 
	#function(x) sd(x) / sqrt(length(x))), 4))
ffgunitcnt <- lapply(ffg, function(x) 
	round(tapply(x[, 'Count'] / x[, 'Unit'], list(x[, 'FFG'], x[, 'Season']), mean), 4))
ffgunitbio <- lapply(ffg, function(x) 
	round(tapply(x[, 'Biomass'] / x[, 'Unit'], list(x[, 'FFG'], x[, 'Season']), mean), 4))
ffgrelcnt <- lapply(ffg, function(x) round(tapply(x[, 'RelCount'], list(x[, 'FFG'], x[, 'Season']), 
	mean), 4))
ffgrelbio <- lapply(ffg, function(x) round(tapply(x[, 'RelBiomass'], list(x[, 'FFG'], x[, 'Season']), 
	mean), 4))
ffgstat <- list(ffgcnt, ffgbio, ffgsize, ffgunitcnt, ffgunitbio, ffgrelcnt, ffgrelbio)
	names(ffgstat) <- c('MeanCount', 'MeanBiomass', 'MeanSize', 'MeanUnitCount', 'MeanUnitBiomass', 'MeanRelCount', 'MeanRelBiomass')
	ffgstat <- lapply(ffgstat, function(x){
		lapply(x, function(y){
			y[is.na(y)] <- 0
			return(y)
		})
	})	


##### Look at differences by FFG over time #####

## Get difference in absolute and relative density/abundance, biomass, and size from Whiting to our study
diff1 <- lapply(ffgstat[c('MeanUnitCount', 'MeanRelCount', 'MeanUnitBiomass', 'MeanRelBiomass', 
	'MeanSize')], 
	function(x) x['Benthic'][[1]] - x['Whiting'][[1]])
diff2 <- lapply(diff1, function(x){
	data.frame(Mean = round(rowMeans(x), 4), 
		SEM = round(apply(x, 1, function(y) sd(y) / sqrt(length(y))), 4))
	})

## Plot differences
plotdiff <- function(){
	par(mfcol = c(2, 3), mar = c(2, 5, 0.2, 1.2), cex = 0.8, oma = c(3, 0, 0.2, 0))
	for(i in 1:5){
		plot(c(1, 6), c(min(diff2[[i]]$Mean - diff2[[i]]$SEM), max(diff2[[i]]$Mean + diff2[[i]]$SEM)), 
		xlab = '', ylab = '', axes = FALSE, type = 'n')
		if(i %in% c(1, 3)){
			axis(1, at = 1:6, padj = 1, mgp = c(3, 0.2, 0), labels = FALSE)
		} else{
			axis(1, las = 2, at = 1:6, labels = c('Shredder', 
				'Scraper/\nGrazer', 'Collector\nFilterer', 'Collector\nGatherer', 'Generalist', 'Predator'))
		}
		axis(2, las = 2)
		if(i == 1){mtext(side = 2, bquote('Difference in density ('*m^-2*')'), line = 3.4, cex = 0.8)}
		if(i == 2){mtext(side = 2, 'Difference in relative abundance', line = 3.6, cex = 0.8)}
		if(i == 3){mtext(side = 2, bquote('Difference in biomass ('*mg*'*'*m^-2*')'), line = 3.4, 
			cex = 0.8)}
		if(i == 4){mtext(side = 2, 'Difference in relative biomass', line = 3.6, cex = 0.8)}
		if(i == 5){mtext(side = 2, 'Difference in lengths (mm)', line = 3.4, cex = 0.8)}
	abline(h = 0, lty = 2)
	box(bty = 'l')
	points(diff2[[i]]$Mean, pch = 16)
	with(diff2[[i]], arrows(x0 = 1:6, y0 = Mean - SEM, y1 = Mean + SEM, code = 3, angle = 90, 
		length = 0.05))
	}
}
plotTypes(plotdiff, 'ObservedDifferences', 'Figures', height = 5.3, width = 9)
	## Panels are pretty similar, which is good.
	
	
##### Build model for FFG change #####

## Combine benthic and Whiting into single dataframe for analysis
ffg1 <- rbind(ffg$Whiting, ffg$Benthic)
	ffg1$Study <- factor(ifelse(year(ffg1$Date) < 2015, 'Pre', 'Post'), levels = c('Pre', 'Post'))

## Build some models
ffgmod <- list()
for(i in 1 : length(levels(ffg1$FFG))){
	mydat <- ffg1[ffg1$FFG == levels(ffg1$FFG)[i],]
	ffgmod[[i]] <- list()
	## Find best data distribution (family) to use
	ffgmod[[i]]$lm1 <- glm(Count ~ 1 + offset(Unit), data = mydat)
	ffgmod[[i]]$po1 <- glm(Count ~ 1 + offset(log(Unit)), family = 'poisson', data = mydat)
	ffgmod[[i]]$nb1 <- glm.nb(Count ~ 1 + offset(log(Unit)), data = mydat)
	## Add a random effect for sample or season
	ffgmod[[i]]$nb2 <- suppressWarnings(glmmTMB(Count ~ 1 + (1 | BarcodeID) + offset(log(Unit)), 
		family = 'nbinom2', data = mydat))
	ffgmod[[i]]$nb3 <- suppressWarnings(glmmTMB(Count ~ 1 + (1 | Season) + offset(log(Unit)), 
		family = 'nbinom2', data = mydat))
	ffgmod[[i]]$nb4 <- suppressWarnings(glmmTMB(Count ~ 1 + (1 | BarcodeID) + (1 | Season) + 
		offset(log(Unit)), family = 'nbinom2', data = mydat))
	## Add pre-post and FFG fixed effects
	ffgmod[[i]]$nb5 <- suppressWarnings(glmmTMB(Count ~ Study + (1 | BarcodeID) + offset(log(Unit)), 
		family = 'nbinom2', data = mydat))
	ffgmod[[i]]$nb6 <- suppressWarnings(glmmTMB(Count ~ Study + (1 | Season) + offset(log(Unit)), 
		family = 'nbinom2', data = mydat))
	ffgmod[[i]]$nb7 <- suppressWarnings(glmmTMB(Count ~ Study + (1 | BarcodeID) + (1 | Season) + 
		offset(log(Unit)), family = 'nbinom2', data = mydat))	
}
	## Note: warnings here are just due to FFGs with NAs that therefore don't fit. Can be ignored.
	names(ffgmod) <- levels(ffg1$FFG)
	
## Compare models using AIC
ffgmodAIC <- AICtable(ffgmod)
	## Negative binomial distribution always the best bet (over linear/Gaussian and Poisson).
	## The choice of random effect or whether to include one is less clear.
		## However, BarcodeID doesn't always fit, and a random effect for season is sensible. Use that.
	## Including a pre-post term improves all models except for Collector-Filterers and -Gatherers.
		## Critical to the analysis though, so use it anyway.

## Get parameters to predict
tlen <- length(levels(ffg1$Study))
flen <- length(levels(ffg1$FFG))
predparm <- data.frame(Study = levels(ffg1$Study), Season = 'Generic', Unit = 1)
	
## Get model-predicted values and confidence intervals
predsl <- lapply(ffgmod, function(x){
	predict(x[['nb6']], predparm, type = 'link', se.fit = TRUE, allow.new.levels = TRUE)
	})
predsr <- lapply(ffgmod, function(x){
	predict(x[['nb6']], predparm, type = 'response', se.fit = TRUE, allow.new.levels = TRUE)
	})
fits1 <- data.frame(FFG = rep(levels(ffg1$FFG), tlen), 
	Study = rep(levels(ffg1$Study), rep(flen, tlen)), Unit = 1)
	fits1$LinkFit <- c(t(sapply(predsl, function(x) x[['fit']])))
	fits1$LinkSE <- c(t(sapply(predsl, function(x) x[['se.fit']])))
	fits1$Density <- c(t(sapply(predsr, function(x) x[['fit']])))
	fits1$SELower <- c(t(sapply(predsr, function(x) x[['fit']]))) - 
		c(t(sapply(predsr, function(x) x[['se.fit']])))
	fits1$SEUpper <- c(t(sapply(predsr, function(x) x[['fit']]))) + 
		c(t(sapply(predsr, function(x) x[['se.fit']])))
fits1$CILower <- exp(fits1$LinkFit - (qnorm(0.975) * fits1$LinkSE))
fits1$CIUpper <- exp(fits1$LinkFit + (qnorm(0.975) * fits1$LinkSE))
fitsw <- fits1[1 : (dim(fits1)[1] / 2), ]
fitsg <- fits1[(1 + (dim(fits1)[1] / 2)) : dim(fits1)[1], ]

## Get model-predicted differences pre-post
diff3 <- data.frame(FFG = levels(ffg1$FFG)) 
	diff3$LinkDiff = unlist(lapply(ffgmod, function(x) fixef(x[['nb6']])$cond[2]))
	diff3$DensityDiff = round(exp(diff3$LinkDiff))
	## The above three lines are in-progress. Try getting SEs using link='response' 
fitdiff <- fitsg$Density - fitsw$Density
	names(fitdiff) <- fitsw$FFG
fitperc <- round(fitdiff / fitsw$Density, 4)

## Plot panel graph of Whiting's and our data
plotpred <- function(){
	par(mfrow = c(2, 1), mar = c(3, 5, 0.2, 1), cex = 1)
	for(i in 1:2){
		mydat <- if(i == 1){fitsw} else{fitsg}
		plot(c(1, 6), c(min(mydat$SELower), max(mydat$SEUpper)), xlab = '', ylab = '', axes = FALSE, 
			type = 'n')
		if(i == 1){axis(1, at = 1:6, padj = 1, mgp = c(3, 0.2, 0), labels = FALSE)
		} else{axis(1, at = 1:6, padj = 1, mgp = c(3, 0.2, 0), labels = c('Shredder', 'Scraper/\nGrazer', 
			'Collector\nFilterer', 'Collector\nGatherer', 'Generalist', 'Predator'))
		}
		axis(2, las = 2)
		mtext(side = 2, bquote('Invertebrate density ('*m^-2*')'), line = 3.5)
		box(bty = 'l')
		points(mydat$Density, pch = 16, cex = 1.5)
		with(mydat, arrows(x0 = 1:6, y0 = SELower, y1 = SEUpper, code = 3, angle = 90, length = 0.05))
		myleg <- ifelse(i == 1, '2011', '2016')
		legend('topright', legend = myleg, bty = 'n')
	}
}
plotTypes(plotpred, 'ModelDensities', 'Figures')

## Plot graph of modeled absolute and % difference
plotpreddiff <- function(){
	par(mfrow = c(2, 1), mar = c(3, 5, 0.2, 1), cex = 1)
	plot(c(1, 6), c(min(fitsg$Density - fitsw$Density), max(fitsg$Density - fitsw$Density)), 
		xlab = '', ylab = '', axes = FALSE, type = 'n')
		axis(1, at = 1:6, padj = 1, mgp = c(3, 0.2, 0), labels = FALSE)
		axis(2, las = 2)
		box(bty = 'l')
		mtext(side = 2, bquote('Difference in density ('*m^-2*')'), line = 3.5)
		abline(h = 0, lty = 2)
		points(1:6, fitsg$Density - fitsw$Density, pch = 16, cex = 1.5)
	plot(c(1, 6), c(-1, .25), xlab = '', ylab = '', axes = FALSE, type = 'n')
	axis(1, at = 1:6, padj = 1, mgp = c(3, 0.2, 0), labels = c('Shredder', 'Scraper/\nGrazer', 
		'Collector\nFilterer', 'Collector\nGatherer', 'Generalist', 'Predator'))
	axis(2, las = 2, at = seq(-1, 1, 0.25), labels = paste0(seq(-1, 1, 0.25) * 100, '%'))
	box(bty = 'l')
	mtext(side = 2, 'Relative difference', line = 3.7)
	abline(h = 0, lty = 2)
	points(1:6, fitperc, pch = 16, cex = 1.5)
}
plotTypes(plotpreddiff, 'ModelDifferences', 'Figures')

## Plot hybrid graph of the above, with model densities and % difference
pchs <- c(17, 16, 18)
cols = c(2, 4, 1)
legs = c('2011', '2016', '2016-2011')
plotpredplusdiff <- function(){
	par(mfrow = c(3, 1), mar = c(3, 5, 0.2, 1), cex = 1)
	for(i in 1:3){
		if(i < 3){
			mydat <- if(i == 1){fitsw} else {fitsg}
			plot(c(1, 6), c(min(mydat$SELower), max(mydat$SEUpper)), xlab = '', ylab = '', axes = FALSE, 
				type = 'n')
			axis(1, at = 1:6, padj = 1, mgp = c(3, 0.2, 0), labels = FALSE)
			mtext(side = 2, bquote('Invertebrate density ('*m^-2*')'), line = 3.5)
			points(mydat$Density, pch = pchs[i], col = cols[i], cex = 1.5)
			with(mydat, arrows(x0 = 1:6, y0 = SELower, y1 = SEUpper, code = 3, angle = 90, length = 0.05))
		} else {
			plot(c(1, 6), c(-1, .25), xlab = '', ylab = '', axes = FALSE, type = 'n')
			axis(1, at = 1:6, padj = 1, mgp = c(3, 0.2, 0), labels = c('Shredder', 'Scraper/\nGrazer', 
				'Collector\nFilterer', 'Collector\nGatherer', 'Generalist', 'Predator'))	
			mtext(side = 2, 'Relative difference', line = 3.5)
			abline(h = 0, lty = 2)
			points(1:6, fitperc, pch = pchs[i], col = cols[i], cex = 1.8)
			#barplot(fitperc)
		}
		axis(2, las = 2)
		box(bty = 'l')
		legend('topleft', legend = paste(LETTERS[i], legs[i], sep = '.  '), bty = 'n')
	}
}
plotTypes(plotpredplusdiff, 'ModelDensitiesAndDifferences', 'Figures', filetype = c('pdf', 'png'))



##### Ordination analysis #####

## Put Whiting and our benthic data in a single dataframe
ord1 <- rbind(dat4$Whiting, dat4$Benthic)
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
	## 90% by row and 232% by columns. Not horrendous.

## Remove rare species
ord3 <- delRare(ord2)
	delspp <- colnames(ord2)[colnames(ord2) %in% colnames(ord3) == FALSE]
	## Removed 3 taxa (DIXL, LIBE, TINL).
cv3 <- cv(ord3)
	## Dropped cv on columns to 219%.

## Relativize by species max
ord4 <- rel(ord3)
cv4 <- cv(ord4)
	## Dropped cv on columns to 56%, but increased cv on rows to 98%.
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
	## 80% of total variation is explained by the 2 axis (66% on Axis 1). Pretty good.

## Build a matrix of grouping variables
unqrow <- match(unique(rownum), rownum)
env1 <- ord1[unqrow, c('Study', 'Season')]
	rownames(env1) <- ord1[unqrow, 'BarcodeID']

## See how groupings load on the ordination
envpts <- envfit(pts2D, env1)
	## 25% for pre-post and 6% for season. Neither is great.

## Plot groupings on ordination
plotord <- function(){
	par(mfrow = c(1, 1), mar = c(4, 4, 0.1, 0.1), cex = 1)
	plot(c(-1.5, 1.5), c(-1.5, 1.5), type = 'n', xlab = 'Axis 1', ylab = 'Axis 2', axes = FALSE)
	axis(1)
	axis(2, las = 2)
	box(bty = 'l')
	points(pts2D, col = as.numeric(env1$Study) + 1, pch = as.numeric(env1$Season) + 14)
	ordiellipse(pts2D, env1$Study, kind = 'sd', conf = 0.95, col = c(2, 3))
	text(spp2D, rownames(spp2D), col = 4, cex = 0.6)
	legend('topright', legend = c('2011', '2016', levels(env1$Season), '95% CI'), 
		col = c(2, 3, rep(1, 4), 8), pch = c(rep(15, 2), 15:18, 1), 
		pt.cex = c(1.5, 1.5, rep(1, 4), 1.5), bty = 'n')
}
plotTypes(plotord, 'Ordination', 'Figures', height = 6.5)
	## Not much there.

## See how taxa load on the ordination
envspp <- envfit(pts2D, ord4)
envspp1 <- data.frame(round(cbind(envspp$vectors[[1]], envspp$vectors[[2]], envspp$vectors[[4]]), 4))
	colnames(envspp1) <- c('Axis1', 'Axis2', 'R2', 'p')
	envspp1 <- envspp1[order(-envspp1$R2),]

## Limit to only taxa with R2 >= 25%
envspp2 <- envspp1[envspp1$R2 >= 0.25,]
	envspp2$FFG <- spp[match(rownames(envspp2), spp$SpeciesID), 'FFG']
	## TABL, HYDE, CHIL, PETL, CHIM, PLAN, and LEPT are most strongly loading on the ordination.
	## No real theme to these in terms of FFG, although most are loading in the 
		## direction of Whiting's samples.


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
	richmod$Richness$r3 <- glmmTMB(Richness ~ Study + (1 | Season), data = rich2)
	richmod$Shannon$r1 <- glmmTMB(Shannon ~ 1, data = rich2)
	richmod$Shannon$r2 <- glmmTMB(Shannon ~ 1 + (1 | Season), data = rich2)
	richmod$Shannon$r3 <- glmmTMB(Shannon ~ Study + (1 | Season), data = rich2)
	richmod$Simpson$r1 <- glmmTMB(Simpson ~ 1, data = rich2)
	richmod$Simpson$r2 <- glmmTMB(Simpson ~ 1 + (1 | Season), data = rich2)
	richmod$Simpson$r3 <- glmmTMB(Simpson ~ Study + (1 | Season), data = rich2)	
	names(richmod) <- colnames(rich2)[1 : 3]
richmodAIC <- AICtable(richmod)
	## Adding a season random effect helps all models, and adding pre-post really helps.

## Get best model for Count, Size, and Biomass
bestmod2 <- bestmod(richmod)

## Get parameters to predict
predparm2 <- data.frame(Study = levels(rich2$Study), Season = 'Generic')
	
## Get model-predicted values and confidence intervals
preds2 <- lapply(bestmod2, predict, newdata = predparm2, type = 'response', se.fit = TRUE, 
	allow.new.levels = TRUE)
fits2 <- list()
for(i in 1 : length(bestmod2)){
	fits2[[i]] <- data.frame(Study = levels(samp1$Study))
	fits2[[i]]$Fit <- round(preds2[[i]]$fit, 2)
	fits2[[i]]$SELower <- round(preds2[[i]]$fit - preds2[[i]]$se.fit, 2)
	fits2[[i]]$SEUpper <- round(preds2[[i]]$fit + preds2[[i]]$se.fit, 2)
	fits2[[i]]$CILower <- round(preds2[[i]]$fit - (qnorm(0.975) * preds2[[i]]$se.fit), 2)
	fits2[[i]]$CIUpper <- round(preds2[[i]]$fit + (qnorm(0.975) * preds2[[i]]$se.fit), 2)	
}
	names(fits2) <- names(bestmod2)
	## Decline in all metrics from pre to post.

## Plot panel graph of overall change in Richness and Shannon's and Simpson's diversity
plotpred2 <- function(){
	par(mfrow = c(3, 1), mar = c(3, 5, 0.2, 1), cex = 1)
	for(i in 1:3){
		mydat <- fits2[[i]]
		plot(c(0.5, 2.5), c(min(mydat$SELower), max(mydat$SEUpper)), xlab = '', ylab = '', axes = FALSE, 
			type = 'n')
		if(i != 3){axis(1, at = 1:2, labels = FALSE)
		} else{axis(1, at = 1:2, labels = c('2011', '2016'))
		}
		axis(2, las = 2)
		if(i == 1){mtext(side = 2, 'Richness (taxa per sample)', line = 3.5)}
		if(i == 2){mtext(side = 2, "Shannon's diversity (H)", line = 3.5)}
		if(i == 3){mtext(side = 2, "Simpson's diversity (D)", line = 3.5)}
		box(bty = 'l')
		points(mydat$Fit, pch = 16, cex = 1.5)
		with(mydat, arrows(x0 = 1:2, y0 = SELower, y1 = SEUpper, code = 3, angle = 90, length = 0.05))
	}
}
plotTypes(plotpred2, 'ModelRichnessDiversity', 'Figures', width = 4)
	

### NEXT STEPS: 
	## Could be interesting to see if predators have gotten bigger, for instance (are large Corydalus and Odonates invulnerable to dace predation?)
	## Figure out how to get confidence intervals on density differences and relative differences (probably a Jeff task).