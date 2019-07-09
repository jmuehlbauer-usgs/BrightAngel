##### Bright Angel trophic cascade analysis #####
## Last updated 8 July 2019 by J.D. Muehlbauer


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
dat <- lapply(paste0(gitdat, gitfiles, '.csv'), read.csv)
		names(dat) <- c('Drift', 'Benthic', 'Whiting')
spp <- read.csv(paste0(gitdat, 'SpeciesList.csv'))


##### Clean up data #####

## Convert dates usable format
dat <- lapply(dat, transform, Date = as.Date(Date))

## Give Whiting Data a barcode name for consistency
dat$Whiting$BarcodeID <- paste(dat$Whiting$Date, dat$Whiting$SampleID)

## Add functional feeding groups
dat <- lapply(dat, transform, FFG = spp[match(SpeciesID, spp$SpeciesID), 'FFG'])

## Add season
dat <- lapply(dat, transform, Season = factor(ifelse(month(Date) == 11, 'November', 
	ifelse(month(Date) == 1, 'January', ifelse(month(Date) == 6, 'June', 'September'))),
	levels = c('November', 'January', 'June', 'September')))
## Add sample area and raw counts to Whiting data
	## Note: Whiting used a 0.086 m^2 Hess, but combined 2 samples per each of his 6 "samples"
dat$Whiting$Area <- 0.086 * 2
dat$Whiting$CountTotal <- dat$Whiting$Density * dat$Whiting$Area
colnames(dat$Benthic)[which(colnames(dat$Benthic) == 'SampleArea')] <- 'Area'

## Keep only columns of interest
dat0 <- lapply(dat, function(x){
	mycols <- c('BarcodeID', 'Date', 'Season', ifelse('Density' %in% colnames(x), 'Area', 'Volume'),
		'SpeciesID', 'FFG', 'CountTotal')
	x[, mycols]
	})
dat0 <- lapply(dat0, setNames, nm = c('BarcodeID', 'Date', 'Season', 'Unit', 'SpeciesID', 'FFG', 'Count'))

##### Clean up taxa #####

## Get taxa list of present taxa
taxa <- rbind(spp[spp$SpeciesID %in% dat0$Drift$SpeciesID,c('SpeciesID', 'Description')], 
	spp[spp$SpeciesID %in% dat0$Benthic$SpeciesID, c('SpeciesID', 'Description')], 
	spp[spp$SpeciesID %in% dat0$Whiting$SpeciesID,c('SpeciesID', 'Description')])
	taxa <- taxa[match(unique(taxa$SpeciesID), taxa$SpeciesID),]
	taxa <- taxa[order(taxa$Description),]

## Convert all taxa in different life stages to same Species ID (e.g., CHIL, CHIP, CHIA all become CHIA)
dat1 <- lapply(dat0, transform, SpeciesID = as.character(SpeciesID))
origt <- c('MCYA', 'CERA', 'CERP', 'CHIA', 'CHIP', 'CULP', 'WIEA', 'SIMA', 'SIMP', 'BASP', 'BAET', 'LEPA', 
	'CAPA', 'TRIA', 'TRIP', 'HYSP', 'HYDA')
newt <- c('MCYL', 'CERL', 'CERL', 'CHIL', 'CHIL', 'CULL', 'WIEL', 'SIML', 'SIML', 'BAEL', 'BAEL', 'LEPL', 
	'CAPL', 'TRIL', 'TRIL', 'HYDE', 'HYDL')
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
	## Whiting probably didn't separate Empidid or Tabanid taxa

## Combine congenerics
origt1 <- c('ELML', 'ELOA', 'PROB', 'HEMR', 'WIEL', 'SILV', 'TABS', 'DICL', 'DRAL', 'DAML', 'COEN', 
	'CAPL', 'TRIL', 'HYOS', 'HYDL', 'HYLA', 'POLY', 'RHCL')
newt1 <- c('MCYL', 'MCYL', 'CERL', 'EMPL', 'EMPL', 'TABL', 'TABL', 'TIPL', 'LIBE', 'ARGI', 'ARGI', 
	'CAPN', 'HYDE', 'HYDE', 'LETR', 'LETR', 'POLL', 'RHYL')
dat2 <- lapply(dat1, transform, SpeciesID = ifelse(SpeciesID %in% origt1, 
	newt1[match(SpeciesID, origt1)], SpeciesID))

## Cut oddballs, delete 0 count rows
dat2 <- lapply(dat2, function(x) x[!(x[, 'SpeciesID'] %in% c('NEMA', 'LYMN', 'RHAG', 'CLAM')),])
dat2 <- lapply(dat2, function(x) x[x[, 'Count'] != 0,])

## Combine rows of same taxa from different life stages
combinefx <- function(oldlist, groupby){
	tlist <- list()
	for(i in 1:length(oldlist)){
		t1 <- oldlist[[i]]
			t1$BarGrp <- paste(t1$BarcodeID, t1[,groupby])
		t2 <- t1[match(unique(t1$BarGrp), t1$BarGrp),]
		t3 <- aggregate(t1$Count, by = list(t1$BarGrp), sum)
		t2$Count <- t3[match(t2$BarGrp, t3[, 1]), 2]
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

## Sort by SampleID then SpeciesID, drop old factor levels
dat4 <- lapply(dat4, function(x) x[order(x[, 'BarcodeID'], x[, 'SpeciesID']),])
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


##### Group specimens by FFG #####

## Create new list, with samples summed by FFG
ffg <- combinefx(dat4, 'FFG')
	ffg <- lapply(ffg, subset, select = -SpeciesID)

## Get relative counts
ffg <- lapply(ffg, function(x){
	totals <- tapply(x[, 'Count'], x[, 'BarcodeID'], sum)
	x[, 'Relative'] <- round(x[, 'Count'] / totals[match(x[, 'BarcodeID'], names(totals))], 4)
	return(x)
	})

## Get mean, standard error, and relative mean by FFG, by season
ffgmn <- lapply(ffg, function(x) round(tapply(x[, 'Count'], list(x[, 'FFG'], x[, 'Season']), mean)))
ffgsem <- lapply(ffg, function(x) round(tapply(x[, 'Count'], list(x[, 'FFG'], x[, 'Season']), 
	function(x) sd(x) / sqrt(length(x))), 4))
ffgunit <- lapply(ffg, function(x) 
	round(tapply(x[, 'Count'] / x[, 'Unit'], list(x[, 'FFG'], x[, 'Season']), mean), 4))
ffgrel <- lapply(ffg, function(x) round(tapply(x[, 'Relative'], list(x[, 'FFG'], x[, 'Season']), mean), 4))

ffgstat <- list(ffgmn, ffgsem, ffgunit, ffgrel)
	names(ffgstat) <- c('MeanCount', 'SEM', 'MeanUnit', 'MeanRelative')
	ffgstat <- lapply(ffgstat, function(x){
		lapply(x, function(y){
			y[is.na(y)] <- 0
			return(y)
		})
	})	


##### Look at differences by FFG over time #####

## Get difference in absolute density and relative abundance from Whiting to our study
diff1 <- lapply(ffgstat[c('MeanUnit', 'MeanRelative')], 
	function(x) x['Benthic'][[1]] - x['Whiting'][[1]])
diff2 <- lapply(diff1, function(x){
	data.frame(Mean = round(rowMeans(x), 4), 
		SEM = round(apply(x, 1, function(y) sd(y) / sqrt(length(y))), 4))
	})

## Plot difference in absolute density and relative abundance
plotdiff <- function(){
	par(mfrow = c(2, 1), mar = c(3, 5, 0.2, 1), cex = 1)
	for(i in 1:2){
		plot(c(1, 6), c(min(diff2[[i]]$Mean - diff2[[i]]$SEM), max(diff2[[i]]$Mean + diff2[[i]]$SEM)), 
		xlab = '', ylab = '', axes = FALSE, type = 'n')
		if(i == 1){
			axis(1, at = 1:6, padj = 1, mgp = c(3, 0.2, 0), labels = FALSE)
			axis(2, las = 2)
			mtext(side = 2, bquote('Difference in density ('~m^-2*')'), line = 3.5)
		} else{
			axis(1, at = 1:6, padj = 1, mgp = c(3, 0.2, 0), labels = c('Shredder', 'Scraper/\nGrazer', 
				'Collector\nFilterer', 'Collector\nGatherer', 'Generalist', 'Predator'))
			axis(2, las = 2, at = seq(-0.1, 0.3, 0.1), labels = paste0(seq(-0.1, 0.3, 0.1) * 100, '%'))
			mtext(side = 2, 'Difference in relative abundance', line = 3.7)
		}
	abline(h = 0, lty = 2)
	box(bty = 'l')
	points(diff2[[i]]$Mean, pch = 16, cex = 1.5)
	with(diff2[[i]], arrows(x0 = 1:6, y0 = Mean - SEM, y1 = Mean + SEM, code = 3, angle = 90, length = 0.05))
	}
}
plotTypes(plotdiff, 'ObservedDifferences', 'Figures')
	## Panels are pretty similar (which is good).
	
	
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
	ffgmod[[i]]$nb2 <- glmmTMB(Count ~ 1 + (1 | BarcodeID) + offset(log(Unit)), 
		family = 'nbinom2', data = mydat)
	ffgmod[[i]]$nb3 <- glmmTMB(Count ~ 1 + (1 | Season) + offset(log(Unit)), 
		family = 'nbinom2', data = mydat)
	ffgmod[[i]]$nb4 <- glmmTMB(Count ~ 1 + (1 | BarcodeID) + (1 | Season) + offset(log(Unit)), 
		family = 'nbinom2', data = mydat)
	## Add pre-post and FFG fixed effects
	ffgmod[[i]]$nb5 <- glmmTMB(Count ~ Study + (1 | BarcodeID) + offset(log(Unit)), 
		family = 'nbinom2', data = mydat)
	ffgmod[[i]]$nb6 <- glmmTMB(Count ~ Study + (1 | Season) + offset(log(Unit)), 
		family = 'nbinom2', data = mydat)
	ffgmod[[i]]$nb7 <- glmmTMB(Count ~ Study + (1 | BarcodeID) + (1 | Season) + offset(log(Unit)), 
		family = 'nbinom2', data = mydat)	
}
	names(ffgmod) <- levels(ffg1$FFG)
	
## Compare models using AIC
ffgmodAIC <- sapply(ffgmod, function(x) sapply(x, AIC))
	## Negative binomial distribution always the best bet (over linear/Gaussian and Poisson).
	## The choice of random effect or whether to include one is less clear.
		## However, BarcodeID doesn't always fit, and a random effect for season is sensible. Use that.
	## Including a pre-post term improves all models except for Collector-Filterers and -Gatherers.
		## Critical to the analysis though, so use it anyway.

## Get parameters to predict
tlen <- length(levels(ffg1$Study))
elen <- length(levels(ffg1$Season))
predparm <- data.frame(Study = levels(ffg1$Study), Season = 'Generic', Unit = 1)
	
## Get model-predicted values and confidence intervals
preds <- lapply(ffgmod, function(x){
	predict(x[['nb6']], predparm, type = 'link', se.fit = TRUE, allow.new.levels = TRUE)
	})
fits1 <- data.frame(FFG = rep(levels(ffg1$FFG), tlen), 
	Study = rep(levels(ffg1$Study), rep(flen, tlen)), Unit = 1)
	fits1$LinkFit <- c(t(sapply(preds, function(x) x[['fit']])))
	fits1$LinkSE <- c(t(sapply(preds, function(x) x[['se.fit']])))
fits1$Density <- round(exp(fits1$LinkFit))
fits1$CILower <- round(exp(fits1$LinkFit - (qnorm(0.975) * fits1$LinkSE)))
fits1$CIUpper <- round(exp(fits1$LinkFit + (qnorm(0.975) * fits1$LinkSE)))
fitsw <- fits1[1 : (dim(fits1)[1] / 2), ]
fitsg <- fits1[(1 + (dim(fits1)[1] / 2)) : dim(fits1)[1], ]

## Get model-predicted differences pre-post
diff1 <- data.frame(FFG = levels(ffg1$FFG)) 
	diff1$LinkDiff = unlist(lapply(ffgmod, function(x) fixef(x[['nb6']])$cond[2]))
	diff1$DensityDiff = round(exp(diff1$LinkDiff))
	## The above three lines are in-progress. 
fitdiff <- fitsg$Density - fitsw$Density
	names(fitdiff) <- fitsw$FFG
fitperc <- round(fitdiff / fitsw$Density, 4)

## Plot panel graph of Whiting's and our data
plotpred <- function(){
	par(mfrow = c(2, 1), mar = c(3, 5, 0.2, 1), cex = 1)
	for(i in 1:2){
		mydat <- if(i == 1){fitsw} else{fitsg}
		plot(c(1, 6), c(min(mydat$CILower), max(mydat$CIUpper)), xlab = '', ylab = '', axes = FALSE, 
			type = 'n')
		if(i == 1){axis(1, at = 1:6, padj = 1, mgp = c(3, 0.2, 0), labels = FALSE)
		} else{axis(1, at = 1:6, padj = 1, mgp = c(3, 0.2, 0), labels = c('Shredder', 'Scraper/\nGrazer', 
			'Collector\nFilterer', 'Collector\nGatherer', 'Generalist', 'Predator'))
		}
		axis(2, las = 2)
		mtext(side = 2, bquote('Invertebrate density ('~m^-2*')'), line = 3.5)
		box(bty = 'l')
		points(mydat$Density, pch = 16, cex = 1.5)
		with(mydat, arrows(x0 = 1:6, y0 = CILower, y1 = CIUpper, code = 3, angle = 90, length = 0.05))
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
		mtext(side = 2, bquote('Difference in density ('~m^-2*')'), line = 3.5)
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


##### Ordination analysis #####

## Put Whiting and our bnethic data in a single dataframe
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
	## Dropped cv on columns to 54%, but increased cv on rows to 97%.
	## Probably OK for analysis from this point.

## Run an NMS with stepdown
NMS(ord4, maxruns = 10000)
0
	## Scree plot is typical; suggests 2-3 dimensions would be best.
	## Going to work with 2D for ease
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
	## 24% for pre-post and 17% for season. Neither is great.

## Plot groupings on ordination
plotord <- function(){
	par(mfrow = c(1, 1), mar = c(4, 4, 0.1, 0.1), cex = 1)
	plot(c(-1.5, 1.5), c(-1.5, 1.5), type = 'n', xlab = 'Axis 1', ylab = 'Axis 2', axes = FALSE)
	axis(1)
	axis(2, las = 2)
	box(bty = 'l')
	points(pts2D, col = as.numeric(env1$Study) + 1, pch = as.numeric(env1$Season) + 14)
	#ordihull(pts2D, env1$Study, draw = 'polygon', col = c(2, 3))
	ordiellipse(pts2D, env1$Study, kind = 'sd', conf = 0.95, col = c(2, 3))
	text(spp2D, rownames(spp2D), col = 4, cex = 0.6)
	legend('topright', legend = c('2011', '2016', levels(env1$Season), '95% CI'), col = c(2, 3, rep(1, 5)), 
		pch = c(rep(15, 2), 15:18, 1), pt.cex = c(1.5, 1.5, rep(1, 4), 1.5), bty = 'n')
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
	## TABL, HYDE, CHIL, PETL, CHIM, and PLAN are most strongly loading on the ordination.
	## No real theme to these in terms of FFG, although most maybe except HYDE),
		## are loading in the direction of Whiting's samples.
		
### NEXT STEPS: 
	## Consider converting Whiting's Biomass to useable format and backcalculating to lengths to pair with our size data.
	## Could be interesting to see if predators have gotten bigger, for instance (are large Corydalus and Odonates invulnerable to dace predation?)
	## Other patterns shown in the model analysis are largely consistent with a trophic cascade where trout removal increases dace, which decrease inverts.
	## Looking at change in Biomass as a surrogate for secondary production/standind stock could also be useful.
	## Figure out how to get confidence intervals on density differences and relative differences (probably a Jeff task).