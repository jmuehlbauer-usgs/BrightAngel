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


##### Get some basic stats #####

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
	
## Get difference in relative abundance and absolute density across seasons
diff1 <- lapply(ffgstat[c('MeanUnit', 'MeanRelative')], 
	function(x) x['Benthic'][[1]] - x['Whiting'][[1]])
diff2 <- lapply(diff1, function(x){
	data.frame(Mean = round(rowMeans(x), 4), 
		SEM = round(apply(x, 1, function(y) sd(y) / sqrt(length(y))), 4))
	})

## Plot differernce in absolute density


## Plot difference in relative abundance from Whiting to our study
plotdiff <- function(){
	par(mfrow = c(2, 1), mar = c(3, 5, 0.1, 0.1))
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
plotdiff()
	## Panels are pretty similar (which is good). 
	## Probably should just use relative abundance moving forward.

##### Plot FFGs by season #####



## Make a rough barplot of relative counts
par(mfrow = c(2, 1))
barplot(ffgmn$Whiting)
barplot(ffgmn$Benthic, legend.text = TRUE)
	### Stopped here. What follows is unverified and probably needs fixing.
	### Might consider doing assumed trophic level instead of ffg.


barplot(ffgmn$Benthic)
lapply(FFG, function(x) tapply(x[, 'Count'], x[, 'SpeciesID'], sum))


## Combine data by FFG for benthics, drift, and Whiting
dffg <- tapply(d1$Concentration, d1$FFG, function(x){sum(x, na.rm = TRUE)})
bffg <- tapply(b1$Density, b1$FFG, function(x){sum(x, na.rm = TRUE)})
wffg <- tapply(w1$Density, w1$FFG, function(x){sum(x, na.rm = TRUE)})

## Combine all data streams, sort by trophic level
ffg <- as.data.frame(cbind(dffg, bffg, wffg))
	rownames(ffg)[1] <- 'Unknown'
ffg1 <- ffg[rownames(ffg) %in% c('Shredder', 'CollectorFilterer', 'CollectorGatherer', 'ScraperGrazer', 'Generalist', 'Predator'),]
ffg2 <- ffg1[c(6, 1, 2, 5, 3, 4),]
	colnames(ffg2) <- c('Drift', 'Benthic', 'Whiting')


panel <- function(){}
par(mfrow = c(2, 1), mar = c(1.5, 5.5, 0.1, 0.1), oma = c(3.2, 0, 0, 0), xpd = FALSE)


wbar <- barplot(ffg2$Whiting, axes = FALSE, xlab = '', ylab = '', names.arg = FALSE, col = 'tomato')
box(bty = 'l')
axis(1, at = wbar, labels = rep('', length(wbar)))
axis(2, las = 2)
mtext(side = 2, expression(paste('Density (# * ', m^-2, ')')), line = 4)
legend('topleft', legend = '2010-2011', bty = 'n')


bbar <- barplot(ffg2$Benthic, axes = FALSE, xlab = '', ylab = '', names.arg = FALSE, col = 'slateblue1')
box(bty = 'l')
axis(1, at = bbar, labels = c('Shredders\n', 'Collector-\nfilterers', 'Collector-\ngatherers', 'Scrapers/\nGrazers', 'Generalists\n', 'Predators\n'), padj = 0.5)
axis(2, las = 2)
mtext(side = 2, expression(paste('Density (# * ', m^-2, ')')), line = 4)
mtext(side = 1, 'Functional feeding group', line = 3.5)
legend('topleft', legend = '2016-2017', bty = 'n')




###Calculating relative densities

ffg2$BenthicRel <- round(ffg2$Benthic / sum(ffg2$Benthic), 4)
ffg2$WhitingRel <- round(ffg2$Whiting / sum(ffg2$Whiting), 4)
#Note that I've combined a couple functions onto a single line of code (rather than calculating the sum on a separate line as I suggested over the phone). The result is the same, this is just a little cleaner. I've also included the round function, which just rounds the results the the specified number of decimal places (4 in this case).

#Subtracting Whitings from ours to see if there was an increase or decrease in densities.
ffg2$BenthicRel - ffg2$WhitingRel
Reldif <- (ffg2$BenthicRel - ffg2$WhitingRel) 


###Barplotting the difference in relative densities 
Relbar <- barplot(Reldif, axes = FALSE, xlab = '', ylab = '', names.arg = FALSE, col = 'slateblue1')
box(bty = 'l')
axis(1, at = bbar, labels = c('Shredders\n', 'Collector-\nfilterers', 'Collector-\ngatherers', 'Scrapers/\nGrazers', 'Generalists\n', 'Predators\n'), padj = 0.5)
axis(2, -1:1, las = 2)
mtext(side = 2, expression(paste('Density differences')), line = 4)
mtext(side = 1, 'Functional Feeding Groups', line = 3.5)
legend('topleft', legend = '', bty = 'n')


###ANOVA test
aov(abundance~ffg)
#error code. Begin here next time. 
 

