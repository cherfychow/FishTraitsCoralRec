

#######
# Responses to fish trait diversity in coral settlement and recruitment
# Data cleanup
# Author: Cher Chow
# Updated: 26 Mar 2021
#######


# Setup -------------------------------------------------------------------

require(tidyverse) # data handling
require(readxl) # import xls files
require(rfishbase) # for retrieving functional trait data

require(ggplot2) # plotting
require(patchwork) # multipanel plotting


require(FD) # calculate dissimilarity matrices + functional diversity indices
require(geometry) # calculates convex hull/trait space for data vis
require(lme4) # for fitting GLMM
require(performance) # diagnostics

set.seed(24) # repeatability :)

# plot aesthetics
quartzFonts(Public=c('Public Sans Regular', 'Public Sans Italic', 'Public Sans Bold', 'Public Sans Bold Italic'))
par(family='Public')
looks <- theme_bw(base_size=13, base_family = 'Public Sans') + theme(panel.grid=element_blank(), axis.ticks=element_line(size=0.3))

# functions
ci.pred <- function(.) predict(., newx, type='response', re.form=NA)
source('convhullvert_function.R') # convex hull vertices for plotting function

## ----Data import'------------------------------------------------------------------------------------------------

# read in bite data from folder
files = list.files(path='/Users/cher/OneDrive - University of St Andrews/Dissertation/DataFish/UBRUV_bite', pattern='*.csv', full.names = T)
bite.sheets = lapply(files, read.csv, header=T) # read in all the xls files in the folder
bite.data <- do.call(rbind, bite.sheets)
head(bite.data, n=10)
rm(bite.sheets)
rm(files)

# read in assemblage data from folder
files = list.files(path='/Users/cher/OneDrive - University of St Andrews/Dissertation/DataFish/UBRUV_spL/', pattern='*.xls', full.names = T)
sp.sheets = lapply(files, read_xlsx, sheet=1, col_names=T, col_types=c('date',rep('text', 6),'numeric',rep('text', 4))) # read in all the xls files in the folder
ubruv.data <- do.call(rbind, sp.sheets)
head(ubruv.data)
rm(sp.sheets)
rm(files)


# ----check and trim-------------------------------------------------------------------------------------------------------------
str(ubruv.data)
summary(ubruv.data)

# do some species ID corrections per Andy Hoey
corrections <- c('Naso brevirostris',
                 'Thalassoma lutescens',
                 'Diagramma pictum',
                 'Stethojulis trilineata')
names(corrections) <- c('^Naso annulatus$', 
                        '^Halichoeres chrysus$',
                        '^Diagramma sp$',
                        '^Bodianus bimaculatus')
ubruv.data$Species <- str_replace_all(ubruv.data$Species, corrections)
length(unique(ubruv.data$Species)) # 117 species recorded
ubruv.trim <- filter(ubruv.data, str_detect(Species, 'idae|unknown', negate = TRUE)) #remove uncertain observations at Family level or higher
# 468 observations left, removed 7 rows


## ----taxon spell check, results='hide'------------------------------------------------------------------------------------------
Sp <- sort(unique(str_replace(ubruv.trim$Species, ' sp$', '')))
Sp <- unique(str_replace(Sp, ' $', '')) # remove accidental duplicates ending with spaces
# e.g. change "Scarus sp" to "Scarus" for WoRMS validation
require(worms) # package for validating species/taxa names
Sp.list <- wormsbymatchnames(Sp, marine_only=T, chunksize = 49) # feed it into worms package to check
Sp.list$ID <- seq(1, length(Sp.list$valid_name), by=1)
invalid <- Sp.list %>% filter(match_type != 'exact' | status != 'accepted') %>% select(ID, valid_name, family, genus)
invalid$before <- Sp[invalid$ID]
ubruv.validated <- ubruv.trim # create a duplicate of the trimmed dataset to replace with validated names

validate <- invalid$valid_name
names(validate) <- invalid$before # make a spell check object
ubruv.validated$SpeciesV <- str_replace_all(ubruv.validated$Species, validate) # replace misspellings

# clean up the length categories
length.convert <- c('2.5', '15', '25', '35', '45', '7.5', '55', '85') # in order as shown by sort(unique(ubruv.validated$Length))
names(length.convert) <- paste('^', sort(unique(ubruv.validated$Length)), '$', sep='')
ubruv.validated$Length <- str_replace_all(ubruv.validated$Length, length.convert)
detach(package:worms)
detach(package:plyr)


## ----species summary------------------------------------------------------------------------------------------------------------
Sp.list <- Sp.list %>% select(valid_name, family) %>% mutate(genus=word(valid_name, 1)) # let's trim the worms table to just taxonomic info, and extract genus from the valid name, not the record that matches the wrong species

Sp.count <- as_tibble(ubruv.validated) %>% group_by(SpeciesV) %>% summarise(n=sum(Count)) %>% filter(str_detect(SpeciesV, ' sp$', negate=T)) # counts per species
Sp.countL <- as_tibble(ubruv.validated) %>% group_by(SpeciesV, Length) %>% summarise(n=sum(Count)) # counts per species and length class
ubruv.validated %>% group_by(Site) %>% summarise(n=sum(Count)) # total number of fish per site
UBRUVspecies <- ubruv.validated %>% group_by(Site, SpeciesV) %>% summarise(n=sum(Count)) # counts per species grouped by site
UBRUVspeciesL <- ubruv.validated %>% group_by(Site, SpeciesV, Length) %>% summarise(n=sum(Count))
# make a separate tibble for grouping by length classes


colnames(Sp.list)[1] <- "SpeciesV" #rename the Species column in the worms dataframe
UBRUVspecies$genus <- word(UBRUVspecies$SpeciesV, 1) # extract genus names to make join easier
UBRUVspecies <- left_join(UBRUVspecies, distinct(Sp.list[,2:3]), by='genus') # add a column of family names to UBRUV datasheet

Sp.count$genus <- word(Sp.count$SpeciesV, 1) # add genus column to the species count summary tibble
Sp.count <- left_join(Sp.count, distinct(Sp.list[,2:3]), by='genus') # and add family taxonomic info too


## ----Species prevalences, echo=FALSE--------------------------------------------------------------------------------------------
ggplot(head(Sp.count[order(Sp.count$family),] %>% filter(n > 5), n=15), aes(y=Species, x=n)) +
  geom_bar(aes(fill=family), color='transparent', stat='identity') +
  labs(y='Species', x='Count') +
  scale_fill_discrete(name='Family') +
  looks + scale_x_continuous(expand=expansion(mult=c(0,.1)))

ggplot(Sp.count %>% filter(n > 5) %>% group_by(family) %>% summarise(totalF=sum(n)), aes(y=family, x=totalF)) +
  geom_bar(fill='white', color='black', stat='identity') +
  labs(y='Family', x='Count') +
  looks + scale_x_continuous(expand=expansion(mult=c(0,.1)))


## ----site composition plot, echo=FALSE------------------------------------------------------------------------------------------
# composition by counts
site.comp <- ggplot(UBRUVspecies %>% group_by(Site, family) %>% summarise(n=sum(n)) %>% filter(n > 2), aes(x=Site, y=n)) +
  geom_bar(aes(fill=family), stat='identity', color='white', size=.6) +
  labs(x=NULL, y='Count') +
  scale_fill_viridis_d(option='plasma', name='Family') +
  scale_x_discrete(labels=c('CB', 'L1', 'N3', 'R', 'SE', 'TB', 'V')) +
  looks + scale_y_continuous(expand=expansion(mult=c(0,.01))) + 
  theme(legend.position='none', axis.text.x=element_blank(), axis.ticks.x = element_blank())

# relative composition
site.compR <- ggplot(UBRUVspecies %>% group_by(Site, family) %>% summarise(n=sum(n)) %>% filter(n > 2), aes(x=Site, y=n)) +
  geom_bar(aes(fill=family), position='fill', color='white', size=.6, stat='identity') +
  labs(x='Study sites', y='Relative abundance') +
  scale_fill_viridis_d(option='plasma', name='Family') +
  scale_x_discrete(labels=c('CB', 'L1', 'N3', 'R', 'SE', 'TB', 'V')) +
  looks + scale_y_continuous(expand=expansion(mult=c(0,.01)))

(site.comp / site.compR) + plot_layout(guides='collect')


## ----Bite data summaries--------------------------------------------------------------------------------------------------------
str(bite.data)
summary(bite.data)


## ----Validate names------------------------------------------------------------------------
# check if there are any species IDs at the genus level
sum(str_detect(bite.data$Species, ' sp$')+0)
# ok, no genus IDs to replace/fix. validate the names from WoRMS.
# make the same corrections per Andy Hoey
corrections <- c('Naso brevirostris',
                 'Thalassoma lutescens',
                 'Diagramma pictum',
                 'Stethojulis trilineata')
names(corrections) <- c('^Naso annulatus$', 
                        '^Halichoeres chrysus$',
                        '^Diagramma sp$',
                        '^Bodianus bimaculatus')
bite.data$Species <- str_replace_all(bite.data$Species, corrections)
library(worms) # reload the WoRMS package
BiteSp <- sort(unique(bite.data$Species))
BiteSp <- unique(str_replace(BiteSp, ' $', '')) # remove accidental duplicates ending with spaces
# e.g. change "Scarus sp" to "Scarus" for WoRMS validation
BiteSp.list <- wormsbymatchnames(BiteSp, marine_only=T) # feed it into worms package to check
BiteSp.list$ID <- seq(1, length(BiteSp.list$valid_name), by=1)
invalid <- BiteSp.list %>% filter(match_type != 'exact' | status != 'accepted') %>% select(ID, valid_name, family, genus)
invalid$before <- BiteSp[invalid$ID]
bite.validated <- bite.data # create a duplicate of the trimmed dataset to replace with validated names

validate <- invalid$valid_name
names(validate) <- invalid$before # make a spell check object
bite.validated$Species <- str_replace_all(bite.validated$Species, validate) # replace misspellings
sum(str_detect(bite.validated$Species, invalid$before)+0) # use str_detect to see if the spell check's worked

# once everything is done, detach and remove, esp the worms package bc of conflicts
detach(package:worms)
detach(package:plyr)
rm(corrections)
rm(validate)
rm(invalid)


## ----Duration data clean--------------------------------------------------------------------------------------------------------
require(lubridate) # use lubridate to deal with time data
bite.validated$BiStart <- with(bite.validated, ms(BiStart))
bite.validated$BiEnd <- with(bite.validated, ms(BiEnd))
bite.validated$Duration <- with(bite.validated, as.numeric(as.duration(BiEnd)-as.duration(BiStart))) # Duration in seconds!
head(bite.validated)
# calculate duration of occurrences/bites and add as new column


## ----Bite + Sp combine----------------------------------------------------------------------------------------------------------
# clean up the length categories
length.convert <- c('15', '2.5', '2.5', '15', '25', '7.5') # in order as shown by sort(unique(bite.validated$Length))
names(length.convert) <- paste('^', sort(unique(bite.validated$Length)), '$', sep='')
bite.validated$Length <- str_replace_all(bite.validated$Length, length.convert)

# join that all up grouped by site and species (this group of data wrangling is NOT BY LENGTH CATEGORIES)
UBRUVbite <- bite.validated %>% group_by(Site, Species) %>% summarise(BiteTotal=sum(BiTotal), DurationTotal=sum(Duration)) %>% filter(BiteTotal > 0) # only those that bite
UBRUVbite$BiteRate <- with(UBRUVbite, (BiteTotal/DurationTotal)*60) # bites per minute
head(UBRUVbite, n=10)
# join the species and bite rate data (NOT CLASSIFIED BY LENGTH)
colnames(UBRUVspecies)[2] = 'Species' # column name align for joining
UBRUV <- full_join(UBRUVspecies, UBRUVbite, by=c('Site', 'Species')) # join it according to site and species
UBRUV$BiteRateAbund <- with(UBRUV, BiteRate/n) # bites per min per individual
UBRUV <- UBRUV %>% replace_na(list(BiteTotal=0, DurationTotal=0, BiteRate=0, BiteRateAbund=0)) # create a bite per second metric for each species

# remove the genus level IDs
UBRUV <- UBRUV %>% filter(str_detect(Species, ' sp$', negate=T))
UBRUVspecies <- UBRUVspecies %>% filter(str_detect(Species, ' sp$', negate=T))
UBRUV <- UBRUV %>% filter(str_detect(Species, ' sp$', negate=T))
head(UBRUV, n=10)
# now that we have the counts per species, we can standardise bite rate per individual


## ----Bite + Species data with Length!-------------------------------------------------------------------------------------------
# a separate tibble for bites grouped by LENGTH, species, and site.
UBRUVbiteL <- bite.validated %>% group_by(Site, Species, Length) %>% summarise(BiteTotal=sum(BiTotal), DurationTotal=sum(Duration), Bitemaxn=max(nInd)) %>% filter(BiteTotal > 0) # only those that bite
UBRUVbiteL$BiteRate <- with(UBRUVbiteL, (BiteTotal/DurationTotal)*60) # bites per minute
head(UBRUVbiteL)
# and join like we did before
colnames(UBRUVspeciesL)[2] = 'Species' # change column name to make sure things line up in joining
UBRUV.L <- full_join(UBRUVspeciesL, UBRUVbiteL, by=c('Site', 'Species', 'Length')) # join it according to site and species
UBRUV.L$BiteRateAbund <- with(UBRUV.L, BiteRate/Bitemaxn) # bites per min per individual
UBRUV.L <- UBRUV.L %>% replace_na(list(BiteTotal=0, DurationTotal=0, BiteRate=0, BiteRateAbund=0)) # create a bite per second metric for each species
head(UBRUV.L, n=10)

# make a column for length classifications that's for treating length as discrete scales
length.convert <- c('10-20', '< 5', '20-30', '5-10') # in order as shown by sort(unique(UBRUVbiteL$Length))
names(length.convert) <- paste('^', sort(unique(UBRUVbiteL$Length)), '$', sep='')
UBRUVbiteL$LengthC <- str_replace_all(UBRUVbiteL$Length, length.convert)
# for UBRUV.L
length.convert <- c('10-20', '< 5', '20-30', '30-40', '40-50', '50-60', '5-10', '80-90')
# in order as shown by sort(unique(UBRUV.L$Length))
names(length.convert) <- paste('^', sort(unique(UBRUV.L$Length)), '$', sep='')
UBRUV.L$LengthC <- str_replace_all(UBRUV.L$Length, length.convert)


## ----Summaries, include=F-------------------------------------------------------------------------------------------------------
str(UBRUV)
summary(UBRUV)


## ----Distributions, echo=FALSE--------------------------------------------------------------------------------------------------
# these don't include the zeroes from the other non-biting fish so that we don't overinflate the distribution
# first a histogram of the bite totals
plot.btotal.hist <- ggplot(UBRUVbite, aes(x=BiteTotal)) +
  geom_histogram(fill='white', color='black', binwidth=25) +
  labs(x='Bite totals', y='Frequency') +
  looks

# histogram of bite rates (not standardised by species abundance)
plot.brate.hist <- ggplot(UBRUVbite, aes(x=BiteRate)) +
  geom_histogram(fill='white', color='black', binwidth=5) +
  labs(x=bquote('Bite rate' ~ '(bites'~ 'min'^-1 ~ ')'), y='Frequency') +
  looks

# histogram and comparison of bite rate (A) across sites
plot.brateA.hist <- ggplot(UBRUV %>% filter(BiteTotal > 0), aes(x=BiteRateAbund)) +
  geom_histogram(fill='white', color='black', binwidth=2.5) +
  labs(x=bquote('Bite rate' ~ '(bites' ~ 'min'^-1 ~ 'indv'^-1 ~')'), y='Frequency') + looks

plot.btotal.hist + plot.brate.hist + plot.brateA.hist
  


## ----Bite rate relationships, echo=FALSE----------------------------------------------------------------------------------------
plot.brateA.dot <- ggplot(UBRUV %>% filter(BiteTotal > 0), aes(x=Site, y=BiteRateAbund)) +
  geom_boxplot() +
  geom_jitter(shape=1, size=2) +
  labs(x='Sites', y=bquote('Bite rate' ~ '(bites' ~ 'min'^-1 ~ 'indv'^-1 ~')')) + looks +
  theme(axis.text.x=element_text(angle=30, hjust=1))

plot.bite_dur <- ggplot(UBRUVbite %>% filter(BiteTotal > 0), aes(x=DurationTotal, y=BiteTotal)) +
  geom_point(aes(color=Site), alpha=0.6, size=3) +
  labs(x='Observation time (s)', y='Total bites') + looks +
  scale_x_continuous(limits=c(min(UBRUVbite$DurationTotal), max(UBRUVbite$DurationTotal))) +
  theme(legend.position='none')

plot.bite_durLog <- ggplot(UBRUVbite %>% filter(BiteTotal > 0), aes(x=log(DurationTotal), y=log(BiteTotal))) +
  geom_smooth(method='lm', size=0.2, color='black') +
  geom_point(aes(color=Site), alpha=0.6, size=3) +
  labs(x='log(Observation time)', y='log(Total bites)') + looks +
  scale_x_continuous(limits=c(min(log(UBRUVbite$DurationTotal)), max(log(UBRUVbite$DurationTotal))))

(plot.bite_dur + plot.bite_durLog)
plot.brateA.dot


## ----Bite rates across taxa, echo=F---------------------------------------------------------------------------------------------

plot.bifamily <- ggplot(UBRUV %>% filter(BiteTotal > 0), aes(x=family, y=BiteRate)) + 
  # *60 for per min rates, more meaningful numbers
  geom_boxplot() +
  labs(x='Family', y=bquote('Bite rate' ~ '(bites' ~ 'min'^-1 ~')')) + looks +
  theme(axis.text.x=element_text(angle=30, hjust=1))

plot.bifamilyL <- ggplot(UBRUVbite %>% filter(BiteTotal > 0) %>% left_join(y=traits[,c('Species', 'ForageMode')], by='Species'), aes(x=ForageMode, y=BiteRate)) +
  geom_jitter(shape=21, size=2.6, width=0.25, aes(fill=Site)) +
  geom_boxplot(fill='transparent') +
  scale_fill_viridis_d(begin=0.1, end=0.95, alpha=0.8) +
  scale_y_log10() +
  labs(x='Foraging mode', y=NULL) + looks +
  theme(axis.text.x=element_text(angle=30, hjust=1), axis.text.y=element_blank(), legend.position='none')

plot.bitroph <- ggplot(UBRUVbite %>% filter(BiteTotal > 0) %>% left_join(y=traits[,c('Species', 'TrophicGroup')], by='Species'), aes(x=TrophicGroup, y=BiteRate)) +
  geom_jitter(shape=21, size=2.6, width=0.2, aes(fill=Site)) +
  geom_boxplot(fill='transparent') +
  scale_fill_viridis_d(begin=0.1, end=0.95, alpha=0.8) +
  scale_y_log10() +
  labs(x='Trophic group', y=bquote('Bite rate' ~ '(bites' ~ 'min'^-1 ~')')) + looks +
  theme(axis.text.x=element_text(angle=30, hjust=1))

plot.bifamily + plot.bifamilyL
(plot.bitroph | plot.bifamilyL) + plot_layout(guides='collect')




## ----Biting across sites in families--------------------------------------------------------------------------------

plot.bite.family <- ggplot(UBRUV %>% filter(BiteTotal > 0) %>% filter(str_detect(family, 'Acanthuridae|Chaetodontidae|Labridae|Scaridae')), aes(x=Species, y=BiteRateAbund*60)) +
  geom_jitter(shape=21, size=2) +
  labs(x='Species', y=bquote('Bite rate' ~ '(bites' ~ 'min'^-1 ~ 'indv'^-1 ~')')) + looks +
  theme(axis.text.x=element_text(angle=30, hjust=1)) +
  facet_grid(~family, scales='free')
  
plot.bite.family

plot.bite.site <- ggplot(UBRUV %>% filter(BiteTotal > 0) %>% filter(str_detect(family, 'Acanthuridae|Chaetodontidae|Labridae|Scaridae|Siganidae')), aes(x=Site, y=BiteRate)) +
  geom_bar(stat='identity', fill='white', color='black') +
  labs(x='Study sites', y=bquote('Bite rate' ~ '(bites' ~ 'min'^-1 ~ ')')) + looks +
  theme(axis.text.x=element_text(angle=30, hjust=1)) +
  facet_grid(~family)

plot.bite.site


## ----Differences in biting family across sites----------------------------------------------------------------------------------
plot.site.bite <- ggplot(UBRUV %>% filter(BiteTotal > 0) %>% group_by(Site, family) %>% summarise(BiteRate=sum(BiteRate)), aes(x=family, y=BiteRate)) +
  geom_bar(color='black', fill='white', stat='identity', width=0.8) +
  labs(x='Family', y=bquote('Total bite rate' ~ '(bites' ~ 'min'^-1~')')) + looks +
  theme(strip.background=element_rect(fill='white', size=0.5), strip.text=element_text(face='bold'), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  facet_wrap(~Site, ncol=4) + scale_y_continuous(expand=expansion(mult=c(0,.05))) +
  scale_x_discrete(labels=c('Ac', 'Ba', 'Ch', 'La', 'Mu', 'Pca', 'Pce', 'Sc', 'Si'))

plot.site.bite


## ----Bite-Length relationship, echo=F-------------------------------------------------------------------------------------------
ggplot(UBRUV.L %>% filter(BiteTotal > 0), aes(x=as.numeric(Length), y=BiteRateAbund)) +
  geom_point(shape=21, size=2, alpha=0.7) +
  labs(x='Length class', y=bquote('Bite rate' ~ '(bites' ~ 'min'^-1 ~ n[max]^-1 ~')')) + looks


## ----Pull Fishbase data---------------------------------------------------------------------------------------------------------
fishbase.sp <- validate_names(sort(unique(UBRUV$Species))) # expect it to remove 5 entries because they're at the genus level
sp.traits <- species(fishbase.sp) # make functional trait data frame for the species list
ecol.traits <- ecology(fishbase.sp) # table of ecological traits


## ----Trait selection------------------------------------------------------------------------------------------------------------
trait.data <- with(ecol.traits, data.frame(Species, FeedingType, FoodTroph, FoodSeTroph))

# deal with the average length column now
length.convert <- c('2.5', '15', '25', '35', '45', '7.5', '55', '85')
names(length.convert) <- sort(unique(Sp.countL$Length))
Sp.countL$Length <- str_replace_all(Sp.countL$Length, length.convert)
Sp.countLength <- Sp.countL %>% group_by(SpeciesV, Length) %>% mutate(sumL=as.numeric(Length)*n) %>% # transitional column
  ungroup() %>% group_by(SpeciesV) %>% summarise(totaln=sum(n), totalL=sum(sumL)) %>% # compress to one row per species
  mutate(meanL=totalL/totaln) # calculate the meanL
colnames(Sp.count)[1] <- 'Species'
colnames(Sp.countLength)[1] <- 'Species'
trait.data <- left_join(trait.data, Sp.countLength %>% select(Species, meanL), by='Species')
food.data <- fooditems(fishbase.sp) # load food data from fishbase
morph.data <- morphometrics(fishbase.sp)
head.data <- morph.data %>% group_by(Species) %>% summarise(POLprop=mean(POL/SL)*100)
# for preorbital length? possible proxy for bite size (but very flawed)

# refine the traits manually outside of R
write_csv(trait.data, na='NA', col_names=T, path='./trait_data.csv')
rm(trait.data)
rm(food.data)
rm(morph.data)
rm(head.data)

# reimport it after supplementing the trait data manually and assemble into table
trait.data <- read_csv('trait_data_manual_v03.csv', col_names=T) # overwrite trait.data with the manual table
trait.range <- read_csv('trait_range_CC.csv', col_names=T)
# merge them, clean up a little. not a whole lot is needed though.
traits <- trait.data %>% select(-Source) %>% left_join(.,trait.range[,4:6], by='Species')
head(traits)
# make sure R reads trait data types properly
traits <- as.data.frame(traits)
str(traits)
traits[,-1] <- traits[,-1] %>% mutate(across(where(is.character), as.factor)) # turn characters into factors
traits$Schooling <- as.ordered(traits$Schooling)
traits$Range <- as.ordered(traits$Range)
str(traits) # confirm!


## ----Pairwise trait check-------------------------------------------------------------------------------------------------------
pairs(traits[,-1], labels=colnames(traits)[-1], col=rgb(0.2,0.2,0.2, 0.5))


## ----Trait space all------------------------------------------------------------------------------------------------------------
# I think I will drop trophicSE as a trait.
# traits$FoodSeTroph <- NULL

# create the abundance matrix for functional evenness and divergence weighting
# species are columns, each site + total = rows
abund <- Sp.count %>% filter(str_detect(Species, ' sp$', negate=T)) %>% select(Species, n) %>% spread(Species, n) # row 1 = all sites
abund$Site <- 'All'
abund <- rbind(UBRUVspecies %>% select(Site, Species, n) %>% spread(Species, n), abund)
# because not all species occur in each site, there are lots of NAs to be filled with 0s
NAlist <- as.list(rep('0',dim(abund)[2]-1)) # list has to have column names for every item = 0
names(NAlist) <- colnames(abund)[2:(dim(abund)[2])]
abund <- abund %>% replace_na(NAlist) # annnd replace NAs
rownames(traits) <- traits$Species # dbFD function requires species in the row names

# abundance data frame is in characters????
abund <- abund %>% mutate(across(where(is.character), as.numeric))
abund <- as.data.frame(abund)
rownames(abund) <- abund$Site
abund$Site <- NULL


## ----Trait space sites, results='hide'------------------------------------------------------------------------------------------
# running the full lizard island traits together is very computationally taxing
# so I'll split all the sites to run the FD package separately
abun <- as.list(rep('0',7))
for (i in 1:7) {
  abun[[i]] <- as.data.frame(UBRUV %>% filter(Site == unique(UBRUV$Site)[i]) %>% ungroup() %>% select(Species, n) %>% spread(Species, n))
}
names(abun) <- unique(UBRUV$Site)
abun[[8]] <- Sp.count %>% filter(str_detect(Species, ' sp$', negate=T)) %>% select(Species, n) %>% spread(Species, n)
names(abun)[8] <- 'all' # add one item for all sites

# also make separate trait tables for each site
traits.sites <- as.list(rep('0', 7))
for (i in 1:7) {
  site.sp <- UBRUV %>% filter(Site == unique(UBRUV$Site)[i]) %>% ungroup() %>% select(Species)
  # filter the tibble by each site and extract the species
  traits.sites[[i]] <- as.data.frame(left_join(site.sp, as_tibble(traits), by='Species'))
  # Use the site specific species list to match with the traits
  rownames(traits.sites[[i]]) <- traits.sites[[i]]$Species # row names bc FD requires it
  traits.sites[[i]][,1] <- NULL # Remove the Species column now that things are joined
}
names(traits.sites) <- unique(UBRUV$Site) # name list items by site
rm(site.sp)
traits.sites[[8]] <- traits[,-1]
names(traits.sites)[8] <- 'all' # all sites again

# create the distance matrix for each site
# Gower distance, Podani method for ordered factors.
trait.dis <- as.list(rep('0', 8))
discheck <- rep('0', 8)
for (i in 1:8) {
  # calculate gower dissimilarities on standardised continuous traits
  # and apply Cailliez correction to make it Euclidean
  trait.dis[[i]] <- cailliez(gowdis(traits.sites[[i]], ord='podani'))
  discheck[i] <- is.euclid(trait.dis[[i]])
}

# check the distance matrices. The FD package will check this too, but good to check before PCoA
discheck
rm(discheck)

# run PCoA and NMDS for all sites, using both because PCoA is returning a lot of dimensions.
PCoA <- as.list(rep('0', 8))
NMDS <- as.list(rep('0', 8))
for (i in 1:8) {
  PCoA[[i]] <-  pcoa(trait.dis[[i]], rn=colnames(trait.dis[[i]])) #PCoA with Cailliez corrections to negative eigenvalues
  NMDS[[i]] <- vegan::metaMDS(trait.dis[[i]], distance='gower', k=4, wascores=T)
}


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------
# Scree plots to check PCoAs, but only 
scree <- as.list(rep('0', 8))
for (i in 1:8) {
  scree[[i]] <- ggplot(data=PCoA[[i]]$values[1:10,], aes(x=1:10, y=Relative_eig/sum(PCoA[[i]]$values$Relative_eig))) +
  geom_line() +
  geom_point(shape=21, fill='white', size=3) + looks +
  labs(x=NULL, y=NULL) +
  scale_x_continuous(breaks=c(1:10)) +
  theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))
}

scree[[1]] + scree[[2]] + scree[[3]] + scree[[4]] + scree[[5]] + scree[[6]] + scree[[7]] + scree[[8]]



## ----Scaled vs observed distance, separate spaces-------------------------------------------------------------------------------
P.dis2 <- as.list(rep(0, 8))
P.dis4 <- as.list(rep(0, 8))
P.stress <- as.list(rep(0, 8))
df_dis <- as.list(rep(0, 8))
P.lm2 <- as.list(rep(0, 8))
P.lm4 <- as.list(rep(0, 8))
for (i in 1:8) {
  P.dis2[[i]] <- cailliez(gowdis(PCoA[[i]]$vectors[,1:2], ord='podani'))
  P.dis4[[i]] <- cailliez(gowdis(PCoA[[i]]$vectors[,3:4], ord='podani'))
  # look at the euclidean distances between species in two 2-D trait spaces
  
  # compress the distance matrices in P.dis into a dataframe
  df_dis[[i]] <- data.frame(Axis12=as.matrix(P.dis2[[i]])[,1], Axis34=as.matrix(P.dis4[[i]])[,1], Observed=as.matrix(trait.dis[[i]])[,1])
  for (j in 2:dim(PCoA[[i]]$vectors)[1]) {
    df_dis[[i]] <- rbind(df_dis[[i]], data.frame(Axis12=as.matrix(P.dis2[[i]])[,j], Axis34=as.matrix(P.dis4[[i]])[,j], Observed=as.matrix(trait.dis[[i]])[,j]))
    df_dis[[i]] <- df_dis[[i]] %>% filter(Observed > 0, Axis12 > 0, Axis34 > 0)
  }
  
  # conduct Mantel tests to get correlations between the observed and scaled dissimilarity matrices
  # Pearson correlation for parametric correlation
  P.lm2[[i]] <- mantel(trait.dis[[i]], P.dis2[[i]], method='pearson', permutations=500)
  P.lm4[[i]] <- mantel(trait.dis[[i]], P.dis4[[i]], method='pearson', permutations=500)
}

    P.stress[[9]] <- ggplot(data=data.frame(x=1:9, y=1:9), aes(x, y)) +
      annotate(geom='text', x=2, y=1:8, hjust=0, label=names(traits.sites), size=3) + looks +
      theme(axis.ticks=element_blank(), axis.text=element_blank()) + labs(x=NULL, y=NULL) + scale_y_continuous(limits=c(0,9))
  # plot all the observed vs scaled dissimilarities with the R2 from the lm we did for every site.
for (i in 1:8) {
  P.stress[[i]] <- ggplot(df_dis[[i]]) + 
    geom_point(aes(x=Observed, y=Axis12), shape=21, alpha=0.1) +
    geom_point(aes(x=Observed, y=Axis34), shape=21, alpha=0.1, color='blue') +
    looks + labs(x=NULL, y=NULL, title=names(traits.sites)[i]) +
    theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))
  P.stress[[9]] <- P.stress[[9]] + 
    annotate(geom='text', x=5, y=i, hjust=0, label=round(P.lm2[[i]]$statistic, 3), size=3) +
    annotate(geom='text', x=7, y=i, hjust=1, label=round(P.lm4[[i]]$statistic, 3), color='blue', size=3)
}

# scaled vs observed dissimilarities in the four dimensional trait space    
P.dis <- as.list(rep(0, 8))
P.stressall <- as.list(rep(0, 8))
P.lm <- as.list(rep(0,8))
df_all <- as.list(rep(0,8))
for (i in 1:8) {
  P.dis[[i]] <- cailliez(gowdis(PCoA[[i]]$vectors[,1:4], ord='podani'))
  # look at the euclidean distances between species in 1st 4 dimensions
  
  # compress the distance matrices in P.dis into a dataframe
  df_all[[i]] <- data.frame(Scaled=as.matrix(P.dis[[i]])[,1], Observed=as.matrix(trait.dis[[i]])[,1])
  for (j in 2:dim(PCoA[[i]]$vectors)[1]) {
    df_all[[i]] <- rbind(df_all[[i]], data.frame(Scaled=as.matrix(P.dis[[i]])[,j], Observed=as.matrix(trait.dis[[i]])[,j]))
    df_all[[i]] <- df_all[[i]] %>% filter(Observed > 0, Scaled > 0)
  }
  # perform a mantel test to test the correlation of the observed vs 4D dissimilarity matrices
  P.lm[[i]] <- mantel(trait.dis[[i]], P.dis[[i]], method='pearson', permutations=500)
}
  
P.stressall[[9]] <- ggplot(data=data.frame(x=1:9, y=1:9), aes(x, y)) +
      annotate(geom='text', x=2, y=1:8, hjust=0, label=names(traits.sites), size=3) + looks +
      theme(axis.ticks=element_blank(), axis.text=element_blank()) + labs(x=NULL, y=NULL) + scale_y_continuous(limits=c(0,9))

for (i in 1:8) {
  # plot all the observed vs scaled dissimilarities in the matrices for every site
  P.stressall[[i]] <- ggplot(df_all[[i]]) + 
    geom_point(aes(x=Observed, y=Scaled), shape=21, alpha=0.1) +
    looks + labs(x=NULL, y=NULL, title=names(traits.sites)[i]) +
    theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))
  P.stressall[[9]] <- P.stressall[[9]] + 
    annotate(geom='text', x=5, y=i, hjust=1, label=round(P.lm[[i]]$statistic, 3), size=3)
}

# print plots

P.stressall[[1]] + P.stressall[[2]] + P.stressall[[3]] + P.stressall[[4]] + 
  P.stressall[[5]] + P.stressall[[6]] + P.stressall[[7]] + P.stressall[[8]] + P.stressall[[9]]
(scree[[8]] + labs(x='Dimensions', y='Relative eigenvalue')) | (P.stressall[[8]] + labs(x='Original dissimilarities', y='PCoA dissimilarities')) & labs(title=NULL)


## ----correlation between eigenvalues + dissimilarities--------------------------------------------------------------------------
EigMant <- data.frame(SumEig=sum(PCoA[[1]]$values$Relative_eig[1:4]), Mantel=P.lm[[1]]$statistic, sig=P.lm[[1]]$signif)
for (i in 2:8) {
  EigMant <- rbind(EigMant, data.frame(SumEig=sum(PCoA[[i]]$values$Relative_eig[1:4]), Mantel=P.lm[[i]]$statistic, sig=P.lm[[i]]$signif))
}
EigMant.lm <- lm(data=EigMant, Mantel~SumEig)
summary(EigMant.lm)

library(ggrepel)
ggplot(EigMant, aes(x=SumEig, y=Mantel)) +
  geom_smooth(method='lm', se=TRUE, color='black', size=0.5) +
  geom_point(shape=21, size=2.5) + 
  geom_text_repel(label=rownames(abund)) + looks


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------
P.pt1 <- as.list(rep(0, 8))
P.pt2 <- as.list(rep(0, 8))
for (i in 1:8) {
P.pt1[[i]] <- ggplot(as.data.frame(PCoA[[i]]$vectors), aes(x=Axis.1, y=Axis.2)) +
    geom_point(shape=21) + looks +
    labs(x=NULL, y=NULL, title=names(traits.sites)[i]) +
    theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))
P.pt2[[i]] <- ggplot(as.data.frame(PCoA[[i]]$vectors), aes(x=Axis.3, y=Axis.4)) +
    geom_point(shape=21) + looks +
    labs(x=NULL, y=NULL, title=names(traits.sites)[i]) +
    theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))
}
(P.pt1[[1]] | P.pt1[[2]] | P.pt1[[3]] | P.pt1[[4]]) / (P.pt2[[1]] | P.pt2[[2]] | P.pt2[[3]] | P.pt2[[4]])
(P.pt1[[5]] | P.pt1[[6]] | P.pt1[[7]] | P.pt1[[8]])/(P.pt2[[5]] | P.pt2[[6]] | P.pt2[[7]] | P.pt2[[8]])


## ----NMDS stress, echo=FALSE----------------------------------------------------------------------------------------------------
par(mfrow=c(2,2))
for (i in 1:8) {
  stressplot(NMDS[[i]], p.col='grey', l.col='black')
  title(main=names(traits.sites)[i])
} # check stress and fit
par(mfrow=c(1,1))
par(mfrow=c(2,4))
for (i in 1:4) {
  with(as.data.frame(NMDS[[i]]$points), plot(x=MDS1, y=MDS2))
  title(main=names(traits.sites)[i], sub=paste('stress =', round(NMDS[[i]]$stress, 3), sep=' '))
}
for(i in 1:4) {
  with(as.data.frame(NMDS[[i]]$points), plot(x=MDS3, y=MDS4))
}
for (i in 5:8) {
  with(as.data.frame(NMDS[[i]]$points), plot(x=MDS1, y=MDS2))
  title(main=names(traits.sites)[i], sub=paste('stress =', round(NMDS[[i]]$stress, 3), sep=' '))
}
for(i in 5:8) {
  with(as.data.frame(NMDS[[i]]$points), plot(x=MDS3, y=MDS4))
}
par(mfrow=c(1,1))


## ----FD calculation-------------------------------------------------------------------------------------------------------------
# calculate FD indices for each site with respect to each site's functional trait space
FD <- as.list(rep(0,8))
for (i in 1:8) {
  FD[[i]] <- dbFD(trait.dis[[i]], abun[[i]], w.abun=T, calc.FRic=T, calc.FDiv=T, m=4, calc.CWM=F, calc.FGR=F, print.pco=T)
  # not standardising FRic because there is no global FRic measure by running each site separately
}
  names(FD) <- names(abun)[[8]]
  
# calculate FD indices for each site with respect to the global functional trait space
# because this is with respect to all study sites, we can standardise FRichness
FD_global <- dbFD(trait.dis[[8]], abund, w.abun=T, calc.FRic=T, calc.FDiv=T, m=4, calc.CWM=F, calc.FGR=F, stand.FRic=T, print.pco=T)


## ----Indices to dataframe-------------------------------------------------------------------------------------------------------
# extract the FD indices from the FD_global object, create a dataframe of them
predictors <- data.frame(Site=names(FD_global$FEve[1:7]),FEve=FD_global$FEve[1:7], FDiv=FD_global$FDiv[1:7], FRic=FD_global$FRic[1:7])
# only 1:7 to ignore the global measures
rownames(predictors) <- names(FD_global$FEve[1:7])


## ----Trait space plots----------------------------------------------------------------------------------------------------------
t.space <- as.list(rep(0,4))
point <- as.data.frame(FD_global$x.axes)[,1:4]
point$Species <- rownames(traits)

source('convhullvert_function.R') # wrote my own function to make plottable convex hull vertices
hull.v <- convhull.vert(point[,1:2])
hull.v2 <- convhull.vert(point[,3:4])

# set up the base layers of the global functional trait space
global1 <- ggplot() + looks +
  geom_polygon(data=hull.v, aes(x=A1, y=A2), color='grey', fill='white') +
  labs(x='Axis 1', y='Axis 2')
global2 <- ggplot() + looks +
  geom_polygon(data=hull.v2, aes(x=A3, y=A4), color='grey', fill='white') +
  labs(x='Axis 3', y='Axis 4')

# palette for sites
palette <- viridis::viridis(n=7, end=0.95, begin=0.1, alpha=0.7)
require(ggrepel)
for (i in 1:7) {
  # Now that the base layers are done, I'll plot each site's points and hulls on
  # first, make a dummy points dataframe only with the species in that site
  s.point <- UBRUVspecies %>% filter(Site == unique(UBRUVspecies$Site)[i]) %>% ungroup() %>% select(Species, n)
  s.point$n <- s.point$n/site.total$site.total[i]
  s.point <- left_join(s.point, point, by='Species')
  s.point <- s.point[order(s.point$n, decreasing=T),order(colnames(s.point))]
  # make abbreviated species names
  shortsp <- as.data.frame(matrix(unlist(strsplit(s.point$Species[1:3], " ")), ncol=2, byrow=TRUE))
  s.point$ShortSp <- ''
  s.point$ShortSp[1:3] <- paste(abbreviate(shortsp[,1], 1, strict = TRUE), ". ", shortsp[,2], sep="")
  
  # now calculate the convex hull for that site only
  # global convex hull/FRic
    s.hull.v <- convhull.vert(s.point[,1:2])
    s.hull.v2 <- convhull.vert(s.point[,3:4])
  
  s.space <- as.list(c(0,0))
  
  s.space[[1]] <- global1 +
    geom_polygon(data=s.hull.v, aes(x=A1, y=A2), fill=palette[i], alpha=0.25) +
    geom_point(data=s.point, aes(x=A1, y=A2, size=n), fill=palette[i], color='#666666', shape=21) +
    # geom_text_repel(data=s.point, aes(label=ShortSp, x=A1, y=A2), force=2,
    #                 size=3.6, fontface='italic', min.segment.length = 0, point.size=5,
    #                 max.overlaps = Inf, box.padding = 1, force_pull=0.1) +
    scale_size_continuous(range=c(1.5,8), guide=F, limits=c(0,0.6)) +
    labs(title=unique(UBRUVspecies$Site)[i]) +
    if (i != 4) {theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position='none')} else {theme()}
  s.space[[2]] <- global2 +
    geom_polygon(data=s.hull.v2, aes(x=A3, y=A4), fill=palette[i], alpha=0.25) +
    geom_point(data=s.point, aes(x=A3, y=A4, size=n), fill=palette[i], color='#666666', shape=21) +
    # geom_text_repel(data=s.point, aes(label=ShortSp, x=A3, y=A4), force=2,
    #                 size=3.6, fontface='italic', min.segment.length = 0, point.size=5,
    #                 max.overlaps = Inf, box.padding = 1, force_pull=0.1) +
    scale_size_continuous(range=c(1.5,8), limits=c(0,0.6))  + 
  if (i != 4) {theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position='none')} else {theme()}
  
  t.space[[i]] <- s.space[[1]] / s.space[[2]] # stack them with patchwork and store
}

rm(s.point)
point$n <- Sp.count$n/624 # put in global abundances to visualise all sites
global1 <- ggplot() + looks +
  geom_polygon(data=hull.v, aes(x=A1, y=A2), fill='#EEEEEE') +
  geom_point(data=point, aes(x=A1, y=A2, size=n), alpha=0.3, shape=21) +
  labs(x='Axis 1', y='Axis 2', title='All sites') +
  theme(plot.title=element_text(face='bold', hjust=0.5, size=13)) +
  scale_size_continuous(range=c(1.5,8), limits=c(0,0.5), guide=F)
  
global2 <- ggplot() + looks +
  geom_polygon(data=hull.v2, aes(x=A3, y=A4), fill='#EEEEEE') +
  geom_point(data=point, aes(x=A3, y=A4, size=n), alpha=0.3, shape=21) +
  labs(x='Axis 3', y='Axis 4') +
  scale_size_continuous(range=c(1.5,8), limits=c(0,0.5)) + theme(legend.title=element_blank())

FD_bar <- ggplot(predictors %>% select(Site, FRic, FEve, FDiv) %>% pivot_longer(c(FRic, FEve, FDiv))) +
  geom_bar(aes(y=value, x=Site, alpha=name, fill=Site), stat='identity', position='dodge', color='black') +
  labs(y='Index measure', x='Study sites') +
  scale_fill_viridis_d(begin=0.1, end=0.95, guide=F) + looks + scale_alpha_discrete(name=NULL, range=c(0.2, 1)) +
  scale_y_continuous(expand=expansion(mult=c(.0,.05)), limits=c(0,1)) +
  scale_x_discrete(labels=c('CB', 'L1', 'N3', 'R', 'SE', 'TB', 'V'))

A_bar <- ggplot(predictors %>% select(Site, HerbProp, BiterProp) %>% pivot_longer(c(HerbProp, BiterProp))) +
  geom_bar(aes(y=value, x=Site, alpha=name, fill=Site), stat='identity', position='dodge', color='black') +
  labs(y='Relative abundance', x='Study sites') +
  scale_fill_viridis_d(begin=0.1, end=0.95, guide=F) + looks + scale_alpha_discrete(name=NULL, range=c(0.5, 1), labels=c('Benthic biters', 'Herbivores')) +
  scale_y_continuous(expand=expansion(mult=c(.0,.05))) +
  scale_x_discrete(labels=c('CB', 'L1', 'N3', 'R', 'SE', 'TB', 'V'))

FD_bar1 <- ggplot(data=predictors, aes(x=Site)) +
  geom_bar(aes(y=FRic, fill=Site), color='black', stat='identity') + looks +
  scale_fill_viridis_d(begin=0.1, end=0.95, alpha=0.7) +
  labs(x=NULL, y=NULL, title='FRic') +
  scale_y_continuous(expand=expansion(mult=c(.0,.05))) +
  theme(legend.position='none', axis.text.x=element_blank(), axis.ticks.x=element_blank())
FD_bar2 <- ggplot(data=predictors, aes(x=Site)) +
  geom_bar(aes(y=FEve, fill=Site), color='black', stat='identity') + looks +
  scale_fill_viridis_d(begin=0.1, end=0.95, alpha=0.7) +
  labs(x=NULL, y='Index measure', title='FEve') +
  scale_y_continuous(expand=expansion(mult=c(.0,.05))) +
  theme(legend.position='none', axis.text.x=element_blank(), axis.ticks.x=element_blank())
FD_bar3 <- ggplot(data=predictors, aes(x=Site)) +
  geom_bar(aes(y=FDiv, fill=Site), color='black', stat='identity') + looks +
  scale_fill_viridis_d(begin=0.1, end=0.95, guide=F, alpha=0.7) +
  labs(x='Study sites', y=NULL, title='FDiv') +
  scale_y_continuous(expand=expansion(mult=c(.0,.05))) +
  scale_x_discrete(labels=c('CB', 'L1', 'N3', 'R', 'SE', 'TB', 'V'))

((global1 / global2) | t.space[[1]]) + plot_layout(guides='collect')
(t.space[[2]] | t.space[[3]]) + plot_layout(guides='collect')
(t.space[[4]] | t.space[[5]]) + plot_layout(guides='collect')
(t.space[[6]] | t.space[[7]]) + plot_layout(guides='collect')
(FD_bar1 + FD_bar2 + FD_bar3) + plot_layout(guides='collect') & theme(axis.text.x=element_blank())

( (global1/global2) | t.space[[1]] | t.space[[2]] | t.space[[3]]) + plot_layout(guides='collect')
(t.space[[4]] | t.space[[5]] | t.space[[6]] | t.space[[7]] ) + plot_layout(guides='collect')
(t.space[[7]] | (global1/global2))  + plot_layout(guides='collect')

## ----Global trait space by trophic group----------------------------------------------------------------------------------------
trait.points <- left_join(traits, point, by='Species')
Troph.2 <- as.list(rep(0,6))
Troph.4 <- as.list(rep(0,6))
for (i in c(2,4,5,6,7,8)) { # draw clusters around the different trophic groups
  Troph.2[[i]] <- convhull.vert(trait.points %>% filter(TrophicGroup == sort(unique(traits$TrophicGroup))[i]) %>% select(A1,A2))
  Troph.4[[i]] <- convhull.vert(trait.points %>% filter(TrophicGroup == sort(unique(traits$TrophicGroup))[i]) %>% select(A3,A4))
}
pal <- viridis::plasma(n=8, end=0.95, begin=0.1, alpha=0.3)
palc <- viridis::plasma(n=8, end=0.95, begin=0.1)
globalgroup1 <- ggplot() + looks +
  geom_polygon(data=Troph.2[[2]], aes(x=A1, y=A2), fill=pal[2], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.2[[4]], aes(x=A1, y=A2), fill=pal[4], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.2[[5]], aes(x=A1, y=A2), fill=pal[5], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.2[[6]], aes(x=A1, y=A2), fill=pal[6], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.2[[7]], aes(x=A1, y=A2), fill=pal[7], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.2[[8]], aes(x=A1, y=A2), fill=pal[8], color='#55555599', size=0.3) +
  geom_point(data=trait.points, aes(x=A1, y=A2, fill=TrophicGroup), shape=21, size=3) +
  labs(x='Axis 1', y='Axis 2', title='Trophic group') + 
  theme(legend.position='none') +
  scale_fill_viridis_d(option='plasma', begin=0.1, end=0.95)
globalgroup2 <- ggplot() + looks +
  geom_polygon(data=Troph.4[[2]], aes(x=A3, y=A4), fill=pal[2], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.4[[4]], aes(x=A3, y=A4), fill=pal[4], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.4[[5]], aes(x=A3, y=A4), fill=pal[5], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.4[[6]], aes(x=A3, y=A4), fill=pal[6], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.4[[7]], aes(x=A3, y=A4), fill=pal[7], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.4[[8]], aes(x=A3, y=A4), fill=pal[8], color='#55555599', size=0.3) +
  geom_point(data=trait.points, aes(x=A3, y=A4, fill=TrophicGroup), shape=21, size=3) +
  labs(x='Axis 3', y='Axis 4') + 
  scale_fill_viridis_d(option='plasma', begin=0.1, end=0.95) +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  theme(legend.title=element_blank())

globalgroup1 / globalgroup2 + plot_layout(guides='collect')


## ----global trait space by foraging---------------------------------------------------------------------------------------------
# Trophic Group
global_group1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_point(aes(fill=TrophicGroup), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_viridis_d(option='plasma', begin=0.1, end=0.95) +
  theme(legend.position='none') +
  labs(title='Trophic Group')
global_group2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_point(aes(fill=TrophicGroup), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_viridis_d(option='plasma', begin=0.1, end=0.95) +
  theme(legend.position='bottom', legend.title=element_blank())

# hulls for everyone
hulls1 <- as.list(rep(0,4))
hulls2 <- as.list(rep(0,4))
for (j in c(5,6,7)) {
  hulls1[[j]] <- data.frame(level=NA, A1=NA, A2=NA)
  hulls2[[j]] <- data.frame(level=NA, A3=NA, A4=NA)
  names(hulls1)[j] <- colnames(traits)[j]
  names(hulls2)[j] <- colnames(traits)[j]
  for (i in 1:length(unique(traits[,j]))) {
    hull1 <- convhull.vert(as_tibble(trait.points) %>% filter(.data[[colnames(trait.points)[j]]] == unique(as.vector(traits[,j]))[i]) %>% select(A1,A2))
    hull2 <- convhull.vert(as_tibble(trait.points) %>% filter(.data[[colnames(trait.points)[j]]] == unique(as.vector(traits[,j]))[i]) %>% select(A3,A4))
    hull1$level <- unique(traits[,j])[i]
    hull2$level <- unique(traits[,j])[i]
    hulls1[[j]] <- rbind(hulls1[[j]], hull1)
    hulls2[[j]] <- rbind(hulls2[[j]], hull2)
  }
}
for (i in c(1,4,6,8,9,10,11)) {
  hull1 <- convhull.vert(as_tibble(trait.points) %>% filter(ForageMode == sort(unique(traits$ForageMode))[i]) %>% select(A1,A2))
  hull2 <- convhull.vert(as_tibble(trait.points) %>% filter(ForageMode == sort(unique(traits$ForageMode))[i]) %>% select(A3,A4))
  hull1$level <- sort(unique(traits$ForageMode))[i]
  hull2$level <- sort(unique(traits$ForageMode))[i]
  hulls1[[4]] <- rbind(hulls1[[4]], hull1)
  hulls2[[4]] <- rbind(hulls2[[4]], hull2)
}


# Forage Mode
global_forage1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_polygon(data=hulls1[[4]], aes(x=A1, y=A2, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=ForageMode), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_viridis_d() +
  theme(legend.position='none') +
  labs(title='Foraging Mode')
global_forage2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_polygon(data=hulls2[[4]], aes(x=A3, y=A4, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=ForageMode), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_viridis_d() +
  theme(legend.title=element_blank())

# Trophic level
global_troph1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_point(aes(fill=FoodTroph), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_gradient(low='black', high='white') +
  theme(legend.position='none') +
  labs(title='Trophic Level')
global_troph2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_point(aes(fill=FoodTroph), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_gradient(low='black', high='white') +
  theme(legend.title=element_blank())

# Water column
global_col1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_polygon(data=hulls1[[5]], aes(x=A1, y=A2, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=ColumnFeed), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_grey() +
  theme(legend.position='none') +
  labs(title='Water Column Feeding')
global_col2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_polygon(data=hulls2[[5]], aes(x=A3, y=A4, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  scale_fill_grey() +
  geom_point(aes(fill=ColumnFeed), alpha=0.7, size=3, shape=21) + looks +
  theme(legend.position='bottom', legend.title=element_blank())

# Range
global_r1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_polygon(data=hulls1[[7]], aes(x=A1, y=A2, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=Range), alpha=0.7, size=3, shape=21) + looks +
  theme(legend.position='none') +
  labs(title='Range') + scale_fill_viridis_d(option='plasma')
global_r2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_polygon(data=hulls2[[7]], aes(x=A3, y=A4, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=Range), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_viridis_d(option='plasma') +
  theme(legend.position='bottom', legend.title=element_blank())

# Schooling
global_s1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_polygon(data=hulls1[[6]], aes(x=A1, y=A2, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=Schooling), alpha=0.7, size=3, shape=21) + looks +
  theme(legend.position='none') + scale_fill_viridis_d() +
  labs(title='Schooling')
global_s2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_polygon(data=hulls2[[6]], aes(x=A3, y=A4, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=Schooling), alpha=0.7, size=3, shape=21) + looks + scale_fill_viridis_d() +
  theme(legend.position='bottom', legend.title=element_blank())

((globalgroup1 / globalgroup2) | (global_forage1 / global_forage2) | (global_troph1 / global_troph2)) + plot_layout(guides='collect')
((global_col1 / global_col2) | (global_r1 / global_r2) | (global_s1 / global_s2))


## ----Calculate species proportion-----------------------------------------------------------------------------------------------

site.total <- data.frame(site.total=rowSums(abund), Site=rownames(abund))

# herbivore abundance
site.herb <- data.frame(Site=predictors$Site, HerbProp=rep(0,7))
for (i in 1:7) {
  herb <- as.data.frame(traits.sites[[i]])
  herb$Species <- rownames(herb)
  herb <-  herb %>% filter(TrophicGroup=='herbivore') %>% ungroup() %>% select(Species)
  herb$Site <- names(traits.sites)[i]
  herb <- left_join(herb, UBRUVspecies, by=c('Site', 'Species')) %>% ungroup() %>% select(Site, Species, n)
  site.herb$HerbProp[i] <- sum(herb$n)/site.total$site.total[i]
}

# biter abundance
site.biter <- data.frame(Site=predictors$Site, HerbProp=rep(0,7))
for (i in 1:7) {
  biter <- as.data.frame(traits.sites[[i]])
  biter$Species <- rownames(biter)
  biter <-  biter %>% filter(ColumnFeed == 'benthic') %>% 
    filter(str_detect(TrophicGroup, 'detritovore|invertivore|omnivore')) %>% ungroup() %>% select(Species)
  biter$Site <- names(traits.sites)[i]
  biter <- left_join(biter, UBRUVspecies, by=c('Site', 'Species')) %>% ungroup() %>% select(Site, Species, n)
  site.biter$BiterProp[i] <- sum(biter$n)/site.total$site.total[i]
}
# add it to our predictors dataframe
predictors$HerbProp <- site.herb$HerbProp
predictors$BiterProp <- site.biter$BiterProp
# predictors$BiProp <- bi.abun$BiProp


# ## ----Taxonomic diversity--------------------------------------------------------------------------------------------------------
# sp.taxa <- worms::wormsbynames(unique(UBRUV$Species), match=T, verbose=T)
# col <- c(2,15:20)
# sp.taxa <- sp.taxa[col]
# rownames(sp.taxa) <- sp.taxa$name
# sp.taxa[,1] <- NULL
# taxdiv <- taxondive(comm=abund, dis=taxa2dist(sp.taxa))
# predictors$TaxDiv <- scale(taxdiv$D[1:7], center=F, scale=T)


## ---- Species contribution calculation loop, eval=FALSE--------------------------------------------------------------------------
# *****I'm doing this species contribution calculation with reference to global Lizard Island trait space
# 1. for each site, remove one species from abundance + traits table
# 2. recalculate FD_global for that species at that site
# 3. dataframe: row for species, column for site, column for contribution

sp.contributions <- data.frame(Site=NA, Species=NA, SpC=NA) # the final dataframe for all sites + species
# calculate an UNWEIGHTED FDiv to use in calculating deltaFDiv
FD <- dbFD(trait.dis[[8]], abund, calc.FRic=F, calc.FDiv=T, m=4, calc.CWM=F, calc.FGR=F, w.abun=F)
FD <- with(FD, c(FDiv, FEve, FRic))
# the global space doesn't change so put them in before the loop
for (i in 1:7) { # for each site
  # data frame for the site's species contributions
  SC.tbl <- as.data.frame(matrix(nrow=nrow(traits.sites[[i]]), ncol=5))
  colnames(SC.tbl) <- c('Site', 'Species', 'DivC', 'EveC', 'RicC')
  SC.tbl$Site <- predictors$Site[i] # put in the site name
  SC.tbl$Species <- rownames(traits.sites[[i]])

  # for the jth species in site i
  for (j in 1:nrow(traits.sites[[i]])) {
    # reconstruct abundance table without species j in site i
    SCabund <- Sp.count %>% filter(str_detect(Species, ' sp$', negate=T)) %>% select(Species)
    dummy <- UBRUVspecies %>% select(Site, Species, n) %>% filter(Site == predictors$Site[i])
    dummy[j,3] <- 0 # species j abundance = 0
    SCabund <- full_join(SCabund, dummy, by='Species')
    SCabund <- UBRUVspecies %>% select(Site, Species, n) %>% filter(Site != predictors$Site[i]) %>% bind_rows(.,SCabund)
    SCabund <- SCabund[,order(colnames(SCabund))] # reorder columns so that it corresponds with study site order
    # instead of using the Sp.count for total sp counts, we're summing it to account for the removal of species j
    SCabund <- SCabund %>% spread(Species, n)
    SCabund <- SCabund[-8,]
    # because not all species occur in each site, there are lots of NAs to be filled with 0s
    NAlist <- as.list(rep('0',dim(SCabund)[2]-1)) # list has to have column names for every item = 0
    names(NAlist) <- colnames(SCabund)[2:(dim(SCabund)[2])]
    SCabund <- SCabund %>% replace_na(NAlist) # annnd replace NAs
    SCabund <- SCabund %>% mutate(across(where(is.character), as.numeric))
    SCabund <- as.data.frame(SCabund)
    SCabund <- rbind(SCabund, colSums(as.data.frame(SCabund[,-1]), na.rm=T))
    rownames(SCabund) <- rownames(abund) # rownames for FD package to run
    SCabund$Site <- NULL

    # recalculate FD without this species
    SC_FD <- dbFD(trait.dis[[8]], SCabund, calc.FRic=T, calc.FDiv=T, m=4, calc.CWM=F, calc.FGR=F, w.abun=F)
    # fill in the site's dataframe with this species' contribution
    SC.tbl[j,3] <- sqrt(((FDiv[i] - SC_FD$FDiv[i]) / FDiv[i])^2) # make it a ratio
    
  }
  sp.contributions <- rbind(sp.contributions, SC.tbl) # add the site's species contributions to the final dataframe
}
sp.contributions <- sp.contributions[-1,] # remove that filler row of NAs
# remove the temporary/dummy objects used in the loop
rm(dummy)
rm(SCabund)
rm(SC.tbl)
rm(NAlist)


## ----Feeding rate calculation, eval=FALSE---------------------------------------------------------------------------------------
## FeedRate <- UBRUV.L %>% ungroup() %>% group_by(Site, Species, Length) %>% select(Site, Species, Length, BiteRate) %>%  mutate(LBi=as.numeric(Length)*BiteRate)
## FeedRate <- FeedRate %>% ungroup() %>% group_by(Site, Species) %>% filter(BiteRate > 0) %>% summarise(sumLBi=sum(LBi)) %>% left_join(., y=sp.contributions)
## predictors <- FeedRate %>% group_by(Site, Species) %>% mutate(SpLBi=sumLBi*SpC) %>% ungroup() %>%  group_by(Site) %>%  summarise(FeedRate=sum(SpLBi)) %>% full_join(., y=predictors, by='Site')
## View(predictors) # yay


## ----Abundance proportion modifier, eval=FALSE----------------------------------------------------------------------------------
## sp.a.prop <- data.frame(Site=NA, Species=NA, SpProp=NA)
## 
## for (i in 1:7) { # for each site i
##   # same setup as before. make an empty dummy df for each site
##   # this df will be used to fill in the final dataframe
##   SA.tbl <- as.data.frame(matrix(nrow=nrow(traits.sites[[i]]), ncol=3))
##   colnames(SA.tbl) <- c('Site', 'Species', 'SpProp') # name the empty df
##   SA.tbl$Site <- predictors$Site[i] # fill the site row
##   SA.tbl$Species <- rownames(traits.sites[[i]]) # fill in the species column for the site
##   site.total <- rowSums(abund[i,]) # calculate the total abundance at site for convenience
##   site.abun <- UBRUVspecies %>% filter(Site==predictors$Site[i])
##   # calculate the species abundance
##   for (j in 1:nrow(SA.tbl)) {
##     SA.tbl[j,3] <- site.abun$n[j]/site.total
##   }
##   sp.a.prop <- rbind(sp.a.prop, SA.tbl)
## }
## sp.a.prop <- sp.a.prop[-1,]
## rm(SA.tbl)
## rm(site.total)
## rm(site.abun)


## ----Option 2 feed rate, eval=FALSE---------------------------------------------------------------------------------------------
## # Calculate the feeding rate with this modifier
## FeedRate <- UBRUV.L %>% ungroup() %>% group_by(Site, Species, Length) %>% select(Site, Species, Length, BiteRate) %>%  mutate(LBi=as.numeric(Length)*BiteRate)
## FeedRate <- FeedRate %>% ungroup() %>% group_by(Site, Species) %>% filter(BiteRate > 0) %>% summarise(sumLBi=sum(LBi))
## predictors$FeedRate2 <- FeedRate %>% left_join(.,sp.a.prop, by=c('Site','Species')) %>% mutate(SpLBi=sumLBi*SpProp) %>% group_by(Site) %>% summarise(FeedRate2=sum(SpLBi)) %>% ungroup() %>% select(FeedRate2) %>% as_vector()


## ----Sp contribution by influence-----------------------------------------------------------------------------------------------
influence <- read.csv('trait_data_influence.csv', header=T)
str(influence)
BiteSp <- unique(UBRUVbite$Species)
influence <- influence %>% filter(str_detect(Species, str_c(BiteSp, collapse='|')))
inf.pca <- princomp(influence[,-1])
loadings(inf.pca)
summary(inf.pca)
screeplot(inf.pca, type='lines') # yes 

sp.inf <-  data.frame(Species=influence$Species, SpInf=inf.pca$scores[,1])
sp.inf$SpInf <- (sp.inf$SpInf + (-range(sp.inf$SpInf)[1]) + 0.1)/2

# Now use it to calculate bite rate
FeedRate3 <- left_join(FeedRate, sp.inf, by='Species') %>% mutate(SpFeed=0.01*SpInf*sumLBi) %>% group_by(Site) %>% summarise(FeedRate3=sum(SpFeed))
predictors$FeedRate3 <- FeedRate3$FeedRate3
predictors$FeedRate3 <- scale(predictors$FeedRate3, scale=T, center=F)

feedF <- feed %>% group_by(Site,TrophicGroup, ForageMode) %>% summarise(Feed=sum(SpFeed), Infl=mean(SpInf))
feedF <- factor(feedF$ForageMode, levels=feed$ForageMode[order(feedF$Infl)])
feedT <- feed %>% group_by(Site,TrophicGroup) %>% summarise(Feed=sum(SpFeed), Infl=mean(SpInf))
feedT <- factor(feedT$TrophicGroup, levels=feedT$TrophicGroup[order(feedT$Infl)])
feedforage <- ggplot(feedF, aes(x=Site, fill=Feed, y=ForageMode)) +
  geom_tile() + looks + 
  scale_fill_gradient(low='#EEEEEE', high='black', name=NULL) + 
  scale_x_discrete(labels=c('CB', 'L1', 'N3', 'R', 'SE', 'TB', 'V')) + 
  labs(x='Site', y='Foraging mode')
feedtroph <- ggplot(feedT, aes(x=Site, fill=Feed, y=TrophicGroup)) +
  geom_tile() + looks + 
  scale_fill_gradient(low='#EEEEEE', high='black', name='Feeding rate', limits=c(0,50)) + 
  scale_x_discrete(labels=c('CB', 'L1', 'N3', 'R', 'SE', 'TB', 'V')) + 
  labs(x='Site', y='Trophic group')
feedtotal <- ggplot(predictors, aes(x=Site, y=FeedRate3)) +
  geom_bar(stat='identity', fill='white', color='black', width=0.8) + looks + 
  scale_x_discrete(labels=c('CB', 'L1', 'N3', 'R', 'SE', 'TB', 'V')) + 
  labs(y='Feeding rate (cm-bites/min)', x='Site') +
  scale_y_continuous(expand=expansion(mult=c(0,.02)), breaks=seq(0,max(predictors$FeedRate3), by=10))

feed_total <- ggplot(feed, aes(x=Site, y=SpFeed)) + 
  geom_boxplot(width=0.75) + looks + 
  scale_x_discrete(labels=c('CB', 'L1', 'N3', 'R', 'SE', 'TB', 'V')) + 
  labs(x='Study site', y='Feeding rate (cm-bites/min)')
(feedtroph + feedforage + feed_total) & theme(legend.position='bottom')
(feedtroph + feedforage) + plot_layout(guides='collect')

## ----Model data structure prep--------------------------------------------------------------------------------------------------
# the following code is extracted from the settlement data exploration .rmd file
settlement <- read_xlsx('../DataCoral/settlement_2019-20_LIRS.xlsx', sheet='Sheet1', col_names=T)
head(settlement)
str(settlement)
# select just the 7 study sites and their spat data from 2019
sett.trim <- settlement %>% group_by(year, site) %>% filter(str_detect(site, 'corner_beach|lagoon_1|resort|southeast|turtle_beach|vickies|north_reef_3'), year=='2019-20')
# create a str_replace renaming object so that the site names are consistent with the predictors data
site.rename <- sort(predictors$Site)
names(site.rename) <- sort(unique(sett.trim$site))
sett.trim$site <- str_replace_all(sett.trim$site, site.rename) # replace

# now consolidate it so that every row is one settlement tile
settlement <- sett.trim %>% group_by(site, tile_number) %>% summarise(Spat=sum(total))
colnames(settlement)[1] <- 'Site'
SpatData <- full_join(settlement, predictors, by='Site') %>% select(-tile_number) %>% as.data.frame()


## ----Boxplot + Cleveland dotplots, echo=FALSE-----------------------------------------------------------------------------------
boxplots <- as.list(rep(0,4))
  boxplots[[1]] <- ggplot(SpatData[seq(1,42,by=6),], aes(x=1, y=FeedRate)) +
    geom_boxplot(width=0.2) + 
    geom_jitter(shape=21, size=2, width=0.05, alpha=0.7) + looks +
    scale_x_continuous(limits=c(0.75,1.25)) +
    scale_y_log10() +
    labs(title='FeedRate', x=NULL, y='log(FeedRate)') + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())
  
  boxplots[[2]] <- ggplot(SpatData, aes(x=1, y=Spat)) +
    geom_boxplot(width=0.2) + 
    geom_jitter(shape=21, size=2, width=0.05, alpha=0.7) + looks +
    scale_x_continuous(limits=c(0.75,1.25)) +
    labs(title='Spat', x=NULL, y='Spat') + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())
  
  boxplots[[3]] <- ggplot(SpatData[seq(1,42,by=6),] %>% pivot_longer(cols=c(FRic, FEve, FDiv)), aes(x=name, y=value)) +
    geom_boxplot() + 
    geom_jitter(shape=21, size=2, alpha=0.7) + looks +
    labs(title='FD indices', x='FD indices', y='Value')
  
  boxplots[[4]] <- ggplot(SpatData, aes(x=1, y=FeedRate3)) +
    geom_boxplot(width=0.2) + 
    geom_jitter(shape=21, size=2, width=0.05, alpha=0.7) + looks +
    scale_x_continuous(limits=c(0.75,1.25)) +
    labs(title='FeedRate3', x=NULL, y='Feed Rate 3') + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())

cleve.plots <- as.list(rep(0,5))
  cleve.plots[[1]] <- ggplot(SpatData, aes(y=1:dim(SpatData)[1], x=Spat)) +
    geom_point(shape=21, size=2) + looks +
    labs(y=NULL, x=NULL, title='Spat')
  cleve.plots[[2]] <- ggplot(SpatData[seq(1,42,by=6),], aes(y=1:7, x=FeedRate3)) +
    geom_point(shape=21, size=2) + looks +
    labs(y=NULL, x=NULL, title='Feeding Rate')
  cleve.plots[[3]] <- ggplot(SpatData[seq(1,42,by=6),], aes(y=1:7, x=FRic)) +
    geom_point(shape=21, size=2) + looks +
    labs(y=NULL, x=NULL, title='FRic')
  cleve.plots[[4]] <- ggplot(SpatData[seq(1,42,by=6),], aes(y=1:7, x=FEve)) +
    geom_point(shape=21, size=2) + looks +
    labs(y=NULL, x=NULL, title='FEve')
  cleve.plots[[5]] <- ggplot(SpatData[seq(1,42,by=6),], aes(y=1:7, x=FDiv)) +
    geom_point(shape=21, size=2) + looks +
    labs(y=NULL, x=NULL, title='FDiv')
  
feed1 <- ggplot(SpatData, aes(x=FeedRate3)) +
  geom_histogram(bins=7, fill='white', color='black') + looks +
  labs(y='Frequency', x='Feeding Rate3')
feed2 <- ggplot(SpatData, aes(x=log(FeedRate3))) +
  geom_histogram(bins=7, fill='white', color='black') + looks +
  labs(y='Frequency', x='Feeding Rate3')
  
# boxplots[[1]] | 
boxplots[[2]] | boxplots[[3]] | boxplots[[4]]
cleve.plots[[1]] + cleve.plots[[2]] + cleve.plots[[3]] + cleve.plots[[4]] + cleve.plots[[5]]
feed1 + feed2


spat18 <- sett.trim %>% group_by(site, tile_number) %>% summarise(s2018=sum(total)) %>% ungroup() %>% select(s2018, site)
colnames(spat18)[2] <- 'Site'
spat18 <- left_join(SpatData[,1:2], spat18, by='Site')
colnames(spat18)[2] <- 's2019'
spat18 <- spat18 %>% pivot_longer(cols=c(s2018,s2019))

spat <- ggplot(spat18, aes(x=Site, y=value, fill=name)) + 
    geom_boxplot() + looks + 
  scale_x_discrete(labels=c('CB', 'L1', 'N3', 'R', 'SE', 'TB', 'V')) + 
  scale_fill_manual(name=NULL, labels=c('2018', '2019'), values=c('#DDDDDD', 'white')) +
  labs(x='Study site', y='Spat counts') + theme(legend.position=c(.87,.87), legend.title=element_blank(), legend.background=element_blank())
recs <- ggplot(RecruitData, aes(x=Site, y=Recruits)) + 
  geom_boxplot() + looks + 
  scale_x_discrete(labels=c('CB', 'L1', 'N3', 'R', 'SE', 'TB', 'V')) + 
  labs(x='Study site', y='Recruit counts')
spat + recs + plot_layout(widths=c(2,1))

## ----pairwise correlations------------------------------------------------------------------------------------------
pairs(predictors[,-1], cex=1.2, col=rgb(0,0,0,0.6))

cor.spat <- as_tibble(round(cor(predictors[-1]),2))
head(cor.spat)
cor.spat$var1 <- colnames(cor.spat)
cor.spat <- pivot_longer(as_tibble(cor.spat), -var1)
ggplot(data = cor.spat %>% filter(value != 1), aes(x=var1, y=name, fill=sqrt(value^2))) + 
  geom_tile() + looks + scale_fill_viridis_c(option='plasma', name=bquote('absolute value'~R^2)) + labs(x=NULL, y=NULL)


## ----Spat GLMM fit--------------------------------------------------------------------------------------------------------------
str(SpatData) # check that all of the data are in the proper data types before fitting
SpatData$Site <- as.factor(SpatData$Site) # Yep, site was not a factor
SpatData$FeedRate3 <- as.numeric(SpatData$FeedRate3)
set.seed(24)
# fit a GLMM with study site as a random effect (Poisson)

SpatModel <- as.list(c(0,0,0))

# SpatModel[[1]] <- glmer(data=SpatData, 
#                    Spat ~ log(FeedRate) + FEve + FDiv + FRic + (1|Site), 
#                    family=poisson)
SpatModel[[1]] <- glmer(data=SpatData, 
                   Spat ~ FeedRate3 + FEve + FDiv + FRic + (1|Site), 
                   family=poisson)
SpatModel[[2]] <- glmer(data=SpatData, 
                   Spat ~ FeedRate3 + FEve + FDiv + FRic + HerbProp + BiterProp + (1|Site), 
                   family=poisson)
SpatModel[[3]] <- glmer(data=SpatData, 
                   Spat ~ FeedRate3 + FDiv + HerbProp + (1|Site), 
                   family=poisson)
# also fit a GLMM with negative binomial because it's very likely that the data is overdispersed
# I would prefer fitting quasipoisson, but there is no function to fit that with mixed models


# look at the model fitting summaries
for (i in 2:5) {
  print(summary(SpatModel[[i]]))
}


## ----Model checks---------------------------------------------------------------------------------------------------------------
compare_performance((SpatModel[[2]]), (SpatModel[[3]]), (SpatModel[[4]]), (SpatModel[[5]]), metrics='all', rank=T) # compare the neg binom model

# look at plots of residuals to check assumptions
for (i in 2:5) {
  print(check_overdispersion(SpatModel[[i]]))
  print(check_singularity(SpatModel[[i]]))
  # diagnostic plots
  print(check_model(SpatModel[[i]]))
}



## ----Overdispersion models------------------------------------------------------------------------------------------------------
SpatModel.nb <- as.list(c(0,0))

SpatModel.nb[[1]] <- glmer.nb(data=SpatData, 
                   Spat ~ FeedRate3 + FEve + FDiv + FRic + HerbProp + BiterProp + (1|Site))
SpatModel.nb[[2]] <- glmer.nb(data=SpatData, 
                   Spat ~ FeedRate3 + HerbProp + BiterProp + (1|Site))
SpatModel.nb[[3]] <- glmer.nb(data=SpatData, 
                   Spat ~ FeedRate3 + FEve + FDiv + FRic + (1|Site))
SpatModel.nb[[4]] <- glmer.nb(data=SpatData, 
                   Spat ~ FeedRate3 + FDiv + FRic + HerbProp + (1|Site))
SpatModel.nb[[5]] <- glmer.nb(data=SpatData, 
                   Spat ~ FeedRate3 + FDiv + FRic + HerbProp + BiterProp + (1|Site))
SpatModel.nb[[6]] <- glmer.nb(data=SpatData, 
                   Spat ~ (1|Site))
SpatModel.nb[[7]] <- glmer.nb(data=SpatData, 
                   Spat ~ FeedRate3 + HerbProp + (1|Site))
SpatModel.nb[[8]] <- glmer.nb(data=SpatData, 
                   Spat ~ FeedRate3 + FEve + FDiv + FRic + BiterProp + (1|Site))
SpatModel.nb[[9]] <- glmer.nb(data=SpatData, 
                   Spat ~ FeedRate3 + FEve + FDiv + FRic + HerbProp + (1|Site))

for (i in 1:length(SpatModel.nb)) {
  print(summary(SpatModel.nb[[i]]))
}


## ----Model comparison summary---------------------------------------------------------------------------------------------------
mod.summ <- data.frame(Model=1:length(SpatModel.nb))
for (i in 1:length(SpatModel.nb)) {
mod.summ$nTerms[i] <- length(rownames(summary(SpatModel.nb[[i]])$coefficients))-1
mod.summ$AICc[i] <- round(MuMIn::AICc(SpatModel.nb[[i]]), 2)
mod.summ$dev[i] <- round(summary(SpatModel.nb[[i]])$AIC[4], 2)
mod.summ$Singularities[i] <- isSingular(SpatModel.nb[[i]])
mod.summ$Dispersion[i] <- round(summary(SpatModel.nb[[i]])$AIC[4]/df.residual(SpatModel.nb[[i]]), 2)
mod.summ$mR2[i] <- as.numeric(r2(SpatModel.nb[[i]])[2])
mod.summ$cR2[i] <- as.numeric(r2(SpatModel.nb[[i]])[1])
}
mod.summ <- mod.summ %>% arrange(AICc)
mod.summ$wAIC <- round(MuMIn::Weights(mod.summ$AICc), 3)
mod.summ$dAIC <- rep(0, length(SpatModel.nb))
for (i in 2:length(SpatModel.nb)) {
  mod.summ$dAIC[i] <- with(mod.summ, AICc[i]-AICc[i-1])
}
rownames(mod.summ) <- mod.summ$Model
print(mod.summ)


## ----Model 9 diagnostic and vis---------------------------------------------------------------------------------------------------
i=3
model_performance(SpatModel.nb[[i]])
  # diagnostic plots
  print(check_model(SpatModel.nb[[i]]))
   # overdispersion check
  # QQ plot
  par(mfrow=c(1,2))
  plot(predict(SpatModel.nb[[i]], type='link', re.form=NA), residuals(SpatModel.nb[[i]], type='deviance'), xlab='Predicted', ylab='Deviance')
    abline(h=0, lty=2)
    qqnorm(residuals(SpatModel.nb[[i]], type='deviance'))
    qqline(residuals(SpatModel.nb[[i]], type='deviance'))
    par(mfrow=c(1,1))

m=9
  l <- length(fixef(SpatModel.nb[[m]]))-1 # how many fixed effect predictors are there
  par.reg9 <- as.list(rep(0,l))
  newfit9 <- as.list(rep(0,l))
  x <- SpatData %>% ungroup() %>% mutate_if(is.numeric, mean)
  conf9 <- as.list(rep(0,l))
  bb <- as.list(rep(0,l))
  # loop inputting the confidence intervals separately in case something goes wrong with bootstrapping.
  
  for (i in 1:l) { # each predictor term
  newx <- x %>% ungroup() %>% select(!matches(names(fixef(SpatModel.nb[[m]]))[i+1])) %>% 
    bind_cols(y=SpatData %>% select(matches(names(fixef(SpatModel.nb[[m]]))[i+1])))
  # simulate new x data by just holding everything at their means
  # take out the predictor term we want to examine, and put the original back in
  newfit9[[i]] <- predict(SpatModel.nb[[m]], newx, type='response', re.form=NA)
  # predict new y values with the new x data
  
  # bootstrap confidence intervals
  bb[[i]] <- bootMer(SpatModel.nb[[m]], FUN=ci.pred, nsim=999) # bootstrap through 999 times
  # we loaded ci.pred at the beginning of the script setup
  bb_se<-apply(bb[[i]]$t,2,function(x) x[order(x)][c(.05*999, .95*999)]) # get the 5% and 95%
  conf9[[i]] <- data.frame(lwr=bb_se[1,], upr=bb_se[2,])
  }

  axis.lab <- c(
    'Functional evenness',
    'Functional richness',
    'Functional divergence',
    'Feeding rate (cm-bites/min)',
    'Herbivore abundance',
    'Benthic feeder abundance'
  )
  names(axis.lab) <- c(
    'FEve', 'FRic', 'FDiv', 'FeedRate3', 'HerbProp', 'BiterProp'
  )

for (i in 1:l) { # loop for partial regression plot panels
  par.reg9[[i]] = ggplot(bind_cols(SpatData, conf9[[i]]) %>% bind_cols(., newy=newfit9[[i]]), aes_string(y='Spat', x=names(fixef(SpatModel.nb[[m]]))[(i+1)])) +
    geom_ribbon(aes(ymax=upr, ymin=lwr), fill='grey', alpha=0.4) +
    geom_point(aes(fill=Site), size=3, alpha=0.8, shape=21, color='#777777') +
    geom_line(aes(y=newy), color='black', size=0.5) +
    labs(y='Spat counts', x=axis.lab[names(fixef(SpatModel.nb[[m]]))[(i+1)]]) + looks +
    scale_fill_viridis_d(begin=0.1, end=0.95) + theme(legend.position='none') +
    coord_cartesian(ylim=c(0,max(SpatData$Spat)+2))
}
  
  
par.reg9[[1]] + par.reg9[[2]] + par.reg9[[3]] + par.reg9[[4]] + par.reg9[[5]]
  


## ----Model 4 diagnostic plots---------------------------------------------------------------------------------------------------
i=4
model_performance(SpatModel.nb[[i]])
  # diagnostic plots
  print(check_model(SpatModel.nb[[i]]))
   # overdispersion check
  # QQ plot
  par(mfrow=c(1,2))
  plot(predict(SpatModel.nb[[i]], type='link', re.form=NA), residuals(SpatModel.nb[[i]], type='deviance'), xlab='Predicted', ylab='Deviance')
    abline(h=0, lty=2)
    qqnorm(residuals(SpatModel.nb[[i]], type='deviance'))
    qqline(residuals(SpatModel.nb[[i]], type='deviance'))
    par(mfrow=c(1,1))

# Partial regression plot
m=4
  l <- length(names(fixef(SpatModel.nb[[m]])))-1
  par.reg4 <- as.list(rep(0,l))
  newfit4 <- as.list(rep(0,l))
  x <- SpatData %>% ungroup() %>% mutate_if(is.numeric, mean)
  conf4 <- as.list(rep(0,l))
  bb <- as.list(rep(0,l))
  # loop inputting the confidence intervals separately in case something goes wrong with bootstrapping.
  
  for (i in 1:l) { # each predictor term
  newx <- x %>% ungroup() %>% select(!matches(names(fixef(SpatModel.nb[[m]]))[i+1])) %>% bind_cols(y=SpatData %>% select(matches(names(fixef(SpatModel.nb[[m]]))[i+1])))
  newfit4[[i]] <- predict(SpatModel.nb[[m]], newx, type='response', re.form=NA)
  # bootstrap confidence intervals
  bb[[i]] <- bootMer(SpatModel.nb[[m]], FUN=ci.pred, nsim=999)
  bb_se<-apply(bb[[i]]$t,2,function(x) x[order(x)][c(.05*999, .95*999)]) # dummy variable to get the 5% and 95%
  conf4[[i]] <- data.frame(lwr=bb_se[1,], upr=bb_se[2,])
  }

  for (i in 1:l) { # loop for partial regression plot panels
    par.reg4[[i]] = ggplot(bind_cols(SpatData, conf4[[i]]) %>% bind_cols(., newy=newfit4[[i]]), aes_string(y='Spat', x=names(fixef(SpatModel.nb[[m]]))[(i+1)])) +
      geom_ribbon(aes(ymax=upr, ymin=lwr), fill='grey', alpha=0.4) +
      geom_point(aes(fill=Site), size=3, alpha=0.6, shape=21, color='#777777') +
      geom_line(aes(y=newy), color='black', size=0.5) +
      labs(y='Spat counts', x=axis.lab[names(fixef(SpatModel.nb[[m]]))[(i+1)]]) + looks +
      scale_fill_viridis_d(begin=0.1, end=0.95) +
      coord_cartesian(ylim=c(0,max(SpatData$Spat)+2)) +
    if (i != l) {theme(legend.position='none')}
      else {theme(legend.position='right')}
  }
  
par.reg4[[1]] + par.reg4[[2]] + par.reg4[[3]] + par.reg4[[4]] + plot_layout(guides='collect')


## ----Recruit Models-------------------------------------------------------------------------------------------------------------

# Read in 2018 spat data
settlement <- read_xlsx('../DataCoral/settlement_2019-20_LIRS.xlsx', sheet='Sheet1', col_names=T)
head(settlement)
str(settlement)

# select just the 7 study sites and their spat data from 2018
sett.trim <- settlement %>% group_by(year, site) %>% filter(str_detect(site, 'corner_beach|lagoon_1|resort|southeast|turtle_beach|vickies|north_reef_3'), year=='2018-19')

# create a str_replace renaming object so that the site names are consistent with the predictors data
site.rename <- sort(predictors$Site)
names(site.rename) <- sort(unique(sett.trim$site))
sett.trim$site <- str_replace_all(sett.trim$site, site.rename) # replace
rm(site.rename)

# now consolidate it so that every row is one settlement tile
settlement <- sett.trim %>% group_by(site) %>% summarise(Spat2018=sum(total))
colnames(settlement)[1] <- 'Site'
predictors$Spat2018 <- settlement$Spat2018

# Import recruit data
recruits <- read_xlsx('../DataCoral/Recruits2019_cleaned.xlsx', sheet='data', col_names=T, na='NA')
str(recruits)

recruits <- recruits %>% drop_na() %>% group_by(site, quadrant) %>% summarise(Recruits=sum(Acrorec,Otherrec,Por,Fa,I)) %>% select(-quadrant)
colnames(recruits)[1] <- 'Site'
# rectify some spelling inconsistencies

predictors$Spat2018 <- scale(predictors$Spat2018, scale=T, center=F)
recruits$Site <- str_replace_all(recruits$Site, 'Vicky', 'Vickis')
recruits$Site <- str_replace_all(recruits$Site, 'Turtle$', 'TurtleBeach')
RecruitData <- na.omit(recruits) %>% full_join(y=predictors, by='Site') %>% mutate_if(is.character, as.factor) # now join it all


# collinearity check
pairs(RecruitData[,-1])
cor.recruit <- as_tibble(round(cor(RecruitData[,-1]),2))
head(cor.recruit)
cor.recruit$var1 <- colnames(cor.recruit)
cor.recruit <- pivot_longer(as_tibble(cor.recruit), -var1)
ggplot(data = cor.recruit %>% filter(value != 1), aes(x=var1, y=name, fill=sqrt(value^2))) + 
  geom_tile() + looks + scale_fill_viridis_c(option='plasma', name=bquote('absolute value'~R^2)) + labs(x=NULL, y=NULL)


## ----Poisson recruit models-----------------------------------------------------------------------------------------------------
RecModels <- as.list(rep(0,2)) # store them all in a list object

RecModels[[1]] <- glmer(data=RecruitData, family=poisson,
                        formula=Recruits ~ Spat2018 + (1|Site)) # null model in a way
RecModels[[2]] <- glmer(data=RecruitData, family=poisson,
                        formula=Recruits ~ Spat2018 + FeedRate3 + FRic + FEve + FDiv + HerbProp + (1|Site))
RecModels[[3]] <- glmer(data=RecruitData, family=poisson,
                        formula=Recruits ~ Spat2018 + FeedRate3 + FRic + FEve + FDiv + HerbProp + BiProp + (1|Site))
RecModels[[4]] <- glmer(data=RecruitData, family=poisson,
                        formula=Recruits ~ Spat2018 + FeedRate3 + HerbProp + BiProp + (1|Site))


## ----Neg binomial recruit models------------------------------------------------------------------------------------------------
RecModels.nb <- as.list(rep(0,0))

RecModels.nb[[1]] <- glmer.nb(data=RecruitData,
                        formula=Recruits ~ Spat2018 + (1|Site)) # null model in a way
RecModels.nb[[2]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + FeedRate3 + FEve + FDiv + FRic + HerbProp + BiterProp + (1|Site))
RecModels.nb[[3]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + FeedRate3 + HerbProp + BiterProp + (1|Site))
RecModels.nb[[4]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + FeedRate3 + FEve + FDiv + FRic + (1|Site))
RecModels.nb[[5]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + FeedRate3 + FDiv + FRic + HerbProp + (1|Site))
RecModels.nb[[6]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + FeedRate3 + FDiv + FRic + HerbProp + BiterProp + (1|Site))
RecModels.nb[[7]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + FeedRate3 + FDiv + HerbProp + (1|Site))
RecModels.nb[[8]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + FeedRate3 + FEve + FDiv + FRic + BiterProp + (1|Site))
RecModels.nb[[9]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + FeedRate3 + HerbProp + (1|Site))
RecModels.nb[[10]] <- glmer.nb(data=RecruitData,
                              formula=Recruits ~ (1|Site)) # null model 2

for (i in 1:length(RecModels.nb)) {
  print(summary(RecModels.nb[[i]]))
}


## ----PModel comparison summary, echo=F------------------------------------------------------------------------------------------
recmod.summ <- data.frame(Model=1:length(RecModels))
for (i in 1:length(RecModels)) {
recmod.summ$nTerms[i] <- length(rownames(summary(RecModels[[i]])$coefficients))-1
recmod.summ$AIC[i] <- round(AIC(RecModels[[i]]), 2)
recmod.summ$Dispersion[i] <- round(deviance(RecModels[[i]])/summary(RecModels[[i]])$AIC[[5]], 2)
recmod.summ$maxVIF[i] <- round(max(check_collinearity(RecModels[[i]])$VIF), 2)
}
print(as.data.frame(arrange(as_tibble(recmod.summ), AIC)))

## ----NB model comparison summary, echo=F----------------------------------------------------------------------------------------
rec.nb.summ <- data.frame(Model=1:length(RecModels.nb))
for (i in 1:length(RecModels.nb)) {
rec.nb.summ$nTerms[i] <- length(rownames(summary(RecModels.nb[[i]])$coefficients))-1
rec.nb.summ$Singularities[i] <- isSingular(RecModels.nb[[i]])
rec.nb.summ$AICc[i] <- round(MuMIn::AICc(RecModels.nb[[i]]), 2)
rec.nb.summ$mR2[i] <- round(as.numeric(r2(RecModels.nb[[i]])[2]),3)
rec.nb.summ$Dispersion[i] <- round(deviance(RecModels.nb[[i]])/summary(RecModels.nb[[i]])$AIC[5], 2)
}
rec.nb.summ <- rec.nb.summ %>% arrange(AICc)
rownames(rec.nb.summ) <- rec.nb.summ$Model
rec.nb.summ$wAIC <- round(MuMIn::Weights(rec.nb.summ$AICc), 3)
rec.nb.summ$dAIC <- rep(0, length(RecModels.nb))
for (i in 2:length(RecModels.nb)) {
  rec.nb.summ$dAIC[i] <- round(with(rec.nb.summ, AICc[i]-AICc[i-1]))
}

rec.nb.summ


## ----R.NB4 diagnostic plots, echo=F---------------------------------------------------------------------------------------------
i=9
model_performance(RecModels.nb[[i]])
check_model(RecModels.nb[[i]])
  # diagnostic plots
  #print(check_model(RecModels.nb[[i]]))
   # overdispersion check
  # QQ plot
  par(mfrow=c(1,2))
  plot(predict(RecModels.nb[[i]], type='link'), residuals(RecModels.nb[[i]], type='deviance'), xlab='Predicted', ylab='Deviance')
    abline(h=0, lty=2)
    qqnorm(residuals(RecModels.nb[[i]], type='deviance'))
    qqline(residuals(RecModels.nb[[i]], type='deviance'))
    par(mfrow=c(1,1))


## ----Partial reg plots R.NB4, echo=F--------------------------------------------------------------------------------------------
  m=9
  l <- length(fixef(RecModels.nb[[m]]))-1
  par.nb9 <- as.list(rep(0,l))
  newfit.r9 <- as.list(rep(0,l))
  conf.r9<- as.list(rep(0,l))
  x <- RecruitData %>% ungroup() %>% mutate_if(is.numeric, mean)
  
  # loop inputting the confidence intervals separately in case something goes wrong with bootstrapping.

  for (i in 1:l) { # each predictor term
  newx <- x %>% ungroup() %>% select(!matches(names(fixef(RecModels.nb[[m]]))[i+1])) %>% 
      bind_cols(y=RecruitData %>% select(matches(names(fixef(RecModels.nb[[m]]))[i+1])))
  newfit.r9[[i]] <- predict(RecModels.nb[[m]], newx, type='response', re.form=NA)
    # bootstrap confidence intervals
  bb[[i]] <- bootMer(RecModels.nb[[m]], FUN=ci.pred, nsim=999)
  bb_se<-apply(bb[[i]]$t,2,function(x) x[order(x)][c(.05*999, .95*999)]) # dummy variable to get the 5% and 95%
  conf.r9[[i]] <- data.frame(lwr=bb_se[1,], upr=bb_se[2,])
  }
  
  axis.lab <- c('2018 Spat counts',
                'Feeding rate (cm-bites/min)',
                'Herbivore abundance')
  names(axis.lab) <- names(fixef(RecModels.nb[[m]])[-1])
  
  for (i in 1:l) { # loop for partial regression plot panels
    par.nb9[[i]] = ggplot(bind_cols(RecruitData, conf.r9[[i]]) %>% bind_cols(., newy=newfit.r9[[i]]), aes_string(y='Recruits', x=names(fixef(RecModels.nb[[m]]))[(i+1)])) +
      geom_ribbon(aes(ymax=upr, ymin=lwr), fill='grey', alpha=0.4) +
      geom_point(aes(fill=Site), size=3, alpha=0.8, shape=21, color='#777777') +
      geom_line(aes(y=newy), color='black', size=0.5) +
      labs(y='Recruit counts', x=axis.lab[names(fixef(RecModels.nb[[m]]))[(i+1)]]) + looks +
      scale_fill_viridis_d(begin=0.1, end=0.95) +
      coord_cartesian(ylim=c(0,max(RecruitData$Recruits)+2)) +
      if (i != l) {theme(legend.position='none')}
    else {theme(legend.position='right')}
  }

  
(par.nb9[[1]] + par.nb9[[2]] + par.nb9[[3]])
