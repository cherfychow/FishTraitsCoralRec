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
ubruv.trim <- filter(ubruv.data, str_detect(Species, 'idae|unknown', negate = TRUE))
# remove uncertain observations at Family level or higher
# 468 observations left, removed 7 rows
# We went from ubruv.data > ubruv.trim to remove uncertain records
rm(ubruv.data) # we can remove this object now


## ----taxon spell check, results='hide'------------------------------------------------------------------------------------------

Sp <- sort(unique(str_replace(ubruv.trim$Species, ' sp$', ''))) # e.g. change "Scarus sp" to "Scarus" for WoRMS validation
Sp <- unique(str_replace(Sp, ' $', '')) # remove accidental duplicates ending with spaces
require(worms) # package for validating species/taxa names

## Warning: WoRMS is computationally heavy! Might take a few minutes
Sp.list <- wormsbymatchnames(Sp, marine_only=T, chunksize = 49) # retrieve WoRMS checks
Sp.list$ID <- seq(1, length(Sp.list$valid_name), by=1) # create an ID column

# create a list of invalid species names to check
invalid <- Sp.list %>% filter(match_type != 'exact' | status != 'accepted') %>% select(ID, valid_name, family, genus)
invalid$before <- Sp[invalid$ID] # retrieve the invalid species names
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
# and now from ubruv.trim > ubruv.validated after name checks
rm(ubruv.trim)
rm(invalid)


## ----species summary------------------------------------------------------------------------------------------------------------
Sp.list <- Sp.list %>% select(valid_name, family) %>% mutate(genus=word(valid_name, 1))
# let's trim the worms table to just taxonomic info, and extract genus from the valid name

Sp.count <- as_tibble(ubruv.validated) %>% 
  group_by(SpeciesV) %>% summarise(n=sum(Count)) %>% 
  filter(str_detect(SpeciesV, ' sp$', negate=T))# counts per species

Sp.countL <- as_tibble(ubruv.validated) %>% 
  group_by(SpeciesV, Length) %>% summarise(n=sum(Count)) # counts per species and length class

ubruv.validated %>% group_by(Site) %>% summarise(n=sum(Count)) # total number of fish per site

UBRUVspecies <- ubruv.validated %>% group_by(Site, SpeciesV) %>% summarise(n=sum(Count)) # counts per species grouped by site
UBRUVspeciesL <- ubruv.validated %>% group_by(Site, SpeciesV, Length) %>% summarise(n=sum(Count))
# make a separate tibble for grouping by length classes


colnames(Sp.list)[1] <- "SpeciesV" #rename the Species column in the worms dataframe
UBRUVspecies$genus <- word(UBRUVspecies$SpeciesV, 1) # extract genus names to make join easier
UBRUVspecies <- left_join(UBRUVspecies, distinct(Sp.list[,2:3]), by='genus') # add a column of family names to UBRUV datasheet

Sp.count$genus <- word(Sp.count$SpeciesV, 1) # add genus column to the species count summary tibble
Sp.count <- left_join(Sp.count, distinct(Sp.list[,2:3]), by='genus') # and add family taxonomic info too


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