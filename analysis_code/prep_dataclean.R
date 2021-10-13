#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Data cleanup
# Author: Cher Chow

#############################################################################


# Setup -------------------------------------------------------------------

require(tidyverse) # data handling
require(readxl) # read xls files
require(worms) # package for validating species/taxa names
require(lubridate) # use lubridate to deal with time data

set.seed(24) # repeatability :)

# Read in data ------------------------------------------------------------

# read in bite data from folder
files = list.files(path='./src/RUV_bite', pattern='*.csv', full.names = T)
bite.sheets = lapply(files, read.csv, header=T) # read in all the xls files in the folder
bite.data <- do.call(rbind, bite.sheets)
head(bite.data, n=10)
rm(bite.sheets)
rm(files)

# read in assemblage data from folder
files = list.files(path='./src/RUV_spL/', pattern='*.xls', full.names = T)
sp.sheets = lapply(files, read_xlsx, sheet=1, col_names=T, col_types=c('date',rep('text', 6),'numeric',rep('text', 4))) # read in all the xls files in the folder
ubruv.data <- do.call(rbind, sp.sheets)
head(ubruv.data)
rm(sp.sheets)
rm(files)


# Check and trim species data-------------------------------------------------------------------------------------------------------------
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


## Taxon checks------------------------------------------------------------------------------------------

Sp <- sort(unique(str_replace(ubruv.trim$Species, ' sp$', ''))) # e.g. change "Scarus sp" to "Scarus" for WoRMS validation
Sp <- unique(str_replace(Sp, ' $', '')) # remove accidental duplicates ending with spaces

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
rm(ubruv.trim, invalid, validate, corrections, Sp)


## ----species summary------------------------------------------------------------------------------------------------------------
Sp.list <- Sp.list %>% select(valid_name, family) %>% mutate(genus=word(valid_name, 1))
# let's trim the worms table to just taxonomic info, and extract genus from the valid name

Sp.count <- as_tibble(ubruv.validated) %>% 
  group_by(SpeciesV) %>% summarise(n=sum(Count)) %>% 
  filter(str_detect(SpeciesV, ' sp$', negate=T))# counts per species
names(Sp.count)[1] = "Species"

site.n <- ubruv.validated %>% group_by(Site) %>% summarise(n=sum(Count)) # total number of fish per site

UBRUVspecies <- ubruv.validated %>% group_by(Site, SpeciesV) %>% summarise(n=sum(Count)) # counts per species grouped by site

colnames(Sp.list)[1] <- "Species" #rename the Species column in the worms dataframe
colnames(UBRUVspecies)[2] <- "Species" #rename the Species column in the worms dataframe

Sp.count$genus <- word(Sp.count$Species, 1) # add genus column to the species count summary tibble
Sp.count <- left_join(Sp.count, distinct(Sp.list[,2:3]), by='genus') # and add family taxonomic info too

# uncomment to save
# write.csv(Sp.count[c(4,3,1,2)], 'src/fish_sptaxonomy.csv')
# write.csv(UBRUVspecies, 'src/fish_assemblage.csv')

## Validate names for bite data------------------------------------------------------------------------

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


## ----Duration data clean--------------------------------------------------------------------------------------------------------
bite.validated$BiStart <- with(bite.validated, ms(BiStart))
bite.validated$BiEnd <- with(bite.validated, ms(BiEnd))
bite.validated$Duration <- with(bite.validated, as.numeric(as.duration(BiEnd)-as.duration(BiStart))) # Duration in seconds!
head(bite.validated)
# calculate duration of occurrences/bites and add as new column

# clean up the length categories
bite.validated$Length <- str_replace_all(bite.validated$Length, ' 10_20$', '10_20')
bite.validated$Length <- str_replace_all(bite.validated$Length, '0_5$', '_5')

# # a separate tibble for bites grouped by LENGTH, species, and site.
UBRUVbiteL <- bite.validated %>% group_by(Site, Species, Length) %>% summarise(BiteTotal=sum(BiTotal), DurationTotal=sum(Duration), Bitemaxn=max(nInd)) %>% filter(BiteTotal > 0) # only keep observations of species that bite
UBRUVbiteL$BiteRate <- with(UBRUVbiteL, (BiteTotal/DurationTotal)*60) # bites per minute
head(UBRUVbiteL)

# uncomment to save
# write.csv(UBRUVbiteL, 'src/fish_bites.csv', row.names=F, col.names=T)

