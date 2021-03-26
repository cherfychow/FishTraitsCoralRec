
#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Generate trait space for trait diversity analysis
# Author: Cher Chow
# Updated: 26 Mar 2021

#############################################################################

require(tidyverse)
require(FD) # calculate dissimilarity matrices + functional diversity indices
require(geometry) # calculates convex hull/trait space for data vis

set.seed(24) # repeatability :)

# Pairwise trait check, just so they're not directly colinear
pairs(traits[,-1], labels=colnames(traits)[-1], col=rgb(0.2,0.2,0.2, 0.5))


## Trait space ------------------------------------------------------------------------------------------------------------

# create the abundance matrix for functional evenness and divergence weighting
# species are columns, each site + total = rows
abund <- Sp.count %>% filter(str_detect(SpeciesV, ' sp$', negate=T)) %>% select(SpeciesV, n) %>% spread(SpeciesV, n) # row 1 = all sites
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