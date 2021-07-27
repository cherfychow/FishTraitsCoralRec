
#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Generate trait space for trait diversity analysis
# Author: Cher Chow

#############################################################################

require(tidyverse)
require(FD) # calculate dissimilarity matrices + functional diversity indices

# plot aesthetics
looks <- theme_bw(base_size=13) + theme(panel.grid=element_blank(), axis.ticks=element_line(size=0.3))

# working directory is repo home
set.seed(24) # repeatability :)

# Pairwise trait check, just so they're not directly colinear
# pairs(traits[,-1], labels=colnames(traits)[-1], col=rgb(0.2,0.2,0.2, 0.5))


## Trait space ------------------------------------------------------------------------------------------------------------

# create the abundance matrix for functional evenness and divergence weighting
# species are columns, each site + total = rows
abund <- Sp.count %>% filter(str_detect(SpeciesV, ' sp$', negate=T)) %>% select(SpeciesV, n) %>% spread(SpeciesV, n) # row 1 = all sites
abund$Site <- 'All'
abund <- rbind(UBRUVspecies %>% select(Site, Species, n) %>% spread(Species, n), abund)
# because not all species occur in each site, there are lots of NAs to be filled with 0s
NAlist <- as.list(rep('0',dim(abund)[2]-1)) # list has to have column names for every item = 0
names(NAlist) <- colnames(abund)[2:ncol(abund)]
abund <- abund %>% replace_na(NAlist) # annnd replace NAs
rownames(traits) <- traits$Species # dbFD function requires species in the row names

# abundance data frame is in characters????
abund <- abund %>% mutate(across(where(is.character), as.numeric))
abund <- as.data.frame(abund)
rownames(abund) <- abund$Site
abund$Site <- NULL
rm(NAlist)

## ----Trait space sites, results='hide'------------------------------------------------------------------------------------------
# running the full lizard island traits together is very computationally taxing
# so I'll split all the sites to run the FD package separately
abun <- as.list(rep('0',7))
for (i in 1:7) {
  abun[[i]] <- as.data.frame(UBRUV %>% filter(Site == unique(UBRUV$Site)[i]) %>% ungroup() %>% select(Species, n) %>% spread(Species, n))
}
names(abun) <- unique(UBRUV$Site)
names(Sp.count)[1] <- 'Species'
names(Sp.countL)[1] <- 'Species'
abun[[8]] <- Sp.count %>% filter(str_detect(Species, ' sp$', negate=T)) %>% select(Species, n) %>% spread(Species, n)
names(abun)[8] <- 'all' # add one item for all sites

trait.dis <- cailliez(gowdis(traits[,-1], ord='podani'))
is.euclid(trait.dis) # check whether distances are Euclidean before running PCoA

# run PCoA
# pcoa <-  pcoa(trait.dis, rn=colnames(trait.dis)) #PCoA with Cailliez corrections to negative eigenvalues
# this is primarily so we have a PCoA object for checking. FD package wraps pcoa internally in the function

## ----FD calculation-------------------------------------------------------------------------------------------------------------

# calculate FD indices for each site with respect to the global functional trait space
# because this is with respect to all study sites, sites are comparable
# except for TRichness because we're using an alternative calculation
FD_global <- dbFD(trait.dis, abund, w.abun=T, calc.FRic=T, calc.FDiv=T, m=4, calc.CWM=F, calc.FGR=F, print.pco=T)
point <- as.data.frame(FD_global$x.axes)[,1:4] # create a data frame of each species position in the 4D trait space
point$Species <- rownames(traits)

# extract the FD indices from the FD_global object, create a dataframe of them
predictors <- data.frame(Site=names(FD_global$FEve[1:7]),TEve=FD_global$FEve[1:7], TDiv=FD_global$FDiv[1:7])
# only 7 for the 7 sites. Ignore the global FD measures
rownames(predictors) <- names(FD_global$FEve[1:7])


# construct data frame of trait space points for each site
# and then that data frame serves as the input to the TOP metric function by Fontana et al. 2015
# Source: doi.org/10.1111/1365-2435.12551
source('analysis_code/TOP_Fontanaetal2015.R')
s.point <- as.list(rep(0,7))
predictors$TOP <- rep(0, 7)
for (i in 1:7) {
  # first, make a dummy points dataframe only with the species in that site
  s.point[[i]] <- UBRUVspecies %>% filter(Site == unique(UBRUVspecies$Site)[i]) %>% ungroup() %>% select(Species)
  s.point[[i]] <- left_join(s.point[[i]], point, by='Species')# use the point df to match the positions from the PCoA
  predictors$TOP[i] <- TOP.index(s.point[[i]][-1])[2] # calculate the TOP index
}

# now we have to standardise this by the global trait space TOP
predictors$TOP <- predictors$TOP/TOP.index(point[1:4])[2]

# Calculate relative abundances herbivore biters --------------------------

site.total <- data.frame(site.total=rowSums(abund), Site=rownames(abund))
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

# herbivore abundance
site.herb <- data.frame(Site=predictors$Site, Herb=rep(0,7))
for (i in 1:7) {
  herb <- as.data.frame(traits.sites[[i]])
  herb$Species <- rownames(herb)
  herb <-  herb %>% filter(TrophicGroup=='herbivore') %>% ungroup() %>% select(Species)
  herb$Site <- names(traits.sites)[i]
  herb <- left_join(herb, UBRUVspecies, by=c('Site', 'Species')) %>% ungroup() %>% select(Site, Species, n)
  site.herb$Herb[i] <- sum(herb$n)/site.total$site.total[i]
}

# biter abundance
site.biter <- data.frame(Site=predictors$Site, Herb=rep(0,7))
for (i in 1:7) {
  biter <- as.data.frame(traits.sites[[i]])
  biter$Species <- rownames(biter)
  biter <-  biter %>% filter(ColumnFeed == 'benthic') %>% 
    filter(str_detect(TrophicGroup, 'detritovore|invertivore|omnivore')) %>% ungroup() %>% select(Species)
  biter$Site <- names(traits.sites)[i]
  biter <- left_join(biter, UBRUVspecies, by=c('Site', 'Species')) %>% ungroup() %>% select(Site, Species, n)
  site.biter$Benthic[i] <- sum(biter$n)/site.total$site.total[i]
}
# add it to our predictors dataframe
predictors$Herb <- site.herb$Herb
predictors$Benthic <- site.biter$Benthic

