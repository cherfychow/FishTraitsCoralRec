
#############################################################################

# FishTraitsCoralRec

# Construct trait space + calculate trait diversity metrics

#############################################################################

require(tidyverse)
require(FD) # calculate dissimilarity matrices + functional diversity indices

# working directory is repo home
set.seed(24) # repeatability :)

# load data
traits <- read.csv('src/fish_traits.csv', header=T)
rownames(traits) <- traits$Species
# we can't run trait analysis on uncertain species IDs
fish_assemblage <- read.csv('src/fish_assemblage.csv', header=T) %>% filter(!str_detect(Species, ' sp$'))
fish_sp <- read.csv('src/fish_sptaxonomy.csv', header=T)[-1]

# Pairwise trait check, just so they're not directly colinear
# pairs(traits[,-1], labels=colnames(traits)[-1], col=rgb(0.2,0.2,0.2, 0.5))

## Trait space set up ------------------------------------------------------------------------------------------------------------

# create the abundance matrix for functional evenness and divergence weighting
# species are columns, each site is a row
# long format to wide
abund <- fish_sp %>% filter(!str_detect(Species, ' sp$')) %>% 
  select(Species, n) %>% 
  pivot_wider(names_from = Species, values_from = n) # row 1 = all sites
abund$Site <- 'All'
abund <- rbind(fish_assemblage %>% filter(!str_detect(Species, ' sp$')) %>% select(Site, Species, n) %>% spread(Species, n), abund)
# because not all species occur in each site, there are lots of NAs to be filled with 0s
NAlist <- as.list(rep(0,dim(abund)[2]-1)) # list has to have column names for every item = 0
names(NAlist) <- colnames(abund)[2:ncol(abund)]
abund <- abund %>% replace_na(NAlist) # annnd replace NAs
rownames(abund) <- abund$Site # FD package reads row names. move over.
abund$Site <- NULL
rm(NAlist)

# Trait space + diversity metrics -------------------------------------------------

# calculate FD indices for each site with respect to the global functional trait space
# because this is with respect to all study sites, sites are comparable
# not calculating TRichness because we're using an alternative calculation
trait.dis <- cailliez(gowdis(traits[,-1], ord='podani'))
is.euclid(trait.dis) # check whether distances are Euclidean before running PCoA

FD_global <- dbFD(trait.dis, abund, w.abun=T, calc.FRic=T, calc.FDiv=T, m=4, calc.CWM=F, calc.FGR=F, print.pco=T)
point <- as.data.frame(FD_global$x.axes)[,1:4] # create a data frame of each species position in the 4D trait space
point$Species <- rownames(traits)

## If you were to run PCoA checks, this is the time ##

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
  s.point[[i]] <- fish_assemblage %>% filter(Site == unique(Site)[i]) %>% ungroup() %>% select(Species)
  s.point[[i]] <- left_join(s.point[[i]], point, by='Species')# use the point df to match the positions from the PCoA
  predictors$TOP[i] <- TOP.index(s.point[[i]][-1])[2] # calculate the TOP index
}

# now we have to standardise this by the global trait space TOP
predictors$TOP <- predictors$TOP/TOP.index(point[1:4])[2]

# the TOP index function generates a vert.txt file but we don't need that anymore
if (file.exists('vert.txt')) {
  #Delete file if it exists
  invisible(file.remove('vert.txt'))
}

# Calculate relative abundances herbivore biters --------------------------


site.total <- data.frame(site.total=rowSums(abund), Site=rownames(abund))
# reference vector of herbivore species
herb_spp <- traits %>% filter(TrophParr == 'herb_mic') %>% pull(Species)
site.herb <- fish_assemblage %>% filter(Species %in% herb_spp) # filter assemblage data by herbivores
site.herb <- site.herb %>% group_by(Site) %>% summarise(Herb = sum(n)) %>% # add them all up
  left_join(., site.total, by="Site") %>% 
  mutate(Herb = Herb/site.total) %>% select(!site.total) # divide by total for rel abundance
# herbivore abundance

# biter abundance
# just pull out the sessile invertivore species
bite_spp <- traits %>% filter(str_detect(TrophParr, 'sess_inv')) %>% pull(Species)

site.biter <- fish_assemblage %>% filter(Species %in% bite_spp) # only 4 sites, but we'll handle the zeroes later
site.biter <- site.biter %>% 
  group_by(Site) %>% summarise(Benthic = sum(n)) %>% # add them all up
  left_join(., site.total, by="Site") %>% 
  mutate(Benthic = Benthic/site.total) %>% select(!site.total) # divide by total for rel abundance
site.biter$site.total <- NULL

# add it to our predictors dataframe
predictors$Herb <- site.herb$Herb
predictors <- left_join(predictors, site.biter, by="Site")
predictors$Benthic[which(is.na(predictors$Benthic))] <- 0 # instead of NAs, I want zeroes

# Species relative contribution -------------------------------------------

# get an idea of species contributions to trait diversity metrics by doing an iterative take one out analysis
spcont <- data.frame(Site = '', Species = '', TEve = '', TDiv = '')
spcont <- list(spcont, spcont, spcont, spcont, spcont, spcont, spcont)
# evenness and divergence first
for (j in 1:ncol(abund)) {
  # run a trait space without species j
  dummydis <-  as.matrix(trait.dis)[-j,-j] %>% as.dist # because it's a symmetric matrix, we remove row and column j
  # doesn't rerun PCoA this way
  traitspace_sp <- dbFD(dummydis, abund[-j], w.abun=T, calc.FRic=T, calc.FDiv=T, m=4, calc.CWM=F, calc.FGR=F, print.pco=F)
  for (i in 1:7) {
    spcont[[i]] <- rbind(spcont[[i]], data.frame(Site = predictors$Site[i],
                                                 Species = traits$Species[j],
                                                 TEve = predictors$TEve[i] - traitspace_sp$FEve[i], 
                                                 TDiv = predictors$TDiv[i] - traitspace_sp$FDiv[i]))
  }
}

for (i in 1:7) { # get rid of the first blank row 
  spcont[[i]] <- spcont[[i]][-1,] %>% 
    filter(Species %in% fish_assemblage$Species[which(fish_assemblage$Site == predictors$Site[i])])
}

# now calculate contributions to TOP measures
for (i in 1:7) {
  dTOP <- ''
  for (j in 1:nrow(s.point[[i]])) {
    dTOP <- c(dTOP, TOP.index(s.point[[i]][-j,-1])[2]) # calculate delta TOP for each species we take out
  }
  spcont[[i]]$TOP <- predictors$TOP[i] - (as.numeric(dTOP[-1]) / TOP.index(point[1:4])[2]) # store the values into the spcont data frames
}

sp_contr <- do.call(rbind, spcont)
rownames(sp_contr) <- 1:nrow(sp_contr)
sp_contr$TEve <- as.numeric(sp_contr$TEve)
sp_contr$TDiv <- as.numeric(sp_contr$TDiv)

# the TOP index function generates a vert.txt file but we don't need that anymore
if (file.exists('vert.txt')) {
  #Delete file if it exists
  invisible(file.remove('vert.txt'))
}

# round the indices contributions
sp_contr$TEve <- round(sp_contr$TEve, 5)
sp_contr$TDiv <- round(sp_contr$TDiv, 5)
sp_contr$TOP <- round(sp_contr$TOP, 5)
sp_contr <- left_join(sp_contr, predictors[1:4], by="Site")
sp_contr_rel <- sp_contr %>% mutate(dTOP = TOP.x / TOP.y, 
                                    dTEve = TEve.x / TEve.y, 
                                    dTDiv = TDiv.x / TDiv.y) %>% 
  select(Species, Site, dTOP, dTEve, dTDiv) %>% arrange(Site, desc(dTDiv))# make it relative 

# export
# write.csv(sp_contr, 'outputs/sp_traitdiv_cont.csv', row.names = F)

# clean up objects that aren't dependencies for downstream sourcing
rm(herb_spp, bite_spp, site.biter, site.herb,
   abund, spcont, dummydis, traitspace_sp, dTOP) # remove these dummy temp objects
