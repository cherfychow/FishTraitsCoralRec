#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Trait weighted foraging rates analysis
# Author: Cher Chow

#############################################################################

require(dplyr)
require(stringr)
set.seed(24)

# load data
influence <- read.csv('./src/fish_traitsinfluence.csv', header=T) # the modified trait influence table
fish_bites <- read.csv('./src/fish_bites.csv', header=T)
str(influence)

sort(unique(fish_bites$Length))
# make a column for length classifications that's for treating length as discrete scales
length.convert <- c('2.5', '15', '25', '7.5') # parallels what is shown with sort(unique(fish_bites$Length))
names(length.convert) <- paste('^', sort(unique(fish_bites$Length)), '$', sep='')
fish_bites$Length <- str_replace_all(fish_bites$Length, length.convert)
# convert length class to the median length to estimate foraging rates

BiteSp <- unique(fish_bites$Species)
influence <- influence %>% filter(Species %in% BiteSp) # only keep the trait influence rows for the biting species


# Trait-weighted factors --------------------------------------------------

# we use the influence version of the trait table to create a trait-weighted scaling coefficient for the bite rates
# run PCA on the influence trait table
inf.pca <- princomp(influence[,-1])
loadings(inf.pca) # loadings show contribution of the traits to the principal components
summary(inf.pca)
screeplot(inf.pca, type='lines') # we want to use just PC1, but check that it captures the majority of the variation
# should see a sharp drop

# PC1 valid
# extract the species' position on PC1
sp.inf <-  data.frame(Species=influence$Species, SpInf=inf.pca$scores[,1])
# do some arbitrary scaling to force it within 0-1
# 0.1 is an arbitrary factor to keep foraging rates on a comparable effect size with other trait indices
sp.inf$SpInf <- (sp.inf$SpInf + (-range(sp.inf$SpInf)[1]) + 0.1)/2


# Nested aggregation ------------------------------------------------------

# Foraging rate_j = sum_ij(Species trait coefficient_i * sum_l(median length_il * bite rate_il))
# site = j
# species = i
# length class = l

# this step calculates sum_l(median length_il * bite rate_il)
ForRate <- fish_bites %>% ungroup() %>% group_by(Site, Species, Length) %>% 
  select(Site, Species, Length, BiteRate) %>%  
  mutate(LBi=as.numeric(Length)*BiteRate) %>%  # median length_il * bite rate_il
  ungroup()
ForRate <- ForRate %>% group_by(Site, Species) %>% 
  summarise(sumLBi=sum(LBi)) # aggregate records by Site and Species, sum_l

# Now use the calculated species influence weighting factors to adjust bite rates to foraging rates
# sum_ij
predictors$ForRate <- left_join(ForRate, sp.inf, by='Species') %>% 
  mutate(SpFeed=0.01*SpInf*sumLBi) %>% group_by(Site) %>% 
  summarise(ForRate=sum(SpFeed)) %>% ungroup() %>% pull(ForRate) %>% as.vector
predictors$ForRate <- scale(predictors$ForRate, scale=T, center=F)

rm(BiteSp, length.convert, sp.inf, influence, inf.pca)
# remove these objects since they're not dependencies for downstream sourcing