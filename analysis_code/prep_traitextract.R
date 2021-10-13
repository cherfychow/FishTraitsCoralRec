
#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Analysis: Trait data extraction
# Author: Cher Chow

#############################################################################

require(tidyverse)
require(rfishbase) # for retrieving functional trait data

# Trait data --------------------------------------------------------------

# pull from Fishbase
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