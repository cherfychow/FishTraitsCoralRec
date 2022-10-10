
#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Analysis: Trait data extraction
# Author: Cher Chow

#############################################################################

require(tidyverse)
require(rfishbase) # for retrieving functional trait data

# assumes these source files are loaded. File = object name
# fish_sptaxonomy.csv = fish_sp (species list)
# fish_assemblage.csv = fish_assemblage (counts by species and site)
# fish_sp <- read.csv('src/fish_sptaxonomy.csv', header = T)[-1]
# fish_assemblage <- read.csv('src/fish_assemblage.csv', header=T) %>% filter(!str_detect(Species, ' sp$'))

# Trait data --------------------------------------------------------------

# pull from Fishbase
fishbase.sp <- validate_names(sort(unique(fish_sp$Species))) # we didn't use fishbase to validate names, so we may have synonyms in the cleaned species list
sp.traits <- species(fishbase.sp) # make functional trait data frame for the species list
ecol.traits <- ecology(fishbase.sp) # table of ecological traits


## ----Trait selection------------------------------------------------------------------------------------------------------------
trait.data <- with(ecol.traits, data.frame(Species, FeedingType, FoodTroph, FoodSeTroph))
food.data <- diet(fishbase.sp) # load food data from fishbase
fooditems <- diet_items()
fooditems <- fooditems %>% arrange(DietCode, desc(DietPercent))
# species, diet code, foodI, foodII merge
diet.traits <- left_join(food.data[2:3], fooditems[c(2,4,5,7)], by="DietCode") %>% 
  select(!DietCode)

diet.traits <- diet.traits %>% filter(!DietPercent < 15) # remove any diet item that doesn't consitute at least 20%

# refine the traits manually outside of R
write_csv(trait.data, na='NA', col_names=T, path='./trait_data.csv')
rm(trait.data)
rm(food.data)
rm(morph.data)
rm(head.data)