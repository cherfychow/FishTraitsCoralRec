
#############################################################################

# FishTraitsCoralRec
# Analysis: Trait data extraction

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
fooditems2 <- fooditems(species_list = fishbase.sp) # not all species have diet data, but do have food items data
fooditems2 <- fooditems2 %>% select(Species, Locality, FoodII, FoodIII, starts_with('Common'), Foodname)
fooditems <- fooditems %>% arrange(DietCode, desc(DietPercent))
# species, diet code, foodI, foodII merge
diet.traits <- left_join(food.data[2:3], fooditems[c(2,4,5,7)], by="DietCode") %>% 
  select(!DietCode)

diet.traits <- diet.traits %>% filter(!DietPercent < 15) %>%  # remove any diet item that doesn't consitute at least 20%
  arrange(Species, desc(DietPercent))

# load Parraviccini et al. (2020) trophic groupings
mod_predators <- read.csv("https://github.com/valerianoparravicini/Trophic_Fish_2020/raw/master/output/results/mod_predators.csv")
troph <- c('1' = 'sess_inv',
           '2' = 'herb_mic_det',
           '3' = 'cor',
           '4' = 'pisc',
           '5' = 'micro_inv',
           '6' = 'macro_inv',
           '7' = 'crust',
           '8' = 'plank') # grouping key from paper
trait.data <- left_join(trait.data, mod_predators, by=c("Species" = "species"))
trait.data$cluster <- str_replace_all(trait.data$cluster, troph)

n_distinct(diet.traits$Species) + nrow(fooditems2 %>% filter(is.na(CommonessII) == F, Species %in% diet.traits$Species == F) %>% distinct(Species)) # no of species with food item dominance info
n_distinct(fish_sp$Species) - 72 # 32 species without dominance data (informed by trophic grouping)

# refine the traits manually outside of R
write_csv(trait.data, na='NA', col_names=T, path='./trait_data.csv')
rm(list = ls()) # clear environment