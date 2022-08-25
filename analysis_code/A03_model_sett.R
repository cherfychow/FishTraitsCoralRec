
#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Coral settlement model selection
# Author: Cher Chow

#############################################################################

# ForRate = trait-weighted foraging rates
# TEve = Trait evenness
# TDiv = Trait divergence
# TOP = Trait richness
# Herb = Relative herbivore abundance
# Benthic = Relative benthic biter abundance

require(tidyverse)
require(lme4)
require(MuMIn)
require(performance)
require(readxl)
set.seed(24)

# predictor collinearity check
cor.sett <- as_tibble(round(cor(predictors[-1]),2))
head(cor.sett)
cor.sett$var1 <- colnames(cor.sett)
cor.sett <- pivot_longer(as_tibble(cor.sett), -var1)
ggplot(data = cor.sett %>% filter(value != 1), aes(x=var1, y=name, fill=sqrt(value^2))) + 
  geom_tile() + looks + scale_fill_viridis_c(name=bquote('absolute value'~R^2)) + labs(x=NULL, y=NULL)

settlement <- read_xlsx('./src/coral_settlement.xlsx', sheet='Sheet1', col_names=T)
# head(settlement)
# str(settlement)
# select just the 7 study sites and their spat data from 2019
settlement <- settlement %>% group_by(year, site) %>% filter(str_detect(site, 'corner_beach|lagoon_1|resort|southeast|turtle_beach|vickies|north_reef_3'), year=='2019-20')
# create a str_replace renaming object so that the site names are consistent with the predictors data
site.rename <- sort(predictors$Site)
names(site.rename) <- sort(unique(settlement$site))
settlement$site <- str_replace_all(settlement$site, site.rename) # replace

# now consolidate it so that every row is one settlement tile
settlement <- settlement %>% group_by(site, tile_number) %>% summarise(Spat=sum(total))
colnames(settlement)[1] <- 'Site'
SpatData <- full_join(settlement, predictors, by='Site') %>% select(-tile_number) %>% as.data.frame()

str(SpatData) # check that all of the data are in the proper data types before fitting
SpatData$Site <- as.factor(SpatData$Site) # Yep, site was not a factor
SpatData$ForRate <- as.numeric(SpatData$ForRate)

sett_global <- glmer.nb(data=SpatData, 
                        Spat ~ ForRate + TEve + TDiv + TOP + Herb + Benthic + (1|Site), na.action = 'na.fail')
# use dredge to run model selection with all possible predictor combinations (with site fixed)
sett.select <- dredge(sett_global, beta = 'none', evaluate = T, rank = 'AICc', fixed = c('ForRate, Site'))
# sett.select2 <- dredge(sett_global, beta = 'none', evaluate = T, rank = 'BIC', fixed = 'Site')
sett.select
sett <- get.models(sett.select, 1)[[1]] # print results for top model
summary(sett)
model_performance(sett)
# diagnostic plots
print(check_model(get.models(sett.select, 1)[[1]]))
# overdispersion check
# QQ plot
par(mfrow=c(1,2))
plot(predict(Model.Sett[[i]], type='link', re.form=NA), residuals(Model.Sett[[i]], type='deviance'), xlab='Predicted', ylab='Deviance')
abline(h=0, lty=2)
qqnorm(residuals(Model.Sett[[i]], type='deviance'))
qqline(residuals(Model.Sett[[i]], type='deviance'))
par(mfrow=c(1,1))

