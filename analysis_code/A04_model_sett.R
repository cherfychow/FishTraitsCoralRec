
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
require(performance)
require(readxl)
set.seed(24)

settlement <- read_xlsx('./src/LIRS_settlement_2019-20.xlsx', sheet='Sheet1', col_names=T)
# head(settlement)
# str(settlement)
# select just the 7 study sites and their spat data from 2019
sett.trim <- settlement %>% group_by(year, site) %>% filter(str_detect(site, 'corner_beach|lagoon_1|resort|southeast|turtle_beach|vickies|north_reef_3'), year=='2019-20')
# create a str_replace renaming object so that the site names are consistent with the predictors data
site.rename <- sort(predictors$Site)
names(site.rename) <- sort(unique(sett.trim$site))
sett.trim$site <- str_replace_all(sett.trim$site, site.rename) # replace

# now consolidate it so that every row is one settlement tile
settlement <- sett.trim %>% group_by(site, tile_number) %>% summarise(Spat=sum(total))
colnames(settlement)[1] <- 'Site'
SpatData <- full_join(settlement, predictors, by='Site') %>% select(-tile_number) %>% as.data.frame()

str(SpatData) # check that all of the data are in the proper data types before fitting
SpatData$Site <- as.factor(SpatData$Site) # Yep, site was not a factor
SpatData$ForRate <- as.numeric(SpatData$ForRate)

SpatModel.nb <- as.list(c(0,0))

SpatModel.nb[[1]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + TEve + TDiv + TOP + Herb + Benthic + (1|Site))
SpatModel.nb[[2]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + Herb + Benthic + (1|Site)) # failed to converge
SpatModel.nb[[3]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + TEve + TDiv + TOP + (1|Site))
SpatModel.nb[[4]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + TDiv + TOP + Herb + (1|Site))
SpatModel.nb[[5]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + TDiv + TOP + Herb + Benthic + (1|Site))
SpatModel.nb[[6]] <- glmer.nb(data=SpatData, 
                              Spat ~ (1|Site))
SpatModel.nb[[7]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + Herb + (1|Site))
SpatModel.nb[[8]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + TEve + TDiv + TOP + Benthic + (1|Site))
SpatModel.nb[[9]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + TEve + TDiv + TOP + Herb + (1|Site))

for (i in 1:length(SpatModel.nb)) {
  print(summary(SpatModel.nb[[i]]))
}

SpatModel.nb <- SpatModel.nb[-2]

## ----Model comparison summary---------------------------------------------------------------------------------------------------
mod.summ <- data.frame(Model=1:length(SpatModel.nb))
for (i in 1:length(SpatModel.nb)) {
  mod.summ$nTerms[i] <- length(rownames(summary(SpatModel.nb[[i]])$coefficients))-1
  mod.summ$AICc[i] <- round(MuMIn::AICc(SpatModel.nb[[i]]), 2)
  mod.summ$dev[i] <- round(summary(SpatModel.nb[[i]])$AIC[4], 2)
  mod.summ$Singularities[i] <- isSingular(SpatModel.nb[[i]])
  mod.summ$Dispersion[i] <- round(summary(SpatModel.nb[[i]])$AIC[4]/df.residual(SpatModel.nb[[i]]), 2)
  mod.summ$mR2[i] <- as.numeric(r2(SpatModel.nb[[i]])[2]) %>% round(.,3)
}
mod.summ <- mod.summ %>% arrange(AICc)
mod.summ$wAIC <- round(MuMIn::Weights(mod.summ$AICc), 3)
mod.summ$dAIC <- mod.summ$AICc-mod.summ$AICc[1] %>% round(., 3)
rownames(mod.summ) <- mod.summ$Model
print(mod.summ)

i=4
model_performance(SpatModel.nb[[i]])
# diagnostic plots
print(check_model(SpatModel.nb[[i]]))
# overdispersion check
# QQ plot
par(mfrow=c(1,2))
plot(predict(SpatModel.nb[[i]], type='link', re.form=NA), residuals(SpatModel.nb[[i]], type='deviance'), xlab='Predicted', ylab='Deviance')
abline(h=0, lty=2)
qqnorm(residuals(SpatModel.nb[[i]], type='deviance'))
qqline(residuals(SpatModel.nb[[i]], type='deviance'))
par(mfrow=c(1,1))

