
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

Model.Sett <- as.list(c(0,0))

Model.Sett[[1]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + TEve + TDiv + TOP + Herb + Benthic + (1|Site))
Model.Sett[[2]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + Herb + Benthic + (1|Site)) # failed to converge
Model.Sett[[3]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + TEve + TDiv + TOP + (1|Site))
Model.Sett[[4]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + TDiv + TOP + Herb + (1|Site))
Model.Sett[[5]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + TDiv + TOP + Herb + Benthic + (1|Site))
Model.Sett[[6]] <- glmer.nb(data=SpatData, 
                              Spat ~ (1|Site))
Model.Sett[[7]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + Herb + (1|Site))
Model.Sett[[8]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + TEve + TDiv + TOP + Benthic + (1|Site))
Model.Sett[[9]] <- glmer.nb(data=SpatData, 
                              Spat ~ ForRate + TEve + TDiv + TOP + Herb + (1|Site))

for (i in 1:length(Model.Sett)) {
  print(summary(Model.Sett[[i]]))
}

Model.Sett <- Model.Sett[-2]

## Model comparison summary
summary.sett <- data.frame(Model=1:length(Model.Sett))
for (i in 1:length(Model.Sett)) {
  summary.sett$nTerms[i] <- length(rownames(summary(Model.Sett[[i]])$coefficients))-1
  summary.sett$AICc[i] <- round(MuMIn::AICc(Model.Sett[[i]]), 2)
  summary.sett$BIC[i] <- round(BIC(Model.Sett[[i]]), 2)
  summary.sett$dev[i] <- round(summary(Model.Sett[[i]])$AIC[4], 2)
  summary.sett$Dispersion[i] <- round(summary(Model.Sett[[i]])$AIC[4]/df.residual(Model.Sett[[i]]), 2)
}
summary.sett <- summary.sett %>% arrange(AICc, BIC)
summary.sett$wAIC <- round(MuMIn::Weights(summary.sett$AICc), 3)
summary.sett$dAIC <- summary.sett$AICc-summary.sett$AICc[1] %>% round(., 3)
summary.sett$dBIC <- summary.sett$BIC-summary.sett$BIC[1] %>% round(., 3)
rownames(summary.sett) <- summary.sett$Model
print(summary.sett)

i=4
model_performance(Model.Sett[[i]])
# diagnostic plots
print(check_model(Model.Sett[[i]]))
# overdispersion check
# QQ plot
par(mfrow=c(1,2))
plot(predict(Model.Sett[[i]], type='link', re.form=NA), residuals(Model.Sett[[i]], type='deviance'), xlab='Predicted', ylab='Deviance')
abline(h=0, lty=2)
qqnorm(residuals(Model.Sett[[i]], type='deviance'))
qqline(residuals(Model.Sett[[i]], type='deviance'))
par(mfrow=c(1,1))

