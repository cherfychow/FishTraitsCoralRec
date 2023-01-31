
#############################################################################

# FishTraitsxCoralRec
# Coral settlement model construction and selection
# Author: Cher Chow

#############################################################################

# ForRate = trait-weighted foraging rates
# TEve = Trait evenness
# TDiv = Trait divergence
# TOP = Trait richness
# Herb = Relative herbivore abundance
# Benthic = Relative benthic biter abundance

require(dplyr)
require(tidyr)
require(ggplot2)
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
looks <- theme_bw(base_size=13) + theme(panel.grid=element_blank(), axis.ticks=element_line(size=0.3))
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
# use dredge to run model selection with all possible predictor combinations 
# fixed variables i.e. constant predictors in all candidate combinations = ForRate + Site
sett.select <- dredge(sett_global, beta = 'none', evaluate = T, rank = 'AICc', fixed = c('ForRate', 'Site'))
sett.select # selection table

## Custom selection table
## Model comparison summary
top10 <- get.models(sett.select, 1:10)
top10[[11]] <- glmer.nb(data=SpatData, 
                        Spat ~ (1|Site), na.action = 'na.fail')
summary.sett <- data.frame(ModelRank=1:11) # just top 10 candidates
for (i in 1:11) {
  summary.sett$nTerms[i] <- getAllTerms(top10[[i]], intercept = F) %>% length
  summary.sett$AICc[i] <- round(MuMIn::AICc(top10[[i]]), 2)
  summary.sett$BIC[i] <- round(BIC(top10[[i]]), 2)
  summary.sett$dev[i] <- round(summary(top10[[i]])$AIC[4], 2)
  summary.sett$Dispersion[i] <- round(summary(top10[[i]])$AIC[4] / df.residual(top10[[i]]), 2)
  summary.sett$mR2[i] <- r.squaredGLMM(top10[[i]])
}
summary.sett <- summary.sett %>% arrange(AICc, BIC)
summary.sett$wAIC <- round(MuMIn::Weights(summary.sett$AICc), 3)
summary.sett$dAIC <- summary.sett$AICc - summary.sett$AICc[1] %>% round(., 3)
summary.sett$dBIC <- summary.sett$BIC - summary.sett$BIC[1] %>% round(., 3)
summary.sett$mR2 <- round(summary.sett$mR2, 3)
print(summary.sett)

# SELECT MODEL
sett <- get.models(sett.select, 1)[[1]] # picked model 1 based on AICc
summary(sett)
model_performance(sett)

# diagnostic plots
print(check_model(sett))

# overdispersion check
# QQ plot
par(mfrow=c(1,2))
plot(predict(sett, type='link', re.form=NA), residuals(sett, type='deviance'), xlab='Predicted', ylab='Deviance')
abline(h=0, lty=2)
qqnorm(residuals(sett, type='deviance'))
qqline(residuals(sett, type='deviance'))
par(mfrow=c(1,1))

# save model outputs
sink(file = 'outputs/settlement_model.txt', append = F) # start file sink

## SETTLEMENT MODEL SELECTION
print(sett.select)
print(summary.sett)

## SELECTED SETTLEMENT MODEL
print(summary(sett))

sink()