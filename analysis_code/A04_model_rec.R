
#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Coral recruitment model analyses (selection and validation)
# Author: Cher Chow

#############################################################################


require(tidyverse)
require(lme4)
require(readxl)
require(MuMIn)
require(performance)
set.seed(24)

# Read in 2018 spat data
settlement.18 <- read_xlsx('./src/coral_settlement.xlsx', sheet='Sheet1', col_names=T)
str(settlement.18)


# Data prep ---------------------------------------------------------------

# select just the 7 study sites and their spat data from 2018
settlement.18 <- settlement.18 %>% group_by(year, site) %>% 
  filter(site %in% c('corner_beach', 'lagoon_1', 'resort', 'southeast', 'turtle_beach', 'vickies', 'north_reef_2'), year=='2018-19')
# north 3 accidentally got labelled as North2

# create a str_replace renaming object so that the site names are consistent with the predictors data
names(site.rename) <- sort(unique(settlement.18$site))
settlement.18$site <- str_replace_all(settlement.18$site, site.rename) # replace
rm(site.rename)

# now consolidate it so that every row is one settlement tile
settlement.18 <- settlement.18 %>% group_by(site) %>% summarise(Spat2018=sum(total))
colnames(settlement.18)[1] <- 'Site'
predictors$Spat2018 <- settlement.18$Spat2018
predictors$ForRate <- as.vector(predictors$ForRate)

# Import recruit data
recruits <- read_xlsx('./src/coral_recruit.xlsx', sheet='data', col_names=T, na='NA')
str(recruits)

recruits <- recruits %>% drop_na() %>% group_by(site, quadrant) %>% summarise(Recruits=sum(Acrorec,Otherrec,Por,Fa,I)) %>% select(-quadrant)
colnames(recruits)[1] <- 'Site'
# rectify some spelling inconsistencies

predictors$Spat2018 <- scale(predictors$Spat2018, scale=T, center=F) %>% as.vector # scale between 0 and 1
recruits$Site <- str_replace_all(recruits$Site, 'Vicky', 'Vickis')
recruits$Site <- str_replace_all(recruits$Site, 'Turtle$', 'TurtleBeach')
RecruitData <- na.omit(recruits) %>% full_join(y=predictors, by='Site') %>% mutate_if(is.character, as.factor) # now join it all

# collinearity check
pairs(RecruitData[,-1])
cor.recruit <- as_tibble(round(cor(RecruitData[,-1]),2))
head(cor.recruit)
cor.recruit$var1 <- colnames(cor.recruit)
cor.recruit <- pivot_longer(as_tibble(cor.recruit), -var1)

ggplot(data = cor.recruit %>% filter(value != 1), aes(x=var1, y=name, fill=sqrt(value^2))) + 
  geom_tile() + scale_fill_viridis_c(option='plasma', name=bquote('absolute value'~R^2)) + labs(x=NULL, y=NULL)

# Construct model candidates and select --------------------------------------------

rec_global <- glmer.nb(data=RecruitData, 
                       Recruits ~ Spat2018 + ForRate + TEve + TDiv + TOP + Herb + Benthic + (1|Site), na.action = 'na.fail')
# use dredge to run model selection with all possible predictor combinations (with site fixed)
rec.select <- dredge(rec_global, beta = 'none', evaluate = T, rank = 'AICc', fixed = c('Spat2018', 'ForRate', 'Site'))
rec.select

## custom model selection table with deviance
top10 <- get.models(rec.select, 1:10)
summary.rec <- data.frame(ModelRank=1:10) # just top 10 candidates
for (i in 1:10) {
  summary.rec$nTerms[i] <- getAllTerms(top10[[i]], intercept = F) %>% length
  summary.rec$AICc[i] <- round(MuMIn::AICc(top10[[i]]), 2)
  summary.rec$BIC[i] <- round(BIC(top10[[i]]), 2)
  summary.rec$dev[i] <- round(summary(top10[[i]])$AIC[4], 2)
  summary.rec$Dispersion[i] <- round(summary(top10[[i]])$AIC[4] / df.residual(top10[[i]]), 2)
}
summary.rec <- summary.rec %>% arrange(AICc, BIC)
summary.rec$wAIC <- round(MuMIn::Weights(summary.rec$AICc), 3)
summary.rec$dAIC <- summary.rec$AICc - summary.rec$AICc[1] %>% round(., 3)
summary.rec$dBIC <- summary.rec$BIC - summary.rec$BIC[1] %>% round(., 3)
print(summary.rec)


recr <- get.models(rec.select, 1)[[1]] # print results for top model
# selected most parsimonious

# Final model diagnostics --------------------------------------------

summary(recr)
model_performance(recr) # AIC and BIC
check_model(recr) # diagnostic plots

# overdispersion check
# QQ plot
par(mfrow=c(1,2))
plot(predict(recr, type='link'), residuals(recr, type='deviance'), xlab='Predicted', ylab='Deviance')
abline(h=0, lty=2)
qqnorm(residuals(recr, type='deviance'))
qqline(residuals(recr, type='deviance'))
par(mfrow=c(1,1))

# save model outputs
sink(file = 'outputs/recruitment_model.txt', append = F) # start file sink

## RECRUITMENT MODEL SELECTION
print(rec.select)
print(summary.rec)

## RECRUITMENT SETTLEMENT MODEL
summary(recr)

sink() # reset back to normal