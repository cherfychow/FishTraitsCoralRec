
#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Coral recruitment model analyses (selection and validation)
# Author: Cher Chow

#############################################################################


require(tidyverse)
require(lme4)
require(readxl)
set.seed(24)

# Read in 2018 spat data
settlement <- read_xlsx('./src/LIRS_settlement_2019-20.xlsx', sheet='Sheet1', col_names=T)
head(settlement)
str(settlement)

# select just the 7 study sites and their spat data from 2018
sett.trim <- settlement %>% group_by(year, site) %>% 
  filter(str_detect(site, 'corner_beach|lagoon_1|resort|southeast|turtle_beach|vickies|north_reef_2'), year=='2018-19')
# north 3 accidentally got labelled as North2

# create a str_replace renaming object so that the site names are consistent with the predictors data
site.rename <- sort(predictors$Site)
names(site.rename) <- sort(unique(sett.trim$site))
sett.trim$site <- str_replace_all(sett.trim$site, site.rename) # replace
rm(site.rename)

# now consolidate it so that every row is one settlement tile
settlement <- sett.trim %>% group_by(site) %>% summarise(Spat2018=sum(total))
colnames(settlement)[1] <- 'Site'
predictors$Spat2018 <- settlement$Spat2018

# Import recruit data
recruits <- read_xlsx('../DataCoral/Recruits2019_cleaned.xlsx', sheet='data', col_names=T, na='NA')
str(recruits)

recruits <- recruits %>% drop_na() %>% group_by(site, quadrant) %>% summarise(Recruits=sum(Acrorec,Otherrec,Por,Fa,I)) %>% select(-quadrant)
colnames(recruits)[1] <- 'Site'
# rectify some spelling inconsistencies

predictors$Spat2018 <- scale(predictors$Spat2018, scale=T, center=F)
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
  geom_tile() + looks + scale_fill_viridis_c(option='plasma', name=bquote('absolute value'~R^2)) + labs(x=NULL, y=NULL)

## ----Construct recruit models------------------------------------------------------------------------------------------------
RecModels <- as.list(rep(0,10))
names(RecModels) <- paste0('Rec', 1:10)

RecModels[[1]] <- glmer.nb(data=RecruitData,
                              formula=Recruits ~ Spat2018 + (1|Site)) # null model. Failed to converge
RecModels[[2]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + TEve + TDiv + TOP + Herb + Benthic + (1|Site))
RecModels[[3]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + Herb + Benthic + (1|Site))
RecModels[[4]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + TEve + TDiv + TOP + (1|Site))
RecModels[[5]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + TDiv + TOP + Herb + (1|Site))
RecModels[[6]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + TDiv + TOP + Herb + Benthic + (1|Site))
RecModels[[7]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + TDiv + Herb + (1|Site))
RecModels[[8]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + TEve + TDiv + TOP + Benthic + (1|Site))
RecModels[[9]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + Herb + (1|Site))
RecModels[[10]] <- glmer.nb(data=RecruitData,
                               formula=Recruits ~ (1|Site)) # null model 2, failed to converge

for (i in 1:length(RecModels)) {
  print(summary(RecModels[[i]]))
}

RecModels <- RecModels[-c(1,10)]

## ----NB model comparison summary, echo=F----------------------------------------------------------------------------------------
rec.summ <- data.frame(Model=c(1:8))
for (i in 1:8) {
  rec.summ$nTerms[i] <- length(rownames(summary(RecModels[[i]])$coefficients))-1
  rec.summ$AICc[i] <- round(MuMIn::AICc(RecModels[[i]]), 2)
  rec.summ$mR2[i] <- round(as.numeric(r2(RecModels[[i]])[2]),3)
  rec.summ$Dispersion[i] <- round(deviance(RecModels[[i]])/summary(RecModels[[i]])$AIC[5], 2)
}
rec.summ <- rec.summ %>% arrange(AICc)
rownames(rec.summ) <- rec.summ$Model
rec.summ$wAIC <- round(MuMIn::Weights(rec.summ$AICc), 3)
rec.summ$dAIC <- rec.summ$AICc - rec.summ$AICc[1] %>% round(., 3)

rec.summ

## ----R4 diagnostic plots, echo=F---------------------------------------------------------------------------------------------
i=8
model_performance(RecModels[[i]])
check_model(RecModels[[i]])
# diagnostic plots
#print(check_model(RecModels[[i]]))
# overdispersion check
# QQ plot
par(mfrow=c(1,2))
plot(predict(RecModels[[i]], type='link'), residuals(RecModels[[i]], type='deviance'), xlab='Predicted', ylab='Deviance')
abline(h=0, lty=2)
qqnorm(residuals(RecModels[[i]], type='deviance'))
qqline(residuals(RecModels[[i]], type='deviance'))
par(mfrow=c(1,1))