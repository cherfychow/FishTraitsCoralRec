
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
settlement.18 <- read_xlsx('./src/coral_settlement.xlsx', sheet='Sheet1', col_names=T)
str(settlement.18)

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

# Import recruit data
recruits <- read_xlsx('./src/coral_recruit.xlsx', sheet='data', col_names=T, na='NA')
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
  geom_tile() + scale_fill_viridis_c(option='plasma', name=bquote('absolute value'~R^2)) + labs(x=NULL, y=NULL)

# Construct model candidates --------------------------------------------

Model.Rec <- as.list(rep(0,10))
names(Model.Rec) <- paste0('Rec', 1:10)

Model.Rec[[1]] <- glmer.nb(data=RecruitData,
                              formula=Recruits ~ Spat2018 + (1|Site)) # null model. Failed to converge
Model.Rec[[2]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + TEve + TDiv + TOP + Herb + Benthic + (1|Site))
Model.Rec[[3]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + Herb + Benthic + (1|Site))
Model.Rec[[4]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + TEve + TDiv + TOP + (1|Site))
Model.Rec[[5]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + TDiv + TOP + Herb + (1|Site))
Model.Rec[[6]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + TDiv + TOP + Herb + Benthic + (1|Site))
Model.Rec[[7]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + TDiv + Herb + (1|Site))
Model.Rec[[8]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + TEve + TDiv + TOP + Benthic + (1|Site))
Model.Rec[[9]] <- glmer.nb(data=RecruitData, 
                              Recruits ~ Spat2018 + ForRate + Herb + (1|Site))
Model.Rec[[10]] <- glmer.nb(data=RecruitData,
                               formula=Recruits ~ (1|Site)) # null model 2, failed to converge

for (i in 1:length(Model.Rec)) {
  print(summary(Model.Rec[[i]]))
}

Model.Rec <- Model.Rec[-c(1,10)] # remove the ones that fail to converge

# Recruitment model comparison --------------------------------------------
summary.rec <- data.frame(Model=c(1:8))
for (i in 1:8) {
  summary.rec$nTerms[i] <- length(rownames(summary(Model.Rec[[i]])$coefficients))-1
  summary.rec$AICc[i] <- round(MuMIn::AICc(Model.Rec[[i]]), 2)
  summary.rec$BIC[i] <- round(BIC(Model.Rec[[i]]), 2)
  summary.rec$dev[i] <- round(summary(Model.Rec[[i]])$AIC[4], 2)
  summary.rec$Dispersion[i] <- round(deviance(Model.Rec[[i]])/summary(Model.Rec[[i]])$AIC[5], 2)
}
summary.rec <- summary.rec %>% arrange(AICc, BIC)
rownames(summary.rec) <- summary.rec$Model
summary.rec$wAIC <- round(MuMIn::Weights(summary.rec$AICc), 3)
summary.rec$dAIC <- summary.rec$AICc - summary.rec$AICc[1] %>% round(., 3)
summary.rec$dBIC <- summary.rec$BIC - summary.rec$BIC[1] %>% round(., 3)

summary.rec

# Final model diagnostics --------------------------------------------
i=8
model_performance(Model.Rec[[i]])
check_model(Model.Rec[[i]])
# diagnostic plots
#print(check_model(Model.Rec[[i]]))
# overdispersion check
# QQ plot
par(mfrow=c(1,2))
plot(predict(Model.Rec[[i]], type='link'), residuals(Model.Rec[[i]], type='deviance'), xlab='Predicted', ylab='Deviance')
abline(h=0, lty=2)
qqnorm(residuals(Model.Rec[[i]], type='deviance'))
qqline(residuals(Model.Rec[[i]], type='deviance'))
par(mfrow=c(1,1))