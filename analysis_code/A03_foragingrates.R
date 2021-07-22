#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Trait weighted foraging rates analysis
# Author: Cher Chow

#############################################################################


set.seed(24)

influence <- read.csv('./src/trait_data_influence.csv', header=T) # the modified trait influence table
str(influence)
BiteSp <- unique(UBRUVbite$Species)
influence <- influence %>% filter(str_detect(Species, str_c(BiteSp, collapse='|')))
inf.pca <- princomp(influence[,-1])
loadings(inf.pca)
summary(inf.pca)
screeplot(inf.pca, type='lines') # yes 

sp.inf <-  data.frame(Species=influence$Species, SpInf=inf.pca$scores[,1])
sp.inf$SpInf <- (sp.inf$SpInf + (-range(sp.inf$SpInf)[1]) + 0.1)/2

ForRate <- UBRUV.L %>% ungroup() %>% group_by(Site, Species, Length) %>% 
  select(Site, Species, Length, BiteRate) %>%  
  mutate(LBi=as.numeric(Length)*BiteRate) %>%  # create the middle variable for the weighted foraging rates
  ungroup()
ForRate <- ForRate %>% group_by(Site, Species) %>% 
  filter(BiteRate > 0) %>% 
  summarise(sumLBi=sum(LBi)) # aggregate records by Site and Species

# Now use the calculated species influence weighting factors to adjust bite rates to foraging rates
predictors$ForRate <- left_join(ForRate, sp.inf, by='Species') %>% 
  mutate(SpFeed=0.01*SpInf*sumLBi) %>% group_by(Site) %>% 
  summarise(ForRate=sum(SpFeed)) %>% ungroup() %>% pull(ForRate)
predictors$ForRate <- scale(predictors$ForRate, scale=T, center=F)

# foraging rates by foraging mode
feed <- left_join(ForRate, sp.inf, by='Species') %>% mutate(SpFeed=0.01*SpInf*sumLBi) %>% 
  left_join(., traits %>% select(Species, TrophicGroup, ForageMode), by="Species")
feedF <- feed %>% group_by(Site,TrophicGroup, ForageMode) %>% summarise(Feed=sum(SpFeed), Infl=mean(SpInf))
feedF$ForageMode <- factor(feedF$ForageMode, levels=feedF %>% ungroup() %>% arrange(Infl) %>% distinct(ForageMode) %>% pull)
# foraging rates by trophic group
feedT <- feed %>% group_by(Site,TrophicGroup) %>% summarise(Feed=sum(SpFeed), Infl=mean(SpInf))
feedT$TrophicGroup <- factor(feedT$TrophicGroup, levels=feedT %>% ungroup() %>% arrange(Infl) %>% distinct(TrophicGroup) %>% pull)