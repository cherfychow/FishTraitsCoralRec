
#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Generate trait space for trait diversity analysis
# This file also contains diagnostic plots for the trait space analysis
# Author: Cher Chow
# Updated: 26 Mar 2021

#############################################################################

require(tidyverse)
require(FD) # calculate dissimilarity matrices + functional diversity indices
require(geometry) # calculates convex hull/trait space for data vis

set.seed(24) # repeatability :)

# Pairwise trait check, just so they're not directly colinear
pairs(traits[,-1], labels=colnames(traits)[-1], col=rgb(0.2,0.2,0.2, 0.5))


## Trait space ------------------------------------------------------------------------------------------------------------

# create the abundance matrix for functional evenness and divergence weighting
# species are columns, each site + total = rows
abund <- Sp.count %>% filter(str_detect(SpeciesV, ' sp$', negate=T)) %>% select(SpeciesV, n) %>% spread(SpeciesV, n) # row 1 = all sites
abund$Site <- 'All'
abund <- rbind(UBRUVspecies %>% select(Site, Species, n) %>% spread(Species, n), abund)
# because not all species occur in each site, there are lots of NAs to be filled with 0s
NAlist <- as.list(rep('0',dim(abund)[2]-1)) # list has to have column names for every item = 0
names(NAlist) <- colnames(abund)[2:ncol(abund)]
abund <- abund %>% replace_na(NAlist) # annnd replace NAs
rownames(traits) <- traits$Species # dbFD function requires species in the row names

# abundance data frame is in characters????
abund <- abund %>% mutate(across(where(is.character), as.numeric))
abund <- as.data.frame(abund)
rownames(abund) <- abund$Site
abund$Site <- NULL
rm(NAlist)

## ----Trait space sites, results='hide'------------------------------------------------------------------------------------------
# running the full lizard island traits together is very computationally taxing
# so I'll split all the sites to run the FD package separately
abun <- as.list(rep('0',7))
for (i in 1:7) {
  abun[[i]] <- as.data.frame(UBRUV %>% filter(Site == unique(UBRUV$Site)[i]) %>% ungroup() %>% select(Species, n) %>% spread(Species, n))
}
names(abun) <- unique(UBRUV$Site)
abun[[8]] <- Sp.count %>% filter(str_detect(Species, ' sp$', negate=T)) %>% select(Species, n) %>% spread(Species, n)
names(abun)[8] <- 'all' # add one item for all sites

# also make separate trait tables for each site
traits.sites <- as.list(rep('0', 7))
for (i in 1:7) {
  site.sp <- UBRUV %>% filter(Site == unique(UBRUV$Site)[i]) %>% ungroup() %>% select(Species)
  # filter the tibble by each site and extract the species
  traits.sites[[i]] <- as.data.frame(left_join(site.sp, as_tibble(traits), by='Species'))
  # Use the site specific species list to match with the traits
  rownames(traits.sites[[i]]) <- traits.sites[[i]]$Species # row names bc FD requires it
  traits.sites[[i]][,1] <- NULL # Remove the Species column now that things are joined
}
names(traits.sites) <- unique(UBRUV$Site) # name list items by site
rm(site.sp)
traits.sites[[8]] <- traits[,-1]
names(traits.sites)[8] <- 'all' # all sites again

# create the distance matrix for each site
# Gower distance, Podani method for ordered factors.
trait.dis <- as.list(rep('0', 8))
discheck <- rep('0', 8)
for (i in 1:8) {
  # calculate gower dissimilarities on standardised continuous traits
  # and apply Cailliez correction to make it Euclidean
  trait.dis[[i]] <- cailliez(gowdis(traits.sites[[i]], ord='podani'))
  discheck[i] <- is.euclid(trait.dis[[i]])
}

# check the distance matrices. The FD package will check this too, but good to check before PCoA
discheck
rm(discheck)

# run PCoA and NMDS for all sites, using both because PCoA is returning a lot of dimensions.
PCoA <- as.list(rep('0', 8))
NMDS <- as.list(rep('0', 8))
for (i in 1:8) {
  PCoA[[i]] <-  pcoa(trait.dis[[i]], rn=colnames(trait.dis[[i]])) #PCoA with Cailliez corrections to negative eigenvalues
  NMDS[[i]] <- vegan::metaMDS(trait.dis[[i]], distance='gower', k=4, wascores=T)
}


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------
# Scree plots to check PCoAs, but only 
scree <- as.list(rep('0', 8))
for (i in 1:8) {
  scree[[i]] <- ggplot(data=PCoA[[i]]$values[1:10,], aes(x=1:10, y=Relative_eig/sum(PCoA[[i]]$values$Relative_eig))) +
    geom_line() +
    geom_point(shape=21, fill='white', size=3) + looks +
    labs(x=NULL, y=NULL) +
    scale_x_continuous(breaks=c(1:10)) +
    theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))
}

scree[[1]] + scree[[2]] + scree[[3]] + scree[[4]] + scree[[5]] + scree[[6]] + scree[[7]] + scree[[8]]



## ----Scaled vs observed distance, separate spaces-------------------------------------------------------------------------------
P.dis2 <- as.list(rep(0, 8))
P.dis4 <- as.list(rep(0, 8))
P.stress <- as.list(rep(0, 8))
df_dis <- as.list(rep(0, 8))
P.lm2 <- as.list(rep(0, 8))
P.lm4 <- as.list(rep(0, 8))
for (i in 1:8) {
  P.dis2[[i]] <- cailliez(gowdis(PCoA[[i]]$vectors[,1:2], ord='podani'))
  P.dis4[[i]] <- cailliez(gowdis(PCoA[[i]]$vectors[,3:4], ord='podani'))
  # look at the euclidean distances between species in two 2-D trait spaces
  
  # compress the distance matrices in P.dis into a dataframe
  df_dis[[i]] <- data.frame(Axis12=as.matrix(P.dis2[[i]])[,1], Axis34=as.matrix(P.dis4[[i]])[,1], Observed=as.matrix(trait.dis[[i]])[,1])
  for (j in 2:dim(PCoA[[i]]$vectors)[1]) {
    df_dis[[i]] <- rbind(df_dis[[i]], data.frame(Axis12=as.matrix(P.dis2[[i]])[,j], Axis34=as.matrix(P.dis4[[i]])[,j], Observed=as.matrix(trait.dis[[i]])[,j]))
    df_dis[[i]] <- df_dis[[i]] %>% filter(Observed > 0, Axis12 > 0, Axis34 > 0)
  }
  
  # conduct Mantel tests to get correlations between the observed and scaled dissimilarity matrices
  # Pearson correlation for parametric correlation
  P.lm2[[i]] <- mantel(trait.dis[[i]], P.dis2[[i]], method='pearson', permutations=500)
  P.lm4[[i]] <- mantel(trait.dis[[i]], P.dis4[[i]], method='pearson', permutations=500)
}

P.stress[[9]] <- ggplot(data=data.frame(x=1:9, y=1:9), aes(x, y)) +
  annotate(geom='text', x=2, y=1:8, hjust=0, label=names(traits.sites), size=3) + looks +
  theme(axis.ticks=element_blank(), axis.text=element_blank()) + labs(x=NULL, y=NULL) + scale_y_continuous(limits=c(0,9))
# plot all the observed vs scaled dissimilarities with the R2 from the lm we did for every site.
for (i in 1:8) {
  P.stress[[i]] <- ggplot(df_dis[[i]]) + 
    geom_point(aes(x=Observed, y=Axis12), shape=21, alpha=0.1) +
    geom_point(aes(x=Observed, y=Axis34), shape=21, alpha=0.1, color='blue') +
    looks + labs(x=NULL, y=NULL, title=names(traits.sites)[i]) +
    theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))
  P.stress[[9]] <- P.stress[[9]] + 
    annotate(geom='text', x=5, y=i, hjust=0, label=round(P.lm2[[i]]$statistic, 3), size=3) +
    annotate(geom='text', x=7, y=i, hjust=1, label=round(P.lm4[[i]]$statistic, 3), color='blue', size=3)
}

# scaled vs observed dissimilarities in the four dimensional trait space    
P.dis <- as.list(rep(0, 8))
P.stressall <- as.list(rep(0, 8))
P.lm <- as.list(rep(0,8))
df_all <- as.list(rep(0,8))
for (i in 1:8) {
  P.dis[[i]] <- cailliez(gowdis(PCoA[[i]]$vectors[,1:4], ord='podani'))
  # look at the euclidean distances between species in 1st 4 dimensions
  
  # compress the distance matrices in P.dis into a dataframe
  df_all[[i]] <- data.frame(Scaled=as.matrix(P.dis[[i]])[,1], Observed=as.matrix(trait.dis[[i]])[,1])
  for (j in 2:dim(PCoA[[i]]$vectors)[1]) {
    df_all[[i]] <- rbind(df_all[[i]], data.frame(Scaled=as.matrix(P.dis[[i]])[,j], Observed=as.matrix(trait.dis[[i]])[,j]))
    df_all[[i]] <- df_all[[i]] %>% filter(Observed > 0, Scaled > 0)
  }
  # perform a mantel test to test the correlation of the observed vs 4D dissimilarity matrices
  P.lm[[i]] <- mantel(trait.dis[[i]], P.dis[[i]], method='pearson', permutations=500)
}

P.stressall[[9]] <- ggplot(data=data.frame(x=1:9, y=1:9), aes(x, y)) +
  annotate(geom='text', x=2, y=1:8, hjust=0, label=names(traits.sites), size=3) + looks +
  theme(axis.ticks=element_blank(), axis.text=element_blank()) + labs(x=NULL, y=NULL) + scale_y_continuous(limits=c(0,9))

for (i in 1:8) {
  # plot all the observed vs scaled dissimilarities in the matrices for every site
  P.stressall[[i]] <- ggplot(df_all[[i]]) + 
    geom_point(aes(x=Observed, y=Scaled), shape=21, alpha=0.1) +
    looks + labs(x=NULL, y=NULL, title=names(traits.sites)[i]) +
    theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))
  P.stressall[[9]] <- P.stressall[[9]] + 
    annotate(geom='text', x=5, y=i, hjust=1, label=round(P.lm[[i]]$statistic, 3), size=3)
}

# print plots

P.stressall[[1]] + P.stressall[[2]] + P.stressall[[3]] + P.stressall[[4]] + 
  P.stressall[[5]] + P.stressall[[6]] + P.stressall[[7]] + P.stressall[[8]] + P.stressall[[9]]
(scree[[8]] + labs(x='Dimensions', y='Relative eigenvalue')) | (P.stressall[[8]] + labs(x='Original dissimilarities', y='PCoA dissimilarities')) & labs(title=NULL)


## ----correlation between eigenvalues + dissimilarities--------------------------------------------------------------------------
EigMant <- data.frame(SumEig=sum(PCoA[[1]]$values$Relative_eig[1:4]), Mantel=P.lm[[1]]$statistic, sig=P.lm[[1]]$signif)
for (i in 2:8) {
  EigMant <- rbind(EigMant, data.frame(SumEig=sum(PCoA[[i]]$values$Relative_eig[1:4]), Mantel=P.lm[[i]]$statistic, sig=P.lm[[i]]$signif))
}
EigMant.lm <- lm(data=EigMant, Mantel~SumEig)
summary(EigMant.lm)

library(ggrepel)
ggplot(EigMant, aes(x=SumEig, y=Mantel)) +
  geom_smooth(method='lm', se=TRUE, color='black', size=0.5) +
  geom_point(shape=21, size=2.5) + 
  geom_text_repel(label=rownames(abund)) + looks


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------
P.pt1 <- as.list(rep(0, 8))
P.pt2 <- as.list(rep(0, 8))
for (i in 1:8) {
  P.pt1[[i]] <- ggplot(as.data.frame(PCoA[[i]]$vectors), aes(x=Axis.1, y=Axis.2)) +
    geom_point(shape=21) + looks +
    labs(x=NULL, y=NULL, title=names(traits.sites)[i]) +
    theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))
  P.pt2[[i]] <- ggplot(as.data.frame(PCoA[[i]]$vectors), aes(x=Axis.3, y=Axis.4)) +
    geom_point(shape=21) + looks +
    labs(x=NULL, y=NULL, title=names(traits.sites)[i]) +
    theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))
}
(P.pt1[[1]] | P.pt1[[2]] | P.pt1[[3]] | P.pt1[[4]]) / (P.pt2[[1]] | P.pt2[[2]] | P.pt2[[3]] | P.pt2[[4]])
(P.pt1[[5]] | P.pt1[[6]] | P.pt1[[7]] | P.pt1[[8]])/(P.pt2[[5]] | P.pt2[[6]] | P.pt2[[7]] | P.pt2[[8]])


## ----NMDS stress, echo=FALSE----------------------------------------------------------------------------------------------------
par(mfrow=c(2,2))
for (i in 1:8) {
  stressplot(NMDS[[i]], p.col='grey', l.col='black')
  title(main=names(traits.sites)[i])
} # check stress and fit
par(mfrow=c(1,1))
par(mfrow=c(2,4))
for (i in 1:4) {
  with(as.data.frame(NMDS[[i]]$points), plot(x=MDS1, y=MDS2))
  title(main=names(traits.sites)[i], sub=paste('stress =', round(NMDS[[i]]$stress, 3), sep=' '))
}
for(i in 1:4) {
  with(as.data.frame(NMDS[[i]]$points), plot(x=MDS3, y=MDS4))
}
for (i in 5:8) {
  with(as.data.frame(NMDS[[i]]$points), plot(x=MDS1, y=MDS2))
  title(main=names(traits.sites)[i], sub=paste('stress =', round(NMDS[[i]]$stress, 3), sep=' '))
}
for(i in 5:8) {
  with(as.data.frame(NMDS[[i]]$points), plot(x=MDS3, y=MDS4))
}
par(mfrow=c(1,1))


## ----FD calculation-------------------------------------------------------------------------------------------------------------
# calculate FD indices for each site with respect to each site's functional trait space
FD <- as.list(rep(0,8))
for (i in 1:8) {
  FD[[i]] <- dbFD(trait.dis[[i]], abun[[i]], w.abun=T, calc.FRic=T, calc.FDiv=T, m=4, calc.CWM=F, calc.FGR=F, print.pco=T)
  # not standardising FRic because there is no global FRic measure by running each site separately
}
names(FD) <- names(abun)[[8]]

# calculate FD indices for each site with respect to the global functional trait space
# because this is with respect to all study sites, we can standardise FRichness
FD_global <- dbFD(trait.dis[[8]], abund, w.abun=T, calc.FRic=T, calc.FDiv=T, m=4, calc.CWM=F, calc.FGR=F, stand.FRic=T, print.pco=T)


## ----Indices to dataframe-------------------------------------------------------------------------------------------------------
# extract the FD indices from the FD_global object, create a dataframe of them
predictors <- data.frame(Site=names(FD_global$FEve[1:7]),FEve=FD_global$FEve[1:7], FDiv=FD_global$FDiv[1:7], FRic=FD_global$FRic[1:7])
# only 1:7 to ignore the global measures
rownames(predictors) <- names(FD_global$FEve[1:7])