
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
require(patchwork) # multipanel plotting

# plot aesthetics
quartzFonts(Public=c('Public Sans Regular', 'Public Sans Italic', 'Public Sans Bold', 'Public Sans Bold Italic'))
par(family='Public')
looks <- theme_bw(base_size=13, base_family = 'Public Sans') + theme(panel.grid=element_blank(), axis.ticks=element_line(size=0.3))

# functions
source('convhullvert_function.R') # convex hull vertices for plotting function

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

trait.dis <- cailliez(gowdis(traits[,-1], ord='podani'))
is.euclid(trait.dis) # check whether distances are Euclidean before running PCoA

# run PCoA
pcoa <-  pcoa(trait.dis, rn=colnames(trait.dis)) #PCoA with Cailliez corrections to negative eigenvalues


# PCoA diagnostics --------------------------------------------------------

# Scree plots to check PCoAs, but only 
scree <- ggplot(data=pcoa$values[1:10,], aes(x=1:10, y=Relative_eig/sum(pcoa$values$Relative_eig))) +
    geom_line() +
    geom_point(shape=21, fill='white', size=3) + looks +
    labs(x=NULL, y=NULL) +
    scale_x_continuous(breaks=c(1:10)) +
    theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))

scree

P.dis2 <- cailliez(gowdis(pcoa$vectors[,1:2], ord='podani'))
P.dis4 <- cailliez(gowdis(pcoa$vectors[,3:4], ord='podani'))
# look at the euclidean distances between species in two 2-D trait spaces
  
# compress the distance matrices in P.dis into a dataframe
df_dis <- data.frame(Axis12=as.matrix(P.dis2)[,1], Axis34=as.matrix(P.dis4)[,1], Observed=as.matrix(trait.dis)[,1])
  for (j in 2:nrow(pcoa$vectors)) {
    df_dis <- rbind(df_dis, data.frame(Axis12=as.matrix(P.dis2)[,j], Axis34=as.matrix(P.dis4)[,j], Observed=as.matrix(trait.dis)[,j]))
    df_dis <- df_dis %>% filter(Observed > 0, Axis12 > 0, Axis34 > 0)
  }
  
  # conduct Mantel tests to get correlations between the observed and scaled dissimilarity matrices
  # Pearson correlation for parametric correlation
  P.lm2 <- mantel(trait.dis, P.dis2, method='pearson', permutations=500)
  P.lm4 <- mantel(trait.dis, P.dis4, method='pearson', permutations=500)

# plot all the observed vs scaled dissimilarities with the R2 from the lm
P.stress <- ggplot(df_dis) + 
    geom_point(aes(x=Observed, y=Axis12), shape=21, alpha=0.1) +
    geom_point(aes(x=Observed, y=Axis34), shape=21, alpha=0.1, color='blue') +
    looks + labs(x=NULL, y=NULL) +
    theme(plot.title=element_text(size=11, hjust=0.5, face='bold')) +
    annotate(geom='text', x=min(df_dis$Observed), y=max(df_dis[1]), label=paste0('R=',round(P.lm2$statistic, 3)), size=3) +
    annotate(geom='text', x=min(df_dis$Observed), y=max(df_dis[1])*.9, label=paste0('R=', round(P.lm4$statistic, 3)), color='blue', size=3)


# scaled vs observed dissimilarities in the four dimensional trait space    
P.dis <- cailliez(gowdis(pcoa$vectors[,1:4], ord='podani'))
# look at the euclidean distances between species in 1st 4 dimensions

# compress the distance matrices in P.dis into a dataframe
df_all <- data.frame(Scaled=as.matrix(P.dis)[,1], Observed=as.matrix(trait.dis)[,1])
for (j in 2:dim(pcoa$vectors)[1]) {
  df_all <- rbind(df_all, data.frame(Scaled=as.matrix(P.dis)[,j], Observed=as.matrix(trait.dis)[,j]))
  df_all <- df_all %>% filter(Observed > 0, Scaled > 0)
}
# perform a mantel test to test the correlation of the observed vs 4D dissimilarity matrices
P.lm <- mantel(trait.dis, P.dis, method='pearson', permutations=500)

P.stressall <- ggplot(df_all) + 
  geom_point(aes(x=Observed, y=Scaled), shape=21, alpha=0.1) +
  looks + labs(x=NULL, y=NULL) +
  theme(plot.title=element_text(size=11, hjust=0.5, face='bold')) + 
  annotate(geom='text', x=min(df_all$Observed), y=max(df_all[1]), hjust=1, label=round(P.lm$statistic, 3), size=3)

# print plots
(scree + labs(x='Dimensions', y='Relative eigenvalue')) | (P.stressall + labs(x='Original dissimilarities', y='PCoA dissimilarities')) & labs(title=NULL)

# look at the PcoA points
P.pt1 <- ggplot(as.data.frame(pcoa$vectors), aes(x=Axis.1, y=Axis.2)) +
  geom_point(shape=21) + looks +
  labs(x=NULL, y=NULL) +
  theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))
P.pt2 <- ggplot(as.data.frame(pcoa$vectors), aes(x=Axis.3, y=Axis.4)) +
  geom_point(shape=21) + looks +
  labs(x=NULL, y=NULL) +
  theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))

P.pt1 | P.pt2


## ----FD calculation-------------------------------------------------------------------------------------------------------------

# calculate FD indices for each site with respect to the global functional trait space
# because this is with respect to all study sites, we can standardise FRichness
FD_global <- dbFD(trait.dis, abund, w.abun=T, calc.FRic=T, calc.FDiv=T, m=4, calc.CWM=F, calc.FGR=F, stand.FRic=T, print.pco=T)

# extract the FD indices from the FD_global object, create a dataframe of them
predictors <- data.frame(Site=names(FD_global$FEve[1:7]),FEve=FD_global$FEve[1:7], FDiv=FD_global$FDiv[1:7], FRic=FD_global$FRic[1:7])
# only 7 for the 7 sites. Ignore the global FD measures
rownames(predictors) <- names(FD_global$FEve[1:7])

