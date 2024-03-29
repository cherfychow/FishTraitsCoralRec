
#############################################################################

# FishTraitsCoralRec
# Trait space quality diagnostic plots and checks

#############################################################################

set.seed(24)
require(tidyverse)
require(patchwork)
require(geometry) # calculates convex hull/trait space for data vis

# PCoA diagnostics --------------------------------------------------------
# plot aesthetics
looks <- theme_bw(base_size=13) + theme(panel.grid=element_blank(), axis.ticks=element_line(size=0.3))
trait_pcoa <- pcoa(trait.dis, correction = 'none')
traits <- traits %>% mutate(across(where(is.character), as.factor))
trait_fit <- envfit(trait_pcoa$vectors, traits[-1], permutations = 999, choices = 1:4)
plot(x = point$A1, y = point$A2, type = 'p') # just the points
plot(trait_fit, choices = c(1,2), axis = FALSE) # with the vectors

# PCo 3 + 4
plot(x = point$A3, y = point$A4, type = 'p') # just the points
plot(trait_fit, choices = c(3,4), axis = FALSE) # with the vectors

# Scree plots to check PCoAs, but only 
scree <- ggplot(data=trait_pcoa$values[1:10,], aes(x=1:10, y=Relative_eig/sum(trait_pcoa$values$Relative_eig))) +
  geom_line() +
  geom_point(shape=21, fill='white', size=3) + looks +
  labs(x=NULL, y=NULL) +
  scale_x_continuous(breaks=c(1:10)) +
  theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))
scree

P.dis2 <- cailliez(gowdis(trait_pcoa$vectors[,1:2], ord='podani'))
P.dis4 <- cailliez(gowdis(trait_pcoa$vectors[,3:4], ord='podani'))
# look at the euclidean distances between species in two 2-D trait spaces

# compress the distance matrices in P.dis into a dataframe
df_dis <- data.frame(Axis12=as.matrix(P.dis2)[,1], Axis34=as.matrix(P.dis4)[,1], Observed=as.matrix(trait.dis)[,1])
for (j in 2:nrow(trait_pcoa$vectors)) {
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
P.dis <- cailliez(gowdis(trait_pcoa$vectors[,1:4], ord='podani'))
# look at the euclidean distances between species in 1st 4 dimensions

# compress the distance matrices in P.dis into a dataframe
df_all <- data.frame(Scaled=as.matrix(P.dis)[,1], Observed=as.matrix(trait.dis)[,1])
for (j in 2:dim(trait_pcoa$vectors)[1]) {
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
FigS1 <- (scree + labs(x='Dimensions', y='Relative eigenvalue')) | (P.stressall + labs(x='Original dissimilarities', y='PCoA dissimilarities')) & labs(title=NULL)
FigS1
# ggsave(filename="./figures/S1_pcoa_traitchange.svg", device='svg', width=25, height=11, units='cm', dpi=300)

# look at the PcoA points
P.pt1 <- ggplot(as.data.frame(trait_pcoa$vectors), aes(x=Axis.1, y=Axis.2)) +
  geom_point(shape=21) + looks +
  labs(x=NULL, y=NULL) +
  theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))
P.pt2 <- ggplot(as.data.frame(trait_pcoa$vectors), aes(x=Axis.3, y=Axis.4)) +
  geom_point(shape=21) + looks +
  labs(x=NULL, y=NULL) +
  theme(plot.title=element_text(size=11, hjust=0.5, face='bold'))

P.pt1 | P.pt2

rm(P.stress, P.pt1, P.pt2, P.dis2, P.dis4, 
   P.dis, P.lm2, P.lm4, df_dis, trait_fit) # get rid of the extraneous objects not needed for FigS1
