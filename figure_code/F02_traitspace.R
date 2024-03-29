
#############################################################################

# FishTraitsCoralRec
# Figure 2: Trait space figures (by site, and also globally showing by each trait)

#############################################################################

require(ggrepel)
require(ggplot2)
require(dplyr)
require(patchwork)
require(viridis)

looks <- theme_bw(base_size=13) + theme(panel.grid=element_blank(), axis.ticks=element_line(size=0.3))

t.space <- as.list(rep(0,4))
source('./analysis_code/function_convhull.R') # make plottable convex hull vertices
source('https://gist.github.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/c9caceef463062c75c71b42607954f7958818ff7/cherulean.R')

hull.v <- convhull.vert(point[,1:2])
hull.v2 <- convhull.vert(point[,3:4])

# palette for sites
palette <- cherulean(7, palette = 'globiceps') %>% paste0(., 'A6')
titles <- c('Corner Beach', 'Lagoon', 'North Reef', 'Resort', 'Southeast', 'Turtle Beach', "Vicki's")

for (i in 1:7) {
  # Now that the base layers are done, I'll plot each site's points and hulls on
  # first, make a dummy points dataframe only with the species in that site
  pointS <- fish_assemblage %>% filter(Site == unique(fish_assemblage$Site)[i]) %>% ungroup() %>% select(Species, n)
  pointS$n <- pointS$n/site.total$site.total[i] # relative abundance of each species
  pointS <- left_join(pointS, point, by='Species')
  pointS <- pointS[order(pointS$n, decreasing=T),order(colnames(pointS))]
  # make abbreviated species names
  shortsp <- as.data.frame(matrix(unlist(strsplit(pointS$Species[1:3], " ")), ncol=2, byrow=TRUE))
  pointS$ShortSp <- ''
  pointS$ShortSp[1:3] <- paste(abbreviate(shortsp[,1], 1, strict = TRUE), ". ", shortsp[,2], sep="")
  
  # now calculate the convex hull for that site only
  # global convex hull/TOP
  s.hull.v <- convhull.vert(pointS[,1:2])
  s.hull.v2 <- convhull.vert(pointS[,3:4])
  
  s.space <- as.list(c(0,0))
  
  s.space[[1]] <- ggplot() + looks +
    geom_polygon(data=hull.v, aes(x=A1, y=A2), color='grey', fill='white') +
    labs(x='PCo1', y='PCo2', title=titles[i]) +
    geom_polygon(data=s.hull.v, aes(x=A1, y=A2), fill=palette[i], alpha=0.5) +
    geom_point(data=pointS %>% left_join(., traits, by="Species"),
               aes(x=A1, y=A2, size=n), fill=palette[i], color='#666666', shape=21) +
    geom_text_repel(data=pointS %>% arrange(n) %>% top_n(5), aes(label=ShortSp, x=A1, y=A2), force=2,
                    size=2.8, fontface='italic', min.segment.length = 0, point.size=5,
                    max.overlaps = Inf, box.padding = 1, force_pull=0.1) +
    scale_size_continuous(range=c(1.5,8), guide="none", limits=c(0,0.6)) +
    if (i != 4) {
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), 
            legend.position='none')
    } 
  else {theme()}
  s.space[[2]] <- ggplot() + looks +
    geom_polygon(data=hull.v2, aes(x=A3, y=A4), color='grey', fill='white') +
    labs(x='PCo3', y='PCo4') +
    geom_polygon(data=s.hull.v2, aes(x=A3, y=A4), fill=palette[i], alpha=0.5) +
    geom_point(data=pointS %>% left_join(., traits, by="Species"), 
               aes(x=A3, y=A4, size=n), fill=palette[i], color='#666666', shape=21) +
    geom_text_repel(data=pointS %>% arrange(n) %>% top_n(5), aes(label=ShortSp, x=A3, y=A4), force=2,
                    size=2.8, fontface='italic', min.segment.length = 0, point.size=5,
                    max.overlaps = Inf, box.padding = 1, force_pull=0.1) +
    scale_size_continuous(range=c(1.5,8), limits=c(0,0.6), name="Relative abundance")  + 
    if (i != 4) {
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), 
            legend.position='none')
    } 
  else {theme()}
  
  t.space[[i]] <- s.space[[1]] / s.space[[2]] # stack them with patchwork and store
  rm(pointS)
}




Sp.count <- fish_assemblage %>% group_by(Species) %>% summarise(n = sum(n)) %>% ungroup
point$n <- Sp.count$n/sum(Sp.count$n) # put in global abundances to visualise all sites
global1 <- ggplot() + looks +
  geom_polygon(data=hull.v, aes(x=A1, y=A2), fill='transparent', color='grey') +
  geom_point(data=left_join(point, traits, by="Species"), aes(x=A1, y=A2, shape=ForageMode), size=2, fill='black', alpha=0.5) +
  labs(x='PCo1', y='PCo2', title='All sites') +
  scale_fill_grey(aesthetics="color",start=0.2, end=0.7) + scale_shape_manual(values=c(1:6,21:25), name="Foraging mode")

global2 <- ggplot() + looks +
  geom_polygon(data=hull.v2, aes(x=A3, y=A4), fill='transparent', color='grey') +
  geom_point(data=left_join(point, traits, by="Species"), aes(x=A3, y=A4, shape=ForageMode), size=2, fill='black', alpha=0.5) +
  labs(x='PCo3', y='PCo4') + theme(legend.title=element_blank()) +
  scale_fill_grey(aesthetics="color", start=0.2, end=0.7) + scale_shape_manual(values=c(1:6,21:25), guide="none")

pred_long <- predictors %>% pivot_longer(cols=!Site, names_to="predictor", values_to="value")

FD_bar <- ggplot(pred_long %>% filter(predictor %in% c('TOP', 'TEve', 'TDiv'))) +
  geom_bar(aes(y=value, x=predictor, fill=Site), stat='identity', position='dodge', color='black') +
  labs(y='Index measure', x=NULL) +
  scale_fill_cherulean(discrete = T, palette = 'globiceps', guide="none") + looks +
  scale_y_continuous(expand=expansion(mult=c(.0,.05)), limits=c(0,1))

A_bar <- ggplot(pred_long %>% filter(predictor == 'Herb' | predictor == 'Benthic')) +
  geom_bar(aes(y=value, x=predictor, group=Site, fill=Site), stat='identity', position='dodge', color='black') +
  labs(y='Relative abundance', x=NULL) +
  scale_fill_cherulean(discrete = T, palette = 'globiceps', guide="none") + looks + 
  scale_y_continuous(expand=expansion(mult=c(.0,.05))) +
  scale_x_discrete(labels=c('Sessile invertivores', 'Herbivores'))

bars <- (A_bar | FD_bar) * plot_layout(widths=c(2,3))

Fig2 <- (bars / ((global1/global2) | t.space[[1]] | t.space[[2]] | t.space[[3]]) / (t.space[[4]] | t.space[[5]] | t.space[[6]] | t.space[[7]])) * plot_layout(guides='collect', heights = c(1,3,3)) * theme(axis.text = element_text(size=9)) & theme(plot.title=element_text(face='bold', hjust=0.5, size=13))
Fig2

# ggsave(plot = Fig2, filename = paste0(fig_dir, "02_traitspaces.svg"), device = "svg", width=30, height=30, units='cm', dpi=300)
rm(FD_bar, A_bar, bars, global1, global2, t.space, s.hull.v, s.hull.v2, s.space, hull.v, hull.v2)
point$n <- NULL # necessary for rerunning
