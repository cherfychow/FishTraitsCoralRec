
#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Figure 2: Trait space figures (by site, and also globally showing by each trait)
# Author: Cher Chow

#############################################################################

t.space <- as.list(rep(0,4))
point <- as.data.frame(FD_global$x.axes)[,1:4]
point$Species <- rownames(traits)

source('./analysis_code/function_convhull.R') # make plottable convex hull vertices

hull.v <- convhull.vert(point[,1:2])
hull.v2 <- convhull.vert(point[,3:4])

# palette for sites
palette <- viridis::viridis(n=7, end=0.95, begin=0.1, alpha=0.7, option="mako")
titles <- c('Corner Beach', 'Lagoon', 'North Reef', 'Resort', 'Southeast', 'Turtle Beach', "Vicki's")
require(ggrepel)
for (i in 1:7) {
  # Now that the base layers are done, I'll plot each site's points and hulls on
  # first, make a dummy points dataframe only with the species in that site
  s.point <- UBRUVspecies %>% filter(Site == unique(UBRUVspecies$Site)[i]) %>% ungroup() %>% select(Species, n)
  s.point$n <- s.point$n/site.total$site.total[i]
  s.point <- left_join(s.point, point, by='Species')
  s.point <- s.point[order(s.point$n, decreasing=T),order(colnames(s.point))]
  # make abbreviated species names
  shortsp <- as.data.frame(matrix(unlist(strsplit(s.point$Species[1:3], " ")), ncol=2, byrow=TRUE))
  s.point$ShortSp <- ''
  s.point$ShortSp[1:3] <- paste(abbreviate(shortsp[,1], 1, strict = TRUE), ". ", shortsp[,2], sep="")
  
  # now calculate the convex hull for that site only
  # global convex hull/FRic
  s.hull.v <- convhull.vert(s.point[,1:2])
  s.hull.v2 <- convhull.vert(s.point[,3:4])
  
  s.space <- as.list(c(0,0))
  
  s.space[[1]] <- ggplot() + looks +
    geom_polygon(data=hull.v, aes(x=A1, y=A2), color='grey', fill='white') +
    labs(x='PCoA1', y='PCoA2', title=titles[i]) +
    geom_polygon(data=s.hull.v, aes(x=A1, y=A2), fill=palette[i], alpha=0.25) +
    geom_point(data=s.point %>% left_join(., traits, by="Species"),
               aes(x=A1, y=A2, size=n), fill=palette[i], color='#666666', shape=21) +
    geom_text_repel(data=s.point %>% arrange(n) %>% top_n(5), aes(label=ShortSp, x=A1, y=A2), force=2,
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
    labs(x='PCoA3', y='PCoA4') +
    geom_polygon(data=s.hull.v2, aes(x=A3, y=A4), fill=palette[i], alpha=0.25) +
    geom_point(data=s.point %>% left_join(., traits, by="Species"), 
               aes(x=A3, y=A4, size=n), fill=palette[i], color='#666666', shape=21) +
    geom_text_repel(data=s.point %>% arrange(n) %>% top_n(5), aes(label=ShortSp, x=A3, y=A4), force=2,
                    size=2.8, fontface='italic', min.segment.length = 0, point.size=5,
                    max.overlaps = Inf, box.padding = 1, force_pull=0.1) +
    scale_size_continuous(range=c(1.5,8), limits=c(0,0.6), name="Relative abundance")  + 
    if (i != 4) {
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), 
            legend.position='none')
      } 
  else {theme()}
  
  t.space[[i]] <- s.space[[1]] / s.space[[2]] # stack them with patchwork and store
}

rm(s.point)
point$n <- Sp.count$n/624 # put in global abundances to visualise all sites
global1 <- ggplot() + looks +
  geom_polygon(data=hull.v, aes(x=A1, y=A2), fill='transparent', color='grey') +
  geom_point(data=left_join(point, traits, by="Species"), aes(x=A1, y=A2, shape=ForageMode), size=2, fill='black', alpha=0.5) +
  labs(x='PCoA1', y='PCoA2', title='All sites') +
  scale_fill_grey(aesthetics="color",start=0.2, end=0.7) + scale_shape_manual(values=c(1:6,21:25), name="Foraging mode")

global2 <- ggplot() + looks +
  geom_polygon(data=hull.v2, aes(x=A3, y=A4), fill='transparent', color='grey') +
  geom_point(data=left_join(point, traits, by="Species"), aes(x=A3, y=A4, shape=ForageMode), size=2, fill='black', alpha=0.5) +
  labs(x='PCoA3', y='PCoA4') + theme(legend.title=element_blank()) +
  scale_fill_grey(aesthetics="color", start=0.2, end=0.7) + scale_shape_manual(values=c(1:6,21:25), guide="none")

pred_long <- predictors %>% pivot_longer(cols=!Site, names_to="predictor", values_to="value")

FD_bar <- ggplot(pred_long %>% filter(predictor %in% c('FRic', 'FEve', 'FDiv'))) +
  geom_bar(aes(y=value, x=predictor, fill=Site), stat='identity', position='dodge', color='black') +
  labs(y='Index measure', x=NULL) +
  scale_fill_viridis_d(begin=0.1, end=0.95, guide="none", option="mako") + looks +
  scale_y_continuous(expand=expansion(mult=c(.0,.05)), limits=c(0,1)) +
  scale_x_discrete(labels=c('TDiv', 'TEve', 'TRic'))

A_bar <- ggplot(pred_long %>% filter(predictor == 'HerbProp' | predictor == 'BiterProp')) +
  geom_bar(aes(y=value, x=predictor, group=Site, fill=Site), stat='identity', position='dodge', color='black') +
  labs(y='Relative abundance', x=NULL) +
  scale_fill_viridis_d(begin=0.1, end=0.95, guide="none", option="mako") + looks + 
  scale_y_continuous(expand=expansion(mult=c(.0,.05))) +
  scale_x_discrete(labels=c('Benthic biters', 'Herbivores'))

bars <- (A_bar | FD_bar) * plot_layout(widths=c(2,3))

Fig2 <- (bars / ((global1/global2) | t.space[[1]] | t.space[[2]] | t.space[[3]]) / (t.space[[4]] | t.space[[5]] | t.space[[6]] | t.space[[7]])) * 
  plot_layout(guides='collect', heights = c(1,3,3)) * 
  theme(axis.text = element_text(size=9)) &
  theme(plot.title=element_text(face='bold', hjust=0.5, size=13))

rm(FD_bar, A_bar, bars, global1, global2, t.space, s.hull.v, s.hull.v2, s.space, hull.v, hull.v2)

ggsave(filename = "./figures/02_traitspaces.svg", device = "svg", width=30, height=30, units='cm', dpi=300)
                                                                                                                                                                         