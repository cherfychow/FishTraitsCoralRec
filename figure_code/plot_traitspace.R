
#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Trait space figures (by site, and also globally showing by each trait)
# Author: Cher Chow
# Updated: 26 Mar 2021

#############################################################################

t.space <- as.list(rep(0,4))
point <- as.data.frame(FD_global$x.axes)[,1:4]
point$Species <- rownames(traits)

source('convhullvert_function.R') # wrote my own function to make plottable convex hull vertices
hull.v <- convhull.vert(point[,1:2])
hull.v2 <- convhull.vert(point[,3:4])

# set up the base layers of the global functional trait space
global1 <- ggplot() + looks +
  geom_polygon(data=hull.v, aes(x=A1, y=A2), color='grey', fill='white') +
  labs(x='Axis 1', y='Axis 2')
global2 <- ggplot() + looks +
  geom_polygon(data=hull.v2, aes(x=A3, y=A4), color='grey', fill='white') +
  labs(x='Axis 3', y='Axis 4')

# palette for sites
palette <- viridis::viridis(n=7, end=0.95, begin=0.1, alpha=0.7)
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
  
  s.space[[1]] <- global1 +
    geom_polygon(data=s.hull.v, aes(x=A1, y=A2), fill=palette[i], alpha=0.25) +
    geom_point(data=s.point, aes(x=A1, y=A2, size=n), fill=palette[i], color='#666666', shape=21) +
    # geom_text_repel(data=s.point, aes(label=ShortSp, x=A1, y=A2), force=2,
    #                 size=3.6, fontface='italic', min.segment.length = 0, point.size=5,
    #                 max.overlaps = Inf, box.padding = 1, force_pull=0.1) +
    scale_size_continuous(range=c(1.5,8), guide=F, limits=c(0,0.6)) +
    labs(title=unique(UBRUVspecies$Site)[i]) +
    if (i != 4) {theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position='none')} else {theme()}
  s.space[[2]] <- global2 +
    geom_polygon(data=s.hull.v2, aes(x=A3, y=A4), fill=palette[i], alpha=0.25) +
    geom_point(data=s.point, aes(x=A3, y=A4, size=n), fill=palette[i], color='#666666', shape=21) +
    # geom_text_repel(data=s.point, aes(label=ShortSp, x=A3, y=A4), force=2,
    #                 size=3.6, fontface='italic', min.segment.length = 0, point.size=5,
    #                 max.overlaps = Inf, box.padding = 1, force_pull=0.1) +
    scale_size_continuous(range=c(1.5,8), limits=c(0,0.6))  + 
    if (i != 4) {theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position='none')} else {theme()}
  
  t.space[[i]] <- s.space[[1]] / s.space[[2]] # stack them with patchwork and store
}

rm(s.point)
point$n <- Sp.count$n/624 # put in global abundances to visualise all sites
global1 <- ggplot() + looks +
  geom_polygon(data=hull.v, aes(x=A1, y=A2), fill='#EEEEEE') +
  geom_point(data=point, aes(x=A1, y=A2, size=n), alpha=0.3, shape=21) +
  labs(x='Axis 1', y='Axis 2', title='All sites') +
  theme(plot.title=element_text(face='bold', hjust=0.5, size=13)) +
  scale_size_continuous(range=c(1.5,8), limits=c(0,0.5), guide=F)

global2 <- ggplot() + looks +
  geom_polygon(data=hull.v2, aes(x=A3, y=A4), fill='#EEEEEE') +
  geom_point(data=point, aes(x=A3, y=A4, size=n), alpha=0.3, shape=21) +
  labs(x='Axis 3', y='Axis 4') +
  scale_size_continuous(range=c(1.5,8), limits=c(0,0.5)) + theme(legend.title=element_blank())

FD_bar <- ggplot(predictors %>% select(Site, FRic, FEve, FDiv) %>% pivot_longer(c(FRic, FEve, FDiv))) +
  geom_bar(aes(y=value, x=Site, alpha=name, fill=Site), stat='identity', position='dodge', color='black') +
  labs(y='Index measure', x='Study sites') +
  scale_fill_viridis_d(begin=0.1, end=0.95, guide=F) + looks + scale_alpha_discrete(name=NULL, range=c(0.2, 1)) +
  scale_y_continuous(expand=expansion(mult=c(.0,.05)), limits=c(0,1)) +
  scale_x_discrete(labels=c('CB', 'L1', 'N3', 'R', 'SE', 'TB', 'V'))

A_bar <- ggplot(predictors %>% select(Site, HerbProp, BiterProp) %>% pivot_longer(c(HerbProp, BiterProp))) +
  geom_bar(aes(y=value, x=Site, alpha=name, fill=Site), stat='identity', position='dodge', color='black') +
  labs(y='Relative abundance', x='Study sites') +
  scale_fill_viridis_d(begin=0.1, end=0.95, guide=F) + looks + scale_alpha_discrete(name=NULL, range=c(0.5, 1), labels=c('Benthic biters', 'Herbivores')) +
  scale_y_continuous(expand=expansion(mult=c(.0,.05))) +
  scale_x_discrete(labels=c('CB', 'L1', 'N3', 'R', 'SE', 'TB', 'V'))

FD_bar1 <- ggplot(data=predictors, aes(x=Site)) +
  geom_bar(aes(y=FRic, fill=Site), color='black', stat='identity') + looks +
  scale_fill_viridis_d(begin=0.1, end=0.95, alpha=0.7) +
  labs(x=NULL, y=NULL, title='FRic') +
  scale_y_continuous(expand=expansion(mult=c(.0,.05))) +
  theme(legend.position='none', axis.text.x=element_blank(), axis.ticks.x=element_blank())
FD_bar2 <- ggplot(data=predictors, aes(x=Site)) +
  geom_bar(aes(y=FEve, fill=Site), color='black', stat='identity') + looks +
  scale_fill_viridis_d(begin=0.1, end=0.95, alpha=0.7) +
  labs(x=NULL, y='Index measure', title='FEve') +
  scale_y_continuous(expand=expansion(mult=c(.0,.05))) +
  theme(legend.position='none', axis.text.x=element_blank(), axis.ticks.x=element_blank())
FD_bar3 <- ggplot(data=predictors, aes(x=Site)) +
  geom_bar(aes(y=FDiv, fill=Site), color='black', stat='identity') + looks +
  scale_fill_viridis_d(begin=0.1, end=0.95, guide=F, alpha=0.7) +
  labs(x='Study sites', y=NULL, title='FDiv') +
  scale_y_continuous(expand=expansion(mult=c(.0,.05))) +
  scale_x_discrete(labels=c('CB', 'L1', 'N3', 'R', 'SE', 'TB', 'V'))

((global1 / global2) | t.space[[1]]) + plot_layout(guides='collect')
(t.space[[2]] | t.space[[3]]) + plot_layout(guides='collect')
(t.space[[4]] | t.space[[5]]) + plot_layout(guides='collect')
(t.space[[6]] | t.space[[7]]) + plot_layout(guides='collect')
(FD_bar1 + FD_bar2 + FD_bar3) + plot_layout(guides='collect') & theme(axis.text.x=element_blank())

( (global1/global2) | t.space[[1]] | t.space[[2]] | t.space[[3]]) + plot_layout(guides='collect')
(t.space[[4]] | t.space[[5]] | t.space[[6]] | t.space[[7]] ) + plot_layout(guides='collect')
(t.space[[7]] | (global1/global2))  + plot_layout(guides='collect')

## ----Global trait space by trophic group----------------------------------------------------------------------------------------
trait.points <- left_join(traits, point, by='Species')
Troph.2 <- as.list(rep(0,6))
Troph.4 <- as.list(rep(0,6))
for (i in c(2,4,5,6,7,8)) { # draw clusters around the different trophic groups
  Troph.2[[i]] <- convhull.vert(trait.points %>% filter(TrophicGroup == sort(unique(traits$TrophicGroup))[i]) %>% select(A1,A2))
  Troph.4[[i]] <- convhull.vert(trait.points %>% filter(TrophicGroup == sort(unique(traits$TrophicGroup))[i]) %>% select(A3,A4))
}
pal <- viridis::plasma(n=8, end=0.95, begin=0.1, alpha=0.3)
palc <- viridis::plasma(n=8, end=0.95, begin=0.1)
globalgroup1 <- ggplot() + looks +
  geom_polygon(data=Troph.2[[2]], aes(x=A1, y=A2), fill=pal[2], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.2[[4]], aes(x=A1, y=A2), fill=pal[4], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.2[[5]], aes(x=A1, y=A2), fill=pal[5], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.2[[6]], aes(x=A1, y=A2), fill=pal[6], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.2[[7]], aes(x=A1, y=A2), fill=pal[7], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.2[[8]], aes(x=A1, y=A2), fill=pal[8], color='#55555599', size=0.3) +
  geom_point(data=trait.points, aes(x=A1, y=A2, fill=TrophicGroup), shape=21, size=3) +
  labs(x='Axis 1', y='Axis 2', title='Trophic group') + 
  theme(legend.position='none') +
  scale_fill_viridis_d(option='plasma', begin=0.1, end=0.95)
globalgroup2 <- ggplot() + looks +
  geom_polygon(data=Troph.4[[2]], aes(x=A3, y=A4), fill=pal[2], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.4[[4]], aes(x=A3, y=A4), fill=pal[4], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.4[[5]], aes(x=A3, y=A4), fill=pal[5], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.4[[6]], aes(x=A3, y=A4), fill=pal[6], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.4[[7]], aes(x=A3, y=A4), fill=pal[7], color='#55555599', size=0.3) +
  geom_polygon(data=Troph.4[[8]], aes(x=A3, y=A4), fill=pal[8], color='#55555599', size=0.3) +
  geom_point(data=trait.points, aes(x=A3, y=A4, fill=TrophicGroup), shape=21, size=3) +
  labs(x='Axis 3', y='Axis 4') + 
  scale_fill_viridis_d(option='plasma', begin=0.1, end=0.95) +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  theme(legend.title=element_blank())

globalgroup1 / globalgroup2 + plot_layout(guides='collect')


## ----global trait space by foraging---------------------------------------------------------------------------------------------
# Trophic Group
global_group1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_point(aes(fill=TrophicGroup), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_viridis_d(option='plasma', begin=0.1, end=0.95) +
  theme(legend.position='none') +
  labs(title='Trophic Group')
global_group2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_point(aes(fill=TrophicGroup), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_viridis_d(option='plasma', begin=0.1, end=0.95) +
  theme(legend.position='bottom', legend.title=element_blank())

# hulls for everyone
hulls1 <- as.list(rep(0,4))
hulls2 <- as.list(rep(0,4))
for (j in c(5,6,7)) {
  hulls1[[j]] <- data.frame(level=NA, A1=NA, A2=NA)
  hulls2[[j]] <- data.frame(level=NA, A3=NA, A4=NA)
  names(hulls1)[j] <- colnames(traits)[j]
  names(hulls2)[j] <- colnames(traits)[j]
  for (i in 1:length(unique(traits[,j]))) {
    hull1 <- convhull.vert(as_tibble(trait.points) %>% filter(.data[[colnames(trait.points)[j]]] == unique(as.vector(traits[,j]))[i]) %>% select(A1,A2))
    hull2 <- convhull.vert(as_tibble(trait.points) %>% filter(.data[[colnames(trait.points)[j]]] == unique(as.vector(traits[,j]))[i]) %>% select(A3,A4))
    hull1$level <- unique(traits[,j])[i]
    hull2$level <- unique(traits[,j])[i]
    hulls1[[j]] <- rbind(hulls1[[j]], hull1)
    hulls2[[j]] <- rbind(hulls2[[j]], hull2)
  }
}
for (i in c(1,4,6,8,9,10,11)) {
  hull1 <- convhull.vert(as_tibble(trait.points) %>% filter(ForageMode == sort(unique(traits$ForageMode))[i]) %>% select(A1,A2))
  hull2 <- convhull.vert(as_tibble(trait.points) %>% filter(ForageMode == sort(unique(traits$ForageMode))[i]) %>% select(A3,A4))
  hull1$level <- sort(unique(traits$ForageMode))[i]
  hull2$level <- sort(unique(traits$ForageMode))[i]
  hulls1[[4]] <- rbind(hulls1[[4]], hull1)
  hulls2[[4]] <- rbind(hulls2[[4]], hull2)
}


# Forage Mode
global_forage1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_polygon(data=hulls1[[4]], aes(x=A1, y=A2, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=ForageMode), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_viridis_d() +
  theme(legend.position='none') +
  labs(title='Foraging Mode')
global_forage2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_polygon(data=hulls2[[4]], aes(x=A3, y=A4, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=ForageMode), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_viridis_d() +
  theme(legend.title=element_blank())

# Trophic level
global_troph1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_point(aes(fill=FoodTroph), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_gradient(low='black', high='white') +
  theme(legend.position='none') +
  labs(title='Trophic Level')
global_troph2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_point(aes(fill=FoodTroph), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_gradient(low='black', high='white') +
  theme(legend.title=element_blank())

# Water column
global_col1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_polygon(data=hulls1[[5]], aes(x=A1, y=A2, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=ColumnFeed), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_grey() +
  theme(legend.position='none') +
  labs(title='Water Column Feeding')
global_col2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_polygon(data=hulls2[[5]], aes(x=A3, y=A4, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  scale_fill_grey() +
  geom_point(aes(fill=ColumnFeed), alpha=0.7, size=3, shape=21) + looks +
  theme(legend.position='bottom', legend.title=element_blank())

# Range
global_r1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_polygon(data=hulls1[[7]], aes(x=A1, y=A2, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=Range), alpha=0.7, size=3, shape=21) + looks +
  theme(legend.position='none') +
  labs(title='Range') + scale_fill_viridis_d(option='plasma')
global_r2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_polygon(data=hulls2[[7]], aes(x=A3, y=A4, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=Range), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_viridis_d(option='plasma') +
  theme(legend.position='bottom', legend.title=element_blank())

# Schooling
global_s1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_polygon(data=hulls1[[6]], aes(x=A1, y=A2, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=Schooling), alpha=0.7, size=3, shape=21) + looks +
  theme(legend.position='none') + scale_fill_viridis_d() +
  labs(title='Schooling')
global_s2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_polygon(data=hulls2[[6]], aes(x=A3, y=A4, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=Schooling), alpha=0.7, size=3, shape=21) + looks + scale_fill_viridis_d() +
  theme(legend.position='bottom', legend.title=element_blank())

((globalgroup1 / globalgroup2) | (global_forage1 / global_forage2) | (global_troph1 / global_troph2)) + plot_layout(guides='collect')
((global_col1 / global_col2) | (global_r1 / global_r2) | (global_s1 / global_s2))