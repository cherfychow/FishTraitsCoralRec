
#############################################################################

# FishTraitsCoralRec
# Supplementary figure: Global (all sites) trait distribution

#############################################################################

require(ggplot2)
require(patchwork)

looks <- theme_bw(base_size=13) + theme(panel.grid=element_blank(), axis.ticks=element_line(size=0.3))
source('analysis_code/function_convhull.R')
# trait space by trophic group
trait.points <- left_join(traits, point, by='Species')
Troph.2 <- as.list(rep(0,9))
Troph.4 <- as.list(rep(0,9))
for (i in c(1,2,4:9)) { # draw clusters around the different trophic groups
  Troph.2[[i]] <- convhull.vert(trait.points %>% filter(TrophParr == sort(unique(traits$TrophParr))[i]) %>% select(A1,A2))
  Troph.4[[i]] <- convhull.vert(trait.points %>% filter(TrophParr == sort(unique(traits$TrophParr))[i]) %>% select(A3,A4))
}
pal <- viridis::plasma(n=9, end=0.95, begin=0.1, alpha=0.3)

globalgroup1 <- ggplot() + looks +
  geom_point(data=trait.points, aes(x=A1, y=A2, fill=TrophParr), shape=21, size=3) +
  labs(x='PCo1', y='PCo2', title='Trophic group') + 
  theme(legend.position='none') +
  scale_fill_viridis_d(option='plasma', begin=0.1, end=0.95)
for (i in c(1,2,4:9)) {
  globalgroup1 <- globalgroup1 +
  geom_polygon(data=Troph.2[[i]], aes(x=A1, y=A2), fill=pal[i], color='#33333399', size=0.3)
}

globalgroup2 <- ggplot() + looks +
  geom_point(data=trait.points, aes(x=A3, y=A4, fill=TrophParr), shape=21, size=3) +
  labs(x='PCo3', y='PCo4') + 
  scale_fill_viridis_d(option='plasma', begin=0.1, end=0.95)
for (i in c(1,2,4:9)) {
  globalgroup2 <- globalgroup2 +
    geom_polygon(data=Troph.4[[i]], aes(x=A3, y=A4), fill=pal[i], color='#33333399', size=0.3)
}

# globalgroup1 / globalgroup2 + plot_layout(guides = "collect")

## ----global trait space by foraging---------------------------------------------------------------------------------------------

# hulls for everyone
hulls1 <- as.list(rep(0,4))
hulls2 <- as.list(rep(0,4))
trait.index <- c(2,4,5,6)
for (j in 2:4) {
  hulls1[[j]] <- data.frame(level=NA, A1=NA, A2=NA)
  hulls2[[j]] <- data.frame(level=NA, A3=NA, A4=NA)
  names(hulls1)[j] <- colnames(traits)[trait.index[j]]
  names(hulls2)[j] <- colnames(traits)[trait.index[j]]
  for (i in 1:length(unique(traits[,trait.index[j]]))) {
    hull1 <- convhull.vert(as_tibble(trait.points) %>% filter(.data[[colnames(trait.points)[trait.index[j]]]] == unique(as.vector(traits[,trait.index[j]]))[i]) %>% select(A1,A2))
    hull2 <- convhull.vert(as_tibble(trait.points) %>% filter(.data[[colnames(trait.points)[trait.index[j]]]] == unique(as.vector(traits[,trait.index[j]]))[i]) %>% select(A3,A4))
    hull1$level <- unique(traits[,trait.index[j]])[i]
    hull2$level <- unique(traits[,trait.index[j]])[i]
    hulls1[[j]] <- rbind(hulls1[[j]], hull1)
    hulls2[[j]] <- rbind(hulls2[[j]], hull2)
  }
}
for (i in c(1,4,6,8,9,10,11)) {
  hull1 <- convhull.vert(as_tibble(trait.points) %>% filter(ForageMode == sort(unique(traits$ForageMode))[i]) %>% select(A1,A2))
  hull2 <- convhull.vert(as_tibble(trait.points) %>% filter(ForageMode == sort(unique(traits$ForageMode))[i]) %>% select(A3,A4))
  hull1$level <- sort(unique(traits$ForageMode))[i]
  hull2$level <- sort(unique(traits$ForageMode))[i]
  hulls1[[1]] <- rbind(hulls1[[1]], hull1)
  hulls2[[1]] <- rbind(hulls2[[1]], hull2)
}

# remove the placeholder dummy rows
for (j in 1:4) {
  hulls1[[j]] <- hulls1[[j]][-1,]
  hulls2[[j]] <- hulls2[[j]][-1,]
}


# Forage Mode
global_forage1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_polygon(data=hulls1[[1]], aes(x=A1, y=A2, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=ForageMode), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_viridis_d() +
  theme(legend.position='none') +
  labs(title='Foraging Mode')
global_forage2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_polygon(data=hulls2[[1]], aes(x=A3, y=A4, fill=level), alpha=0.3, color='#55555599', size=0.3) +
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
  theme(legend.position = 'bottom') +
  scale_fill_gradient(low='black', high='white') +
  theme(legend.title=element_blank())

# Water column
global_col1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_polygon(data=hulls1[[2]], aes(x=A1, y=A2, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=ColumnFeed), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_manual(values = c('lightblue', 'navy', 'dodgerblue')) +
  theme(legend.position='none') +
  labs(title='Water Column Feeding')
global_col2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_polygon(data=hulls2[[2]], aes(x=A3, y=A4, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  scale_fill_manual(values = c('lightblue', 'navy', 'dodgerblue')) +
  geom_point(aes(fill=ColumnFeed), alpha=0.7, size=3, shape=21) + looks +
  theme(legend.position='bottom', legend.title=element_blank())

# Schooling
trait.points$Schooling <- as.factor(trait.points$Schooling)
hulls1[[3]]$level <- as.factor(hulls1[[3]]$level)
hulls2[[3]]$level <- as.factor(hulls2[[3]]$level)
global_s1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_polygon(data=hulls1[[3]], aes(x=A1, y=A2, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=Schooling), alpha=0.7, size=3, shape=21) + looks +
  theme(legend.position='none') + scale_fill_viridis_d() +
  labs(title='Schooling')
global_s2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_polygon(data=hulls2[[3]], aes(x=A3, y=A4, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=Schooling), alpha=0.7, size=3, shape=21) + looks + scale_fill_viridis_d() +
  theme(legend.position='bottom', legend.title=element_blank())

# Range
trait.points$Range <- as.factor(trait.points$Range)
hulls1[[4]]$level <- as.factor(hulls1[[4]]$level)
hulls2[[4]]$level <- as.factor(hulls2[[4]]$level)
global_r1 <- ggplot(data=trait.points, aes(x=A1, y=A2)) +
  geom_polygon(data=hulls1[[4]], aes(x=A1, y=A2, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=Range), alpha=0.7, size=3, shape=21) + looks +
  theme(legend.position='none') +
  labs(title='Range') + scale_fill_viridis_d(option='plasma')
global_r2 <- ggplot(data=trait.points, aes(x=A3, y=A4)) +
  geom_polygon(data=hulls2[[4]], aes(x=A3, y=A4, fill=level), alpha=0.3, color='#55555599', size=0.3) +
  geom_point(aes(fill=Range), alpha=0.7, size=3, shape=21) + looks +
  scale_fill_viridis_d(option='plasma') +
  theme(legend.position='bottom', legend.title=element_blank())

((globalgroup1 / globalgroup2) | (global_forage1 / global_forage2) | (global_troph1 / global_troph2))
((global_col1 / global_col2) | (global_r1 / global_r2) | (global_s1 / global_s2))
