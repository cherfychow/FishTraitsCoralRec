
#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Supplementary figure: Global trait distribution
# Author: Cher Chow

#############################################################################

# trait space by trophic group
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