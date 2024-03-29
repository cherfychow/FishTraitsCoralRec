#############################################################################

# FishTraitsCoralRec
# Figure 3: Foraging rates by trophic group and functional group
# Author: Cher Chow

#############################################################################

require(tidyverse)
require(patchwork)
require(viridis)

looks <- theme_bw(base_size=13) + theme(panel.grid=element_blank(), axis.ticks=element_line(size=0.3))
source('https://gist.github.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/c9caceef463062c75c71b42607954f7958818ff7/cherulean.R')

# foraging rates by foraging mode
feed <- left_join(ForRate, sp.inf, by='Species') %>% mutate(SpFeed=0.01*SpInf*sumLBi) %>% 
  left_join(., traits %>% select(Species, TrophParr, ForageMode), by="Species")
feedF <- feed %>% group_by(Site,TrophParr, ForageMode) %>% summarise(Feed=sum(SpFeed), Infl=mean(SpInf))
feedF$ForageMode <- factor(feedF$ForageMode, levels=feedF %>% ungroup() %>% arrange(Infl) %>% distinct(ForageMode) %>% pull)

# foraging rates by trophic group
feedT <- feed %>% group_by(Site,TrophParr) %>% summarise(Feed=sum(SpFeed), Infl=mean(SpInf))
feedT$TrophParr <- factor(feedT$TrophParr, levels=feedT %>% ungroup() %>% arrange(Infl) %>% distinct(TrophParr) %>% pull)

feedforage <- ggplot(feedF, aes(x=Site, fill=Feed, y=ForageMode)) +
  geom_tile() + looks + 
  scale_fill_gradient(low='#EEEEEE', high='black', name=NULL) + 
  scale_x_discrete(labels=c('CB', 'L', 'NR', 'R', 'SE', 'TB', 'V')) + 
  labs(x='Site', y='Foraging mode', title="b")
feedtroph <- ggplot(feedT, aes(x=Site, fill=Feed, y=TrophParr)) +
  geom_tile() + looks + 
  scale_fill_gradient(low='#EEEEEE', high='black', name='Foraging rate', limits=c(0,50)) + 
  scale_x_discrete(labels=c('CB', 'L', 'NR', 'R', 'SE', 'TB', 'V')) + 
  labs(x='Site', y='Trophic group', title="a")
feed_total <- ggplot(feed, aes(x=Site, y=SpFeed)) + 
  geom_boxplot(width=0.6, aes(fill = Site)) + looks + 
  geom_jitter(shape = 21, color = "#222222", fill = 'transparent', size = 1.5, width = 0.1) +
  scale_x_discrete(labels=c('CB', 'L', 'NR', 'R', 'SE', 'TB', 'V')) + 
  scale_fill_cherulean(alpha = 0.6, palette = "globiceps", discrete = T, guide="none") +
  labs(x='Study site', y='Foraging rate (cm-bites/min)', title="c")

Fig3 <- (feedtroph + feedforage + feed_total) & theme(legend.position='bottom')
# (feedtroph + feedforage) + plot_layout(guides='collect')

Fig3
ggsave(filename = paste0(fig_dir, '03_foragingrates.svg'), device='svg', width=27, height=10, units='cm', dpi=300)
rm(feedtroph, feedforage, feed_total)
