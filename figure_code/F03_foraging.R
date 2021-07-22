#############################################################################

# FishTraitsxCoralRec
# Responses to fish trait diversity in coral settlement and recruitment
# Figure 3: Foraging rates by trophic group and functional group
# Author: Cher Chow

#############################################################################

require(tidyverse)
require(patchwork)
require(viridis)

feedforage <- ggplot(feedF, aes(x=Site, fill=Feed, y=ForageMode)) +
  geom_tile() + looks + 
  scale_fill_gradient(low='#EEEEEE', high='black', name=NULL) + 
  scale_x_discrete(labels=c('CB', 'L', 'NR', 'R', 'SE', 'TB', 'V')) + 
  labs(x='Site', y='Foraging mode', title="b")
feedtroph <- ggplot(feedT, aes(x=Site, fill=Feed, y=TrophicGroup)) +
  geom_tile() + looks + 
  scale_fill_gradient(low='#EEEEEE', high='black', name='Foraging rate', limits=c(0,50)) + 
  scale_x_discrete(labels=c('CB', 'L', 'NR', 'R', 'SE', 'TB', 'V')) + 
  labs(x='Site', y='Trophic group', title="a")
feed_total <- ggplot(feed, aes(x=Site, y=SpFeed, fill=Site)) + 
  geom_boxplot(width=0.75) + looks + 
  scale_x_discrete(labels=c('CB', 'L', 'NR', 'R', 'SE', 'TB', 'V')) + 
  scale_fill_viridis(alpha=0.6, option="mako", end=0.95, begin=0.1, discrete = T, guide="none") +
  labs(x='Study site', y='Foraging rate (cm-bites/min)', title="c")

Fig3 <- (feedtroph + feedforage + feed_total) & theme(legend.position='bottom')
# (feedtroph + feedforage) + plot_layout(guides='collect')

Fig3
rm(feedtroph, feedforage, feed_total)
ggsave(filename = './figures/03_foragingrates.svg', device='svg', width=27, height=10, units='cm', dpi=300)
