
# Coral settlement model plotting
# Author: Cher Chow
# Date: 11 Aug 2020

require(tidyverse)
require(patchwork)
require(lme4)

set.seed(24)
ci.pred <- function(.) predict(., newx, type='response', re.form=NA)
looks <- theme_bw(base_size = 13) + theme(panel.grid = element_blank())

# Partial regression plot
l <- length(names(fixef(sett)))-1
par.reg <- as.list(rep(0,l))
newfit <- as.list(rep(0,l))
x <- SpatData %>% ungroup() %>% mutate_if(is.numeric, mean)
conf <- as.list(rep(0,l))
bb <- as.list(rep(0,l))
# loop inputting the confidence intervals separately in case something goes wrong with bootstrapping.

for (i in 1:l) { # each predictor term
  newx <- x %>% ungroup() %>% select(!matches(names(fixef(sett))[i+1])) %>% bind_cols(y=SpatData %>% select(matches(names(fixef(sett))[i+1])))
  newfit[[i]] <- predict(sett, newx, type='response', re.form=NA)
  # bootstrap confidence intervals
  bb[[i]] <- bootMer(sett, FUN=ci.pred, nsim=999)
  bb_se<-apply(bb[[i]]$t,2,function(x) x[order(x)][c(.05*999, .95*999)]) # dummy variable to get the 5% and 95%
  conf[[i]] <- data.frame(lwr=bb_se[1,], upr=bb_se[2,])
}

axis.lab <- c(
  'Herbivore abundance',
  'TDiv',
  'TOP',
  'Foraging rate (cm-bites/min)'
)
names(axis.lab) <- names(fixef(sett)[-1])
sites <- c('Corner Beach', 'Lagoon', 'North Reef', 'Resort', 'Southeast', 'Turtle Beach', "Vicki's")

for (i in 1:l) { # loop for partial regression plot panels
  par.reg[[i]] = ggplot(bind_cols(SpatData, conf[[i]]) %>% bind_cols(., newy=newfit[[i]]), aes_string(y='Spat', x=names(fixef(sett))[(i+1)])) +
    geom_ribbon(aes(ymax=upr, ymin=lwr), fill='#6dc7a9', color='transparent', alpha = 0.4) +
    geom_point(aes(fill=Site, shape=Site), size=2.5) +
    geom_line(aes(y=newy), color='#075c40', size=0.5) +
    labs(y=if (i == 4) {'Spat count'} else {NULL}, x=axis.lab[names(fixef(sett))[(i+1)]]) + looks +
    scale_fill_viridis_d(begin=0.1, end=0.95, option='mako', name='Study site', labels=sites) +
    scale_shape_manual(values=c(3,4,21:25), name='Study site', labels=sites) +
    coord_cartesian(ylim=c(0,max(SpatData$Spat)+2)) +
    theme(legend.position=ifelse(i!=l,yes='none', no='right'), axis.text.y=if (i != 4) {element_blank()})
}

Fig4 <- (par.reg[[4]] | par.reg[[1]] | par.reg[[2]] | par.reg[[3]]) + plot_layout(guides='collect')
Fig4
ggsave(plot = Fig4, filename=paste0(fig_dir, '/04_sett.svg'), device='svg', width=30, height=7.5, units='cm', dpi=300)

rm(par.reg, bb, bb_se, conf, newx, newfit)
