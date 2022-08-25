
# Coral settlement model plotting
# Author: Cher Chow
# Date: 11 Aug 2020

require(tidyverse)
require(patchwork)
require(lme4)

set.seed(24)
ci.pred <- function(.) predict(., newx, type='response', re.form=NA)
looks <- theme_classic(base_size = 13) + theme(panel.grid = element_blank())

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
  'Benthic forager abundance',
  'Foraging rate (cm-bites/min)',
  'Herbivore abundance',
  'TDiv',
  'TEve',
  'TOP'
)
names(axis.lab) <- names(fixef(sett)[-1])
sites <- c('Corner Beach', 'Lagoon', 'North Reef', 'Resort', 'Southeast', 'Turtle Beach', "Vicki's")

for (i in 1:l) { # loop for partial regression plot panels
  par.reg[[i]] = ggplot(bind_cols(SpatData, conf[[i]]) %>% bind_cols(., newy=newfit[[i]]), aes_string(y='Spat', x=names(fixef(sett))[(i+1)])) +
    geom_ribbon(aes(ymax=upr, ymin=lwr), fill='transparent', color='black', linetype='dashed', size=0.3) +
    geom_point(aes(fill=Site, shape=Site), size=3) +
    geom_line(aes(y=newy), color='black', size=0.5) +
    labs(y=NULL, x=axis.lab[names(fixef(sett))[(i+1)]]) + looks +
    scale_fill_viridis_d(begin=0.1, end=0.95, option='mako', name='Study site', labels=sites) +
    scale_shape_manual(values=c(3,4,21:25), name='Study site', labels=sites) +
    coord_cartesian(ylim=c(0,max(SpatData$Spat)+2)) +
    theme(legend.position=ifelse(i!=l,yes='none', no='right'), axis.text.y=if (i %in% c(2,4)) {element_blank()})
}

Fig4 <- par.reg[[1]] + par.reg[[2]] + par.reg[[3]] + par.reg[[4]] + par.reg[[5]] + par.reg[[6]] + plot_layout(guides='collect')
Fig4
ggsave(filename='./figures/04_sett_benth.eps', device='eps', width=17, height=14, units='cm', dpi=300)

rm(par.reg, bb, bb_se, conf, newx, newfit)
