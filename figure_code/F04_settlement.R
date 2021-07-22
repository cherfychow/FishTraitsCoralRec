
# Coral settlement model plotting
# Author: Cher Chow
# Date: 11 Aug 2020

require(tidyverse)
require(patchwork)
require(lme4)

set.seed(24)
ci.pred <- function(.) predict(., newx, type='response', re.form=NA)

# Partial regression plot
m=4
l <- length(names(fixef(SpatModel.nb[[m]])))-1
par.reg4 <- as.list(rep(0,l))
newfit4 <- as.list(rep(0,l))
x <- SpatData %>% ungroup() %>% mutate_if(is.numeric, mean)
conf4 <- as.list(rep(0,l))
bb <- as.list(rep(0,l))
# loop inputting the confidence intervals separately in case something goes wrong with bootstrapping.

for (i in 1:l) { # each predictor term
  newx <- x %>% ungroup() %>% select(!matches(names(fixef(SpatModel.nb[[m]]))[i+1])) %>% bind_cols(y=SpatData %>% select(matches(names(fixef(SpatModel.nb[[m]]))[i+1])))
  newfit4[[i]] <- predict(SpatModel.nb[[m]], newx, type='response', re.form=NA)
  # bootstrap confidence intervals
  bb[[i]] <- bootMer(SpatModel.nb[[m]], FUN=ci.pred, nsim=999)
  bb_se<-apply(bb[[i]]$t,2,function(x) x[order(x)][c(.05*999, .95*999)]) # dummy variable to get the 5% and 95%
  conf4[[i]] <- data.frame(lwr=bb_se[1,], upr=bb_se[2,])
}

axis.lab <- c(
  'Trait evenness',
  'Trait richness',
  'Trait divergence',
  'Foraging rate (cm-bites/min)',
  'Herbivore abundance',
  'Benthic forager abundance'
)
names(axis.lab) <- names(fixef(SpatModel.nb[[m]])[-1])
sites <- c('Corner Beach', 'Lagoon', 'North Reef', 'Resort', 'Southeast', 'Turtle Beach', "Vicki's")

for (i in 1:l) { # loop for partial regression plot panels
  par.reg4[[i]] = ggplot(bind_cols(SpatData, conf4[[i]]) %>% bind_cols(., newy=newfit4[[i]]), aes_string(y='Spat', x=names(fixef(SpatModel.nb[[m]]))[(i+1)])) +
    geom_ribbon(aes(ymax=upr, ymin=lwr), fill='transparent', color='black', linetype='dashed', size=0.3) +
    geom_point(aes(fill=Site, shape=Site), size=3) +
    geom_line(aes(y=newy), color='black', size=0.5) +
    labs(y=NULL, x=axis.lab[names(fixef(SpatModel.nb[[m]]))[(i+1)]]) + looks +
    scale_fill_viridis_d(begin=0.1, end=0.95, option='mako', name='Study site', labels=sites) +
    scale_shape_manual(values=c(3,4,21:25), name='Study site', labels=sites) +
    coord_cartesian(ylim=c(0,max(SpatData$Spat)+2)) +
    theme(legend.position=ifelse(i!=l,yes='none', no='right'), axis.text.y=if (i %in% c(2,4)) {element_blank()})
}

Fig4 <- par.reg4[[1]] + par.reg4[[2]] + par.reg4[[3]] + par.reg4[[4]] + plot_layout(guides='collect')
Fig4
ggsave(filename='./figures/04_settlement.eps', device='eps', width=17, height=14, units='cm', dpi=300)

rm(par.reg4, bb, bb_se, conf4, newx, newfit4)