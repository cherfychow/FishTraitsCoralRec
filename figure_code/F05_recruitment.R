
# Coral recruitment model partial regression plot
# Author: Cher Chow
# Date: 11 Aug 2020

require(tidyverse)
require(patchwork)
require(lme4)

set.seed(24)
ci.pred <- function(.) predict(., newx, type='response', re.form=NA)

l <- length(fixef(recr))-1
recr_par <- as.list(rep(0,l))
recr_newfit <- as.list(rep(0,l))
recr_conf<- as.list(rep(0,l))
bb <- as.list(rep(0,l))
x <- RecruitData %>% ungroup() %>% mutate_if(is.numeric, mean)

# loop inputting the confidence intervals separately in case something goes wrong with bootstrapping.

for (i in 1:l) { # each predictor term
  newx <- x %>% ungroup() %>% select(!matches(names(fixef(recr))[i+1])) %>% 
    bind_cols(y=RecruitData %>% ungroup() %>% select(matches(names(fixef(recr))[i+1]))) # everything at its mean except for the predictor term i
  recr_newfit[[i]] <- predict(recr, newx, type='response', re.form=NA)
  # bootstrap confidence intervals
  bb[[i]] <- bootMer(recr, FUN=ci.pred, nsim=999)
  bb_se<-apply(bb[[i]]$t,2,function(x) x[order(x)][c(.05*999, .95*999)]) # dummy variable to get the 5% and 95%
  recr_conf[[i]] <- data.frame(lwr=bb_se[1,], upr=bb_se[2,])
}

axis.lab <- c('Herbivore abundance',
              'TOP',
              'Foraging rate (cm-bites/min)',
              '2018 Spat counts')
names(axis.lab) <- names(fixef(recr)[-1])
sites <- c('Corner Beach', 'Lagoon', 'North Reef', 'Resort', 'Southeast', 'Turtle Beach', "Vicki's")

for (i in 1:l) { # loop for partial regression plot panels
  recr_par[[i]] <- ggplot(bind_cols(RecruitData, recr_conf[[i]]) %>% bind_cols(., newy=recr_newfit[[i]]), aes_string(y='Recruits', x=names(fixef(recr))[(i+1)])) +
    geom_ribbon(aes(ymax=upr, ymin=lwr), fill='#a7a9eb', color='transparent', alpha = 0.4) +
    geom_point(aes(fill=Site, shape=Site), size=2.5) +
    geom_line(aes(y=newy), color='#242780', size=0.5) +
    labs(y=if (i == 4) {'Recruit count'} else {NULL}, x=axis.lab[names(fixef(recr))[(i+1)]]) + looks +
    scale_fill_viridis_d(begin=0.1, end=0.95, option='mako', name='Study site', labels=sites) +
    scale_shape_manual(values=c(3,4,21:25), name='Study site', labels=sites) +
    coord_cartesian(ylim=c(0,max(RecruitData$Recruits)+2)) +
    theme(legend.position=ifelse(i!=l,yes='none', no='right'), axis.text.y=if (i != 4) {element_blank()})
}

Fig5 <- (recr_par[[4]] + recr_par[[3]] + recr_par[[1]] + recr_par[[2]]) + plot_layout(guides='collect', nrow = 1)
Fig5
ggsave(plot = Fig5, filename=paste0(fig_dir, '05_recruitment.svg'), device='svg', width=30, height=7.5, units='cm', dpi=300)
