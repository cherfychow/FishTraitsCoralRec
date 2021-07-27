
# Coral recruitment model partial regression plot
# Author: Cher Chow
# Date: 11 Aug 2020

require(tidyverse)
require(patchwork)
require(lme4)

set.seed(24)
ci.pred <- function(.) predict(., newx, type='response', re.form=NA)

m=8
l <- length(fixef(RecModels[[m]]))-1
par9 <- as.list(rep(0,l))
newfit.r9 <- as.list(rep(0,l))
conf.r9<- as.list(rep(0,l))
bb <- as.list(rep(0,l))
x <- RecruitData %>% ungroup() %>% mutate_if(is.numeric, mean)

# loop inputting the confidence intervals separately in case something goes wrong with bootstrapping.

for (i in 1:l) { # each predictor term
  newx <- x %>% ungroup() %>% select(!matches(names(fixef(RecModels[[m]]))[i+1])) %>% 
    bind_cols(y=RecruitData %>% select(matches(names(fixef(RecModels[[m]]))[i+1])))
  newfit.r9[[i]] <- predict(RecModels[[m]], newx, type='response', re.form=NA)
  # bootstrap confidence intervals
  bb[[i]] <- bootMer(RecModels[[m]], FUN=ci.pred, nsim=999)
  bb_se<-apply(bb[[i]]$t,2,function(x) x[order(x)][c(.05*999, .95*999)]) # dummy variable to get the 5% and 95%
  conf.r9[[i]] <- data.frame(lwr=bb_se[1,], upr=bb_se[2,])
}

axis.lab <- c('2018 Spat counts',
              'Foraging rate (cm-bites/min)',
              'Herbivore abundance')
names(axis.lab) <- names(fixef(RecModels[[m]])[-1])
sites <- c('Corner Beach', 'Lagoon', 'North Reef', 'Resort', 'Southeast', 'Turtle Beach', "Vicki's")

# for (i in 1:l) { # loop for partial regression plot panels
#   par9[[i]] = ggplot(bind_cols(RecruitData, conf.r9[[i]]) %>% bind_cols(., newy=newfit.r9[[i]]), aes_string(y='Recruits', x=names(fixef(RecModels[[m]]))[(i+1)])) +
#     geom_ribbon(aes(ymax=upr, ymin=lwr), fill='grey', alpha=0.4) +
#     geom_point(aes(fill=Site), size=3, alpha=0.8, shape=21, color='#777777') +
#     geom_line(aes(y=newy), color='black', size=0.5) +
#     labs(y='Recruit counts', x=axis.lab[names(fixef(RecModels[[m]]))[(i+1)]]) + looks +
#     scale_fill_viridis_d(begin=0.1, end=0.95) +
#     coord_cartesian(ylim=c(0,max(RecruitData$Recruits)+2)) +
#     if (i != l) {theme(legend.position='none')}
#   else {theme(legend.position='right')}
# }

for (i in 1:l) { # loop for partial regression plot panels
  par9[[i]] = ggplot(bind_cols(RecruitData, conf.r9[[i]]) %>% bind_cols(., newy=newfit.r9[[i]]), aes_string(y='Recruits', x=names(fixef(RecModels[[m]]))[(i+1)])) +
    geom_ribbon(aes(ymax=upr, ymin=lwr), fill='transparent', color='black', linetype='dashed', size=0.3) +
    geom_point(aes(fill=Site, shape=Site), size=3) +
    geom_line(aes(y=newy), color='black', size=0.5) +
    labs(y=ifelse(i==1, 'Recruit counts', ''), x=axis.lab[names(fixef(RecModels[[m]]))[(i+1)]]) + looks +
    scale_fill_viridis_d(begin=0.1, end=0.95, option='mako', name='Study site', labels=sites) +
    scale_shape_manual(values=c(3,4,21:25), name='Study site', labels=sites) +
    coord_cartesian(ylim=c(0,max(RecruitData$Recruits)+2)) +
    theme(legend.position=ifelse(i!=l,yes='none', no='right'), axis.text.y=if (i != 1) {element_blank()})
}

Fig5 <- (par9[[1]] + par9[[2]] + par9[[3]]) & plot_layout(guides='collect')
Fig5
ggsave(filename='./figures/05_recruitment.eps', device='eps', width=25, height=7.5, units='cm', dpi=300)
