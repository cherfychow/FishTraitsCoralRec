
# FishTraitsCoralRec
# RUV species accumulation curves
# Assemblage species accumulation over time per site


# Setup and wrangling -----------------------------------------------------

require(dplyr)
require(readxl)
require(stringr)

ruv1 <- read.csv('src/UBRUV_species_2019.csv', header=T)
ruv2 <- read_excel('src/Video2_SE-NR-TB.xlsx')

sites <- c('North3' = 'NR_2', 'Southeast' = 'SE_2', 'TurtleBeach' = 'TB_2')
ruv2$Site <- str_replace_all(ruv2$Site, sites)

sum((colnames(ruv1) == colnames(ruv2)) + 0) # check that columns match
ruv.data <- bind_rows(ruv1, ruv2) # combine them
rm(ruv1, ruv2)

library(lubridate)
ruv.data$TimeElapsed <- ms(ruv.data$Timestamp)
ruv.data <- ruv.data %>% arrange(Site, Date, Camera, VidFile, TimeElapsed)
# sort by timestamp

#  timestamps are by video but we want them to be timestamps relative to total observation time, not video time.
cameras <- ruv.data %>% distinct(Site, Camera, VidFile)
ncam <- cameras %>% group_by(Site, Camera) %>% summarise(ncam = n_distinct(VidFile)) %>% pull(ncam)
seq <- ''
  for (i in 1:length(ncam)) {
    seq <- c(seq, 1:ncam[i])
  } # generate a sequence vector to identify the nominal order of video files
seq <- seq[-1] # trim that dummy start

cameras$VidSeq <- seq
rm(i, seq, ncam)

# add seq to the RUV observation data
ruv.data <- left_join(ruv.data, cameras, by=c('Site', 'Camera', 'VidFile'))
ruv.data$VidSeq <- as.numeric(ruv.data$VidSeq)

# use lubridate to parse timestamps as minute seconds
ruv.data$TimeElapsed[ which(ruv.data$VidSeq == 2 & ruv.data$Site != 'Vickis') ] <- ruv.data$TimeElapsed[ which(ruv.data$VidSeq == 2 & ruv.data$Site != 'Vickis')] + minutes(12)
# for the second video files, add 12 minutes (all except Vicki's)
ruv.data$TimeElapsed[ which(ruv.data$VidSeq == 3 & ruv.data$Site != 'Vickis') ] <- ruv.data$TimeElapsed[ which(ruv.data$VidSeq == 3 & ruv.data$Site != 'Vickis')] + minutes(24)
# and 24 minutes for the third videos
ruv.data$TimeElapsed[ which(ruv.data$VidSeq == 4 & ruv.data$Site != 'Vickis') ] <- ruv.data$TimeElapsed[ which(ruv.data$VidSeq == 4 & ruv.data$Site != 'Vickis')] + minutes(36)
# and 24 minutes for the third videos

# Vickis exception
for (i in 2:6) {
  ruv.data$TimeElapsed[which(ruv.data$VidSeq == i & ruv.data$Site == 'Vickis')] <- ruv.data$TimeElapsed[which(ruv.data$VidSeq == i & ruv.data$Site == 'Vickis')] + minutes( (i-1) * 5 )
}

# do a manual visual check
# some sites might have time elapsed go over 30 because observations had to start in the middle of the first video file etc

# have TimeElapsed reflect time stamps relative to first observation
ruv.data$TimeElapsed <- (seconds(ruv.data$TimeElapsed) %>% as.numeric)/60 # period objects to minutes

# create an object to index rows for each site
site.index <- as.list(rep('', 10))
for (i in 1:10) {
  site.index[[i]] <- which(ruv.data$Site == unique(ruv.data$Site)[i])
}

for (i in 1:10) {
  for (j in site.index[[i]]) {
    ruv.data$TElapsed[j] <- ruv.data$TimeElapsed[j] - ruv.data$TimeElapsed[site.index[[i]][1]]
  }
}

# Calculate SAC -----------------------------------------------------------

# make a dataframe to populate species accumulation per timestamp
SAC <- ruv.data %>% select(Site, TElapsed) # matches the row index of ruv.data

# for every row at every site, calculate the species number accumulation
for (i in 1:10) {
  for (j in site.index[[i]]) {
    if (j == site.index[[i]][1]) { # the first row from each site represents the first species record, so these will always start at 1
      SAC$NSpecies[j] <- 1
    }
    else {
      SAC$NSpecies[j] <- ruv.data$Species[site.index[[i]][1]:j] %>% n_distinct # number of unique species from the first row of the site to row j
    }
  }
}


# Resort has a big gap between minute 1 and 7, which will make model fitting difficult
# for times where species number stayed the same, fill them in
SAC <- SAC %>% bind_rows(., data.frame(Site = rep('Resort', 6), TElapsed = 2:7, NSpecies = rep(4, 6))) %>% arrange(Site, TElapsed)

# redo site.index
site.index <- as.list(rep('', 10))
for (i in 1:10) {
  site.index[[i]] <- which(SAC$Site == unique(SAC$Site)[i])
}

# check the SAC calculations by creating a species summary table
sp.count <- ruv.data %>% group_by(Site) %>% summarise(SpTotal=n_distinct(Species))
# does the final species count in SAC match the total we calculated?
for (i in 1:10) { print(SAC$NSpecies[site.index[[i]][length(site.index[[i]])]] == sp.count$SpTotal[i]) }

# model the SAC trend
require(nlstools)
set.seed(24)
SAC_all <- nls(formula = NSpecies ~ a - (a-b) * exp(-c * TElapsed), data=SAC, start = list(a=45, b=1, c=1))
SAC_pred <- nlsBootPredict(nlsBoot(SAC_all, niter=500), newdata = data.frame(TElapsed = SAC$TElapsed), interval="confidence")
colnames(SAC_pred)[2:3] <- c('lwr', 'upr')

# SAC curves for species accumulation for each site due to richness differences?

# model the overall trend again
SAC_s <- as.list(rep(0,10))
SACs_pred <- as.list(rep(0,10))
set.seed(24)
for (i in c(1:4,6:10)) {
  SAC_s[[i]] <- nls(formula = NSpecies ~ a - (a-b) * exp(-c* TElapsed), data=SAC[site.index[[i]],] %>% as.data.frame(), start = list(a=9, b=1, c=0.07))
  SACs_pred[[i]] <- nlsBootPredict(nlsBoot(SAC_s[[i]], niter=500), newdata = data.frame(TElapsed = SAC$TElapsed[site.index[[i]]]), interval="confidence") %>% as.data.frame
  colnames(SACs_pred[[i]])[2:3] <- c('lwr', 'upr')
}

# see if Gompertz function will fit Resort better
# b = slope around the inflection point
# c = lower asymptotte
# d = higher asymptote
# e = the time value producing a response halfway betewen d and c
SAC_s[[5]] <- nls(formula = NSpecies ~ c + (d-c) * exp( -exp(b * TElapsed - e)), data=SAC[site.index[[5]],] %>% as.data.frame(), start = list(b=0.2, c=1, d=18, e=6))
SACs_pred[[5]] <- nlsBootPredict(nlsBoot(SAC_s[[5]], niter=500), newdata = data.frame(TElapsed = SAC$TElapsed[site.index[[5]]]), interval="confidence") %>% as.data.frame
colnames(SACs_pred[[5]])[2:3] <- c('lwr', 'upr')

# Visualise ---------------------------------------------------------------

require(ggplot2)
require(stringr)

ggplot() +
  geom_ribbon(data=as.data.frame(SAC_pred) %>% bind_cols(., TElapsed=SAC$TElapsed), aes(ymin=lwr, ymax=upr, x=TElapsed), fill='grey80') +
  geom_line(data=SAC, aes(x=TElapsed, y=NSpecies, color=Site), size=0.5) +
  geom_line(data=data.frame(x=SAC$TElapsed, y=SAC_pred[,1]), aes(x, y), color='black', size=0.5) +
  theme_classic(base_size=13) +
  labs(x='Time elapsed (min)', y='Number of species observed') +
  theme(axis.line = element_line(size=0.2)) +
  scale_color_viridis_d(end=0.9, begin=0.1)


palette <- viridis::viridis(begin=0.2, end=0.9, n=10)
plot_SACs <- ggplot() +
  theme_classic(base_size=13) +
  geom_point(data=SAC %>% filter(!str_detect(Site, '\\_2')), aes(x=TElapsed, y=NSpecies, color=Site), size=0.5, shape=21, alpha=0.5) +
  labs(x='Time elapsed (min)', y='Cumulative number of species observed') +
  theme(axis.line = element_line(size=0.2)) +
  scale_color_viridis_d(end=0.9, begin=0.1) +
  coord_cartesian(ylim=c(0,45))

for (i in c(1:3,5,7,9,10)) {
 plot_SACs <- plot_SACs +
    geom_ribbon(data=SACs_pred[[i]] %>% bind_cols(., TElapsed=SAC$TElapsed[site.index[[i]]]), aes(ymin=lwr, ymax=upr, x=TElapsed), fill=palette[i], alpha=0.2) +
    geom_line(data=data.frame(x=SAC$TElapsed[site.index[[i]]], y=SACs_pred[[i]][,1]), aes(x, y), color=palette[i], size=0.5)
}
  
plot_SACs
