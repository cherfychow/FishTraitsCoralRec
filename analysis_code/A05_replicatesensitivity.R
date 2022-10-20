
# FishTraitsCoralRec
# RUV sensitivity analysis
# Assemblage sampling sensitivity


# Setup and data clean--------------------------------------------------------

require(dplyr) # data handling
require(readxl) # read xls files
require(worms) # package for validating species/taxa names
require(stringr)

set.seed(24) # repeatability :)

# needed functions
source('analysis_code/TOP_Fontanaetal2015.R')
source('./analysis_code/function_convhull.R') # make plottable convex hull vertices, plotting purposes

# read in assemblage data from folder
ruv2.data <- read_excel('src/Video2_SE-NR-TB.xlsx')
head(ruv2.data) # check read in
detach(package:readxl)

## Check and trim species data------------------------------------------------
str(ruv2.data)
# fix column data types first
ruv2.data$Date <- as.POSIXct(ruv2.data$Date, tryFormats='%Y-%m-%d')
ruv2.data$Site <- as.factor(ruv2.data$Site)

summary(ruv2.data)
length(unique(ruv2.data$Species)) # 72 species recorded

ruv2.data %>% filter(str_detect(Species, 'idae$| sp.$| sp$')) %>% nrow # how many uncertain IDs
# if acceptable number, just trim the data to exclude uncertain IDs
ruv2.data <- ruv2.data %>% filter(!str_detect(Species, 'idae$| sp.$| sp$'))
# 215 observations left, removed 1 row1

## Taxon checks and clean------------------------------------------------------

ruv2.data$Species %>% unique %>% sort # visual check
ruv2.data$Species <- str_replace_all(ruv2.data$Species, 'annulatus', 'brevirostris') # species ID correction
ruv2.data <- ruv2.data %>% filter(!str_detect(Species, 'juvenile')) # remove the unknown juvenile IDs
Sp <- sort(unique(ruv2.data$Species))

## Warning: WoRMS is computationally heavy! Might take a few minutes
Sp.list <- wormsbymatchnames(Sp, marine_only=T) # retrieve WoRMS checks
Sp.list$ID <- seq(1, length(Sp.list$valid_name), by=1) # create an ID column

# create a list of invalid species names to check
invalid <- Sp.list %>% filter(match_type != 'exact' | status != 'accepted') %>% select(ID, valid_name, family, genus)
invalid$before <- Sp[invalid$ID] # retrieve the invalid species names

validate <- invalid$valid_name
names(validate) <- invalid$before # make a spell check object
ruv2.data$Species <- str_replace_all(ruv2.data$Species, validate) # replace misspellings
rm(invalid)
detach(package:worms)


## Clean species lists and tables----------------------------------------------

require(dplyr)
require(stringr)

Sp.list <- Sp.list %>% select(valid_name, family) %>% mutate(genus=word(valid_name, 1))
# let's trim the worms table to just taxonomic info, and extract genus from the valid name

# at this point, restarting the R session might be needed. Unknown package interactions.
Sp.count <- ruv2.data %>% 
  dplyr::group_by(Species) %>% dplyr::summarise(n=sum(Count))  # counts per species

site.n <- ruv2.data %>% group_by(Site) %>% dplyr::summarise(n=sum(Count)) # total number of fish per site

fish.ruv <- ruv2.data %>% group_by(Site, Species) %>% dplyr::summarize(n=sum(Count)) # aggregate species counts

colnames(Sp.list)[1] <- "Species" #rename the Species column in the worms dataframe

Sp.count$genus <- word(Sp.count$Species, 1) # add genus column to the species count summary tibble
Sp.count <- left_join(Sp.count, distinct(Sp.list[,2:3]), by='genus') # and add family taxonomic info too


# Trait space construction, new only ------------------------------------------------

require(FD)
require(tidyverse)
traits2 <- read.csv('src/fish_traits2.csv', header=T) # read in the traits table for the second set of videos
traits <- read.csv('src/fish_traits.csv', header=T)

# create the abundance matrix for functional evenness and divergence weighting
# species are columns, each site is a row
# long format to wide
abund <- Sp.count %>% filter(!str_detect(Species, ' sp$')) %>% 
  select(Species, n) %>% 
  pivot_wider(names_from = Species, values_from = n) # row 1 = all sites
abund$Site <- 'All'
abund <- rbind(fish.ruv %>% filter(!str_detect(Species, ' sp$')) %>% select(Site, Species, n) %>% spread(Species, n), abund)
# because not all species occur in each site, there are lots of NAs to be filled with 0s
NAlist <- as.list(rep(0,ncol(abund)-1)) # list has to have column names for every item = 0
names(NAlist) <- colnames(abund)[2:ncol(abund)]
abund <- abund %>% replace_na(NAlist) # annnd replace NAs

# abundance data frame is in characters????
abund <- as.data.frame(abund)
rownames(abund) <- abund$Site
abund$Site <- NULL
abund <- abund %>% dplyr::mutate(across(everything(), as.numeric))
abund <- abund %>% select(sort(current_vars()))
rm(NAlist)
rownames(traits2) <- traits2$Species

# calculate FD indices for each site with respect to the global functional trait space
# because this is with respect to all study sites, sites are comparable
# not calculating TRichness because we're using an alternative calculation
trait.dis2 <- cailliez(gowdis(traits2[-1], ord='podani'))
is.euclid(trait.dis2) # check whether distances are Euclidean before running PCoA

set.seed(24)
TS2 <- dbFD(trait.dis2, abund, w.abun=T, calc.FRic=T, calc.FDiv=T, m=5, calc.CWM=F, calc.FGR=F, print.pco=T)
tr.point <- as.data.frame(TS2$x.axes)[,1:4] # create a data frame of each species position in the 4D trait space
tr.point$Species <- rownames(traits2)

# extract the FD indices from the TS2 object, create a dataframe of them
predictors <- data.frame(Site=names(TS2$FEve[1:3]),TEve=TS2$FEve[1:3], TDiv=TS2$FDiv[1:3])
# only 7 for the 7 sites. Ignore the global FD measures
rownames(predictors) <- names(TS2$FEve[1:3])

# construct data frame of trait space points for each site
# and then that data frame serves as the input to the TOP metric function by Fontana et al. 2015
# Source: doi.org/10.1111/1365-2435.12551

require(geometry)
s.point <- as.list(rep(0,3))
predictors$TOP <- rep(0,3)
for (i in 1:3) {
  # first, make a dummy points dataframe only with the species in that site
  s.point[[i]] <- fish.ruv %>% filter(Site == unique(fish.ruv$Site)[i]) %>% ungroup() %>% select(Species)
  s.point[[i]] <- left_join(s.point[[i]], tr.point, by='Species')# use the point df to match the positions from the PCoA
  predictors$TOP[i] <- TOP.index(s.point[[i]][-1])[2] # calculate the TOP index
}

# now we have to standardise this by the global trait space TOP
predictors$TOP <- predictors$TOP/TOP.index(tr.point[1:4])[2]

# the TOP index function generates a vert.txt file but we don't need that anymore
if (file.exists('vert.txt')) {
  #Delete file if it exists
  invisible(file.remove('vert.txt'))
}

## Calculate relative abundances herbivore biters --------------------------

# reference vector of herbivore species
herb_spp <- traits2 %>% filter(TrophParr == 'herb_mic') %>% pull(Species)
site.herb <- fish.ruv %>% filter(Species %in% herb_spp) # filter assemblage data by herbivores
site.herb <- site.herb %>% group_by(Site) %>% summarise(Herb = sum(n)) %>% # add them all up
  left_join(., site.n, by="Site") %>% 
  mutate(Herb = Herb/n) %>% select(!n) # divide by total for rel abundance
# herbivore abundance

# biter abundance
# just pull out the sessile invertivore species
bite_spp <- traits %>% filter(str_detect(TrophParr, 'sess_inv')) %>% pull(Species)
site.biter <- fish.ruv %>% filter(Species %in% bite_spp) # filter assemblage data by herbivores
site.biter <- site.biter %>% group_by(Site) %>% summarise(Benthic = sum(n)) %>% # add them all up
  left_join(., site.n, by="Site") %>% 
  mutate(Benthic = Benthic/n) %>% select(!n) # divide by total for rel abundance

# add it to our predictors dataframe
predictors$Herb <- site.herb$Herb
predictors <- left_join(predictors, site.biter, by="Site")
predictors$Benthic[which(is.na(predictors$Benthic))] <- 0 # because some sites could have none, the NAs need to be zeroes

# clean up objects that aren't dependencies for downstream sourcing
rm(herb_spp, biter_spp, site.biter, site.herb, abund)

## Visualise trait space 2 only --------------------------

require(ggrepel)
require(patchwork)
require(viridis)
require(ggplot2)

t.space <- as.list(rep(0,4))

hull.v <- convhull.vert(tr.point[,1:2])
hull.v2 <- convhull.vert(tr.point[,3:4])

# palette for sites
palette <- viridis::viridis(n=3, end=0.95, begin=0.1, alpha=0.7, option="mako")
titles <- c('North Reef', 'Southeast', 'Turtle Beach')
looks <- theme_bw(base_size=12) + theme(panel.grid = element_blank())

# object for total fishes in each site


for (i in 1:3) {
  # Now that the base layers are done, I'll plot each site's points and hulls on
  # first, make a dummy points dataframe only with the species in that site
  pointS <- fish.ruv %>% filter(Site == unique(fish.ruv$Site)[i]) %>% ungroup() %>% select(Species, n)
  pointS$n <- pointS$n/site.n$n[i] # relative abundance of each species
  pointS <- left_join(pointS, tr.point[1:5], by='Species')
  pointS <- pointS[order(pointS$n, decreasing=T),order(colnames(pointS))]
  # make abbreviated species names
  shortsp <- as.data.frame(matrix(unlist(strsplit(pointS$Species[1:3], " ")), ncol=2, byrow=TRUE))
  pointS$ShortSp <- ''
  pointS$ShortSp[1:3] <- paste(abbreviate(shortsp[,1], 1, strict = TRUE), ". ", shortsp[,2], sep="")
  
  # now calculate the convex hull for that site only
  # global convex hull/TOP
  s.hull.v <- convhull.vert(pointS[,1:2])
  s.hull.v2 <- convhull.vert(pointS[,3:4])
  
  s.space <- as.list(c(0,0))
  
  s.space[[1]] <- ggplot() + looks +
    geom_polygon(data=hull.v, aes(x=A1, y=A2), color='grey', fill='white') +
    labs(x='PCo1', y='PCo2', title=titles[i]) +
    geom_polygon(data=s.hull.v, aes(x=A1, y=A2), fill=palette[i], alpha=0.25) +
    geom_point(data=pointS %>% left_join(., traits2, by="Species"),
               aes(x=A1, y=A2, size=n), fill=palette[i], color='#666666', shape=21) +
    geom_text_repel(data=pointS %>% arrange(n) %>% top_n(5), aes(label=ShortSp, x=A1, y=A2), force=2,
                    size=2.8, fontface='italic', min.segment.length = 0, point.size=5,
                    max.overlaps = Inf, box.padding = 1, force_pull=0.1) +
    scale_size_continuous(range=c(1.5,8), guide="none", limits=c(0,0.7)) +
    if (i != 4) {
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), 
            legend.position='none')
    } 
  else {theme()}
  s.space[[2]] <- ggplot() + looks +
    geom_polygon(data=hull.v2, aes(x=A3, y=A4), color='grey', fill='white') +
    labs(x='PCo3', y='PCo4') +
    geom_polygon(data=s.hull.v2, aes(x=A3, y=A4), fill=palette[i], alpha=0.25) +
    geom_point(data=pointS %>% left_join(., traits2, by="Species"), 
               aes(x=A3, y=A4, size=n), fill=palette[i], color='#666666', shape=21) +
    geom_text_repel(data=pointS %>% arrange(n) %>% top_n(5), aes(label=ShortSp, x=A3, y=A4), force=2,
                    size=2.8, fontface='italic', min.segment.length = 0, point.size=5,
                    max.overlaps = Inf, box.padding = 1, force_pull=0.1) +
    scale_size_continuous(range=c(1.5,8), limits=c(0,0.7), name="Relative abundance", guide="legend") + 
    if (i != 4) {
      theme(axis.title.y=element_blank(), axis.text.y=element_blank())
    }
  else {theme()}
  
  t.space[[i]] <- s.space[[1]] / s.space[[2]] # stack them with patchwork and store
}

rm(pointS)
tr.point$n <- Sp.count$n/sum(Sp.count$n) # put in global abundances to visualise all sites
global1 <- ggplot() + looks +
  geom_polygon(data=hull.v, aes(x=A1, y=A2), fill='transparent', color='grey') +
  geom_point(data=left_join(tr.point, traits2, by="Species"), aes(x=A1, y=A2, shape=ForageMode), size=2, fill='black', alpha=0.5) +
  labs(x='PCo1', y='PCo2', title='All sites') +
  scale_fill_grey(aesthetics="color",start=0.2, end=0.7) + scale_shape_manual(values=c(1:6,21:25), name="Foraging mode")

global2 <- ggplot() + looks +
  geom_polygon(data=hull.v2, aes(x=A3, y=A4), fill='transparent', color='grey') +
  geom_point(data=left_join(tr.point, traits2, by="Species"), aes(x=A3, y=A4, shape=ForageMode), size=2, fill='black', alpha=0.5) +
  labs(x='PCo3', y='PCo4') + theme(legend.title=element_blank()) +
  scale_fill_grey(aesthetics="color", start=0.2, end=0.7) + scale_shape_manual(values=c(1:6,21:25), guide="none")

pred_long <- predictors %>% pivot_longer(cols=!Site, names_to="predictor", values_to="value")

FD_bar <- ggplot(pred_long %>% filter(predictor %in% c('TOP', 'TEve', 'TDiv'))) +
  geom_bar(aes(y=value, x=predictor, fill=Site), stat='identity', position='dodge', color='black') +
  labs(y='Index measure', x=NULL) +
  scale_fill_viridis_d(begin=0.1, end=0.95, guide="none", option="mako") + looks +
  scale_y_continuous(expand=expansion(mult=c(.0,.05)), limits=c(0,1))

A_bar <- ggplot(pred_long %>% filter(predictor == 'Herb' | predictor == 'Benthic')) +
  geom_bar(aes(y=value, x=predictor, group=Site, fill=Site), stat='identity', position='dodge', color='black') +
  labs(y='Relative abundance', x=NULL) +
  scale_fill_viridis_d(begin=0.1, end=0.95, guide="none", option="mako") + looks + 
  scale_y_continuous(expand=expansion(mult=c(.0,.05))) +
  scale_x_discrete(labels=c('Benthic biters', 'Herbivores'))

bars <- (A_bar | FD_bar) * plot_layout(widths=c(2,3))

Fig2 <- (bars / ((global1/global2) | t.space[[1]] | t.space[[2]] | t.space[[3]])) * plot_layout(guides='collect', heights = c(1,3)) * theme(axis.text = element_text(size=9)) & theme(plot.title=element_text(face='bold', hjust=0.5, size=13))
Fig2
ggsave(filename = "../MS_CoralReefs/rev2/figures/traitcompare_indp.svg", device = "svg", width=30, height=16, units='cm', dpi=300)

rm(bars, global1, global2, t.space, A_bar, FD_bar, pred_long, s.space, s.hull.v, s.hull.v2, shortsp, site.n, looks, hull.v, hull.v2, Sp, palette, i, validate, titles) # clear figure objects



# Trait space comparison,  together ---------------------------------------

# so the section above compares the second videos in isolation
# see how the trait space construction differs from the original
# this comparison maps the original assemblage data with the new assemblage data from the second videos onto the same trait space

fish_assemblage <- read.csv('src/fish_assemblage.csv', header=T) %>% filter(!str_detect(Species, ' sp$'))
fish_sp <- read.csv('src/fish_sptaxonomy.csv', header=T)[-1]
ruv.sp <- Sp.count
rm(Sp.count, Sp.list, ruv2.data, s.point, tr.point, validate, Sp, trait.dis2) # clear the objects from previous section

# make a joined species total abundance table
sp.combined <- full_join(fish_sp, ruv.sp, by=c("family", "genus", "Species")) %>% arrange(Species)
# n.x = original RUV
# n.y = second video RUV
colnames(sp.combined)[4:5] <- c('RUV1', 'RUV2')

# now, making a joined assemblage data table
sites <- c('NR_2', 'SE_2', 'TB_2')
names(sites) <- unique(fish.ruv$Site)
fish.ruv$Site <- str_replace_all(fish.ruv$Site, sites)
ruv.combined <- bind_rows(fish_assemblage, fish.ruv)
str(ruv.combined)

# create the abundance matrix for functional evenness and divergence weighting
# species are columns, each site is a row
# long format to wide
abund <- sp.combined %>% filter(!str_detect(Species, ' sp$')) %>% 
  select(Species, RUV1) %>% 
  pivot_wider(names_from = Species, values_from = RUV1)
abund <- rbind(sp.combined %>% filter(!str_detect(Species, ' sp$')) %>% select(Species, RUV2) %>% pivot_wider(names_from = Species, values_from = RUV2), abund)
abund$Site <- c('RUV1', 'RUV2')
abund <- rbind(ruv.combined %>% filter(!str_detect(Species, ' sp$')) %>% select(Site, Species, n) %>% pivot_wider(names_from = Species, values_from = n), abund)
# because not all species occur in each site, there are lots of NAs to be filled with 0s
NAlist <- as.list(rep(0,dim(abund)[2]-1)) # list has to have column names for every item = 0
names(NAlist) <- colnames(abund)[2:ncol(abund)]
abund <- abund %>% replace_na(NAlist) # annnd replace NAs
colnames(abund)[1] <- '0Site'
abund <- abund %>% select(sort(tidyselect::peek_vars()))
abund[-1] <- abund[-1] %>% dplyr::mutate(across(everything(), as.numeric))
rm(NAlist)

# combined trait table
traits.combined <- traits2 %>% filter(!Species %in% traits$Species) %>% bind_rows(traits)
traits.combined <- traits.combined %>% arrange(Species)
rownames(traits.combined) <- traits.combined$Species

## Calculate trait space ---------------------------------------------------

# calculate FD indices for each site, each video group, with respect to the global functional trait space
dis.combined <- cailliez(gowdis(traits.combined[-1], ord='podani'))
is.euclid(dis.combined) # check whether distances are Euclidean before running PCoA

set.seed(24)
require(FD)
TS.combined <- dbFD(dis.combined, abund[-1], w.abun=T, calc.FRic=T, calc.FDiv=T, m=5, calc.CWM=F, calc.FGR=F, print.pco=T)
sp_points <- as.data.frame(TS.combined$x.axes)[,1:4] # create a data frame of each species position in the 4D trait space
sp_points$Species <- rownames(traits.combined)

N <- nrow(abund)-2 # number of sites total
# extract the FD indices from the object created by FD, create a dataframe of them
names(abund)[1] <- 'Site'
new_vars <- data.frame(Site = abund$Site[1:N], TEve = TS.combined$FEve[1:N], TDiv = TS.combined$FDiv[1:N])
rownames(predictors) <- names(TS.combined$FEve[1:N])


## Calculate TOP -----------------------------------------------------------

require(geometry)
site.point <- as.list(rep(0,N))
new_vars$TOP <- rep(0,N)
for (i in 1:N) {
  # first, make a dummy points dataframe only with the species in that site
  site.point[[i]] <- ruv.combined %>% filter(Site == unique(ruv.combined$Site)[i]) %>% ungroup() %>% select(Species)
  site.point[[i]] <- left_join(site.point[[i]], sp_points, by='Species')# use the point df to match the positions from the PCoA
  new_vars$TOP[i] <- TOP.index(site.point[[i]][-1])[2] # calculate the TOP index
}

# now we have to standardise this by the global trait space TOP
new_vars$TOP <- new_vars$TOP/TOP.index(sp_points[1:4])[2]

# the TOP index function generates a vert.txt file but we don't need that anymore
if (file.exists('vert.txt')) {
  #Delete file if it exists
  invisible(file.remove('vert.txt'))
}


# ## Comparison indices ---------------------------------------------------

# Bray-Curtis dissimilarity
# North Reef
vegdist(abund[c(3,8),-1], method="bray", binary=F)
# Southeast
vegdist(abund[c(5,9),-1], method="bray", binary=F)
# Turtle Beach
vegdist(abund[c(6,10),-1], method="bray", binary=F)

# trait diversity differences
tdelta <- unlist(new_vars[3,-1] - new_vars[8,-1])
tdelta <- c(unlist(new_vars[5,-1] - new_vars[9,-1]), tdelta)
tdelta <- c(unlist(new_vars[6,-1] - new_vars[10,-1]), tdelta)
mean(tdelta)


## Visualise assemblage trait space differences ----------------------------

require(ggrepel)
require(tidyverse)
require(patchwork)
require(viridis)

t.space <- as.list(rep(0,3))
source('./analysis_code/function_convhull.R') # make plottable convex hull vertices

hull.v <- convhull.vert(sp_points[,1:2]) # calculate convex hull vertices for the first 2 PCoA dimensions
hull.v2 <- convhull.vert(sp_points[,3:4]) # calculate verts for the next 2 dimensions

# global figure aesthetics
palette <- viridis::viridis(n=7, end=0.95, begin=0.1, alpha=0.7, option="mako") # palette only for original data
looks <- theme_bw(base_size=12) + theme(panel.grid = element_blank())
detach(package:plyr)
site.n <- ruv.combined %>% group_by(Site) %>% summarise(total=sum(n))

sites = list(c(3,8), c(5,9), c(6,10))
for (i in 1:3) {
  # first, make a dummy points dataframe only with the species in the original site
  pointS <- ruv.combined %>% filter(Site %in% abund$Site[sites[[i]][1]]) %>% ungroup() %>% select(Species, n)
  pointS$n <- pointS$n/site.n$total[sites[[i]][1]] # relative abundance of each species
  pointS <- left_join(pointS, sp_points, by='Species') # add their PCoA coordinates
  pointS <- pointS[order(pointS$n, decreasing=T),order(colnames(pointS))] 
  # make abbreviated species names just for the 3 most abundant species
  shortsp <- as.data.frame(matrix(unlist(strsplit(pointS$Species[1:3], " ")), ncol=2, byrow=TRUE))
  # split genus and species into different columns
  pointS$ShortSp <- ''
  pointS$ShortSp[1:3] <- paste(abbreviate(shortsp[,1], 1, strict = TRUE), ". ", shortsp[,2], sep="") 
  # shorten into G. species format
  
  # now calculate the convex hull
  hull <- convhull.vert(pointS[,1:2])
  hull_2 <- convhull.vert(pointS[,3:4])
  
  # rinse and repeat for the second site
  pointS2 <- ruv.combined %>% filter(Site %in% abund$Site[sites[[i]][2]]) %>% ungroup() %>% select(Species, n)
  pointS2$n <- pointS2$n/site.n$total[sites[[i]][2]] # relative abundance of each species
  pointS2 <- left_join(pointS2, sp_points, by='Species') # add their PCoA coordinates
  pointS2 <- pointS2[order(pointS2$n, decreasing=T),order(colnames(pointS2))] 
  # make abbreviated species names just for the 3 most abundant species
  shortsp <- as.data.frame(matrix(unlist(strsplit(pointS2$Species[1:3], " ")), ncol=2, byrow=TRUE))
  # split genus and species into different columns
  pointS2$ShortSp <- ''
  pointS2$ShortSp[1:3] <- paste(abbreviate(shortsp[,1], 1, strict = TRUE), ". ", shortsp[,2], sep="") 
  # shorten into G. species format
  
  # now calculate the convex hull
  hull2 <- convhull.vert(pointS2[,1:2])
  hull2_2 <- convhull.vert(pointS2[,3:4])
  
  s.space <- as.list(c(0,0))
  
  s.space[[1]] <- ggplot() + looks +
    geom_polygon(data=hull.v, aes(x=A1, y=A2), color='grey', fill='white') +
    labs(x='PCo1', y='PCo2', title=abund$Site[sites[[i]][1]]) +
    geom_polygon(data=hull, aes(x=A1, y=A2), fill=palette[sites[[i]][1]], alpha=0.25) +
    geom_point(data=pointS %>% left_join(., traits.combined, by="Species"),
               aes(x=A1, y=A2), fill=palette[sites[[i]][1]], color='#666666', shape=21) +
    geom_text_repel(data=pointS %>% arrange(n) %>% top_n(5), aes(label=ShortSp, x=A1, y=A2), force=2,
                    size=2.8, fontface='italic', min.segment.length = 0, point.size=5,
                    max.overlaps = Inf, box.padding = 1, force_pull=0.1) +
    geom_polygon(data=hull2, aes(x=A1, y=A2), fill='orange', alpha=0.25) +
    geom_point(data=pointS2 %>% left_join(., traits.combined, by="Species"),
               aes(x=A1, y=A2), fill='orange', color='#666666', shape=21) +
    geom_text_repel(data=pointS2 %>% arrange(n) %>% top_n(5), aes(label=ShortSp, x=A1, y=A2), force=2,
                    size=2.8, fontface='italic', min.segment.length = 0, point.size=5,
                    max.overlaps = Inf, box.padding = 1, force_pull=0.1)
    
  s.space[[2]] <- ggplot() + looks +
    geom_polygon(data=hull.v2, aes(x=A3, y=A4), color='grey', fill='white') +
    labs(x='PCo3', y='PCo4') +
    geom_polygon(data=hull_2, aes(x=A3, y=A4), fill=palette[sites[[i]][1]], alpha=0.25) +
    geom_point(data=pointS %>% left_join(., traits.combined, by="Species"), 
               aes(x=A3, y=A4), fill=palette[sites[[i]][1]], color='#666666', shape=21) +
    geom_text_repel(data=pointS %>% arrange(n) %>% top_n(5), aes(label=ShortSp, x=A3, y=A4), force=2,
                    size=2.8, fontface='italic', min.segment.length = 0, point.size=5,
                    max.overlaps = Inf, box.padding = 1, force_pull=0.1) +
    geom_polygon(data=hull2_2, aes(x=A3, y=A4), fill='orange', alpha=0.25) +
    geom_point(data=pointS2 %>% left_join(., traits.combined, by="Species"), 
               aes(x=A3, y=A4), fill='orange', color='#666666', shape=21) +
    geom_text_repel(data=pointS2 %>% arrange(n) %>% top_n(5), aes(label=ShortSp, x=A3, y=A4), force=2,
                    size=2.8, fontface='italic', min.segment.length = 0, point.size=5,
                    max.overlaps = Inf, box.padding = 1, force_pull=0.1)
    
  t.space[[i]] <- s.space[[1]] / s.space[[2]] # stack them with patchwork and store
}

rm(pointS, pointS2, hull, hull_2, hull2, hull2_2)

# do something similar for the global trait spaces to compare RUV1 with RUV2
pointS <- sp.combined %>% select(Species, n=RUV1) %>% na.omit()
pointS$n <- pointS$n/624 # relative abundance of each species
pointS <- left_join(pointS, sp_points, by='Species') # add their PCoA coordinates
pointS <- pointS[order(pointS$n, decreasing=T),order(colnames(pointS))] 
# make abbreviated species names just for the 3 most abundant species
shortsp <- as.data.frame(matrix(unlist(strsplit(pointS$Species[1:3], " ")), ncol=2, byrow=TRUE))
# split genus and species into different columns
pointS$ShortSp <- ''
pointS$ShortSp[1:3] <- paste(abbreviate(shortsp[,1], 1, strict = TRUE), ". ", shortsp[,2], sep="") 
# shorten into G. species format

hull <- convhull.vert(pointS[,1:2])
hull_2 <- convhull.vert(pointS[,3:4])

# rinse and repeat for the second site
pointS2 <- sp.combined %>% select(Species, n=RUV2) %>% na.omit()
pointS2$n <- pointS2$n/392 # relative abundance of each species
pointS2 <- left_join(pointS2, sp_points, by='Species') # add their PCoA coordinates
pointS2 <- pointS2[order(pointS2$n, decreasing=T),order(colnames(pointS2))] 
# make abbreviated species names just for the 3 most abundant species
shortsp <- as.data.frame(matrix(unlist(strsplit(pointS2$Species[1:3], " ")), ncol=2, byrow=TRUE))
# split genus and species into different columns
pointS2$ShortSp <- ''
pointS2$ShortSp[1:3] <- paste(abbreviate(shortsp[,1], 1, strict = TRUE), ". ", shortsp[,2], sep="") 
# shorten into G. species format

hull2 <- convhull.vert(pointS2[,1:2])
hull2_2 <- convhull.vert(pointS2[,3:4])

global1 <- ggplot() + looks +
  geom_polygon(data=hull, aes(x=A1, y=A2), fill='grey80', color='grey20') +
  geom_point(data=left_join(pointS, traits2, by="Species"), aes(x=A1, y=A2), size=2, shape=21, fill='black', alpha=0.5) +
  geom_polygon(data=hull2, aes(x=A1, y=A2), fill='orange', color='orange', alpha=0.5) +
  geom_point(data=left_join(pointS2, traits2, by="Species"), aes(x=A1, y=A2), size=2, shape=21, fill='orange', alpha=0.5) +
  labs(x='PCo1', y='PCo2', title='RUV1 v RUV2')

global2 <- ggplot() + looks +
  geom_polygon(data=hull_2, aes(x=A3, y=A4), fill='grey80', color='grey20') +
  geom_point(data=left_join(pointS, traits2, by="Species"), aes(x=A3, y=A4), size=2, shape=21, fill='black', alpha=0.5) +
  geom_polygon(data=hull2_2, aes(x=A3, y=A4), fill='orange', color='orange', alpha=0.5) +
  geom_point(data=left_join(pointS2, traits2, by="Species"), aes(x=A3, y=A4), size=2, shape=21, fill='orange', alpha=0.5) +
  labs(x='PCo3', y='PCo4')

((global1 / global2) | t.space[[1]] | t.space[[2]] | t.space[[3]])
ggsave(filename = "../MS_CoralReefs/rev2/figures/traitcompare_together.svg", device = "svg", width=28, height=12, units='cm', dpi=300)

pred_long <- new_vars %>% pivot_longer(cols=!Site, names_to="predictor", values_to="value")
pred_long <- pred_long %>% filter(Site %in% abund$Site[c(3,5,6,8:10)])
pred_long$Site <- factor(pred_long$Site, levels = c('North3', 'Southeast', 'TurtleBeach', 'NR_2', 'SE_2', 'TB_2'))
ggplot(pred_long) +
  geom_bar(aes(y=value, x=predictor, fill=Site), stat='identity', position='dodge', color='black') +
  labs(y='Index measure', x=NULL) + theme_bw() +
  scale_y_continuous(expand=expansion(mult=c(.0,.05)), limits=c(0,1))

ggsave(filename = "../MS_CoralReefs/rev2/figures/traitindex_compare.eps", device = "eps", width=21, height=7, units='cm', dpi=300)


rm(bars, global1, global2, t.space, A_bar, FD_bar, pred_long, s.space, s.hull.v, s.hull.v2, shortsp, looks, hull.v, hull.v2, Sp, palette, i, validate, titles) # clear figure objects

