
# Execution script for all analysis and figure generation components
# Author: Cher Chow


# Analysis ----------------------------------------------------------------

source('./analysis_code/A01_data_cleanmerge.R') # sets up data tables
source('./analysis_code/A02_traitspace.R') # construct trait spaces for each site
# also calculates trait diversity predictors for models
source('./analysis_code/A03_foragingrates.R') # processing of weighted foraging rates
source('./analysis_code/A04_model_sett.R') # settlement model constructions + selection
source('./analysis_code/A05_model_rec.R') # recruitment model constructions + selection


# Figure generation -------------------------------------------------------

source('./figure_code/F02_traitspace.R') # trait space visualisation with barplot of predictors
source('./figure_code/F03_foraging.R') # foraging rates by trophic and functional group
source('./figure_code/F04_settlement.R') # settlement partial regression plot
source('./figure_code/F05_recruitment.R') # recruitment partial regression plot