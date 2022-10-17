
# Execution script for all analysis and figure generation components
# Author: Cher Chow
# Last version check: 13 Oct 2021


# Analysis ----------------------------------------------------------------

source('./analysis_code/A01_traitspace.R') # construct trait spaces for each site
# also calculates trait diversity predictors for models
# expect computation time of a few minutes, longest part of analysis
source('./analysis_code/A02_foragingrates.R') # processing of weighted foraging rate predictors
source('./analysis_code/A03_model_sett.R') # settlement model constructions + selection
source('./analysis_code/A04_model_rec.R') # recruitment model constructions + selection


# Figure generation -------------------------------------------------------

source('./figure_code/F02_traitspace.R') # trait space visualisation with barplot of predictors
source('./figure_code/F03_foraging.R') # foraging rates by trophic and functional group
source('./figure_code/F04_settlement.R') # settlement partial regression plot
source('./figure_code/F05_recruitment.R') # recruitment partial regression plot