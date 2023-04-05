
# FishTraitsCoralRec
### _Coral settlement and recruitment are negatively related to reef fish trait diversity_ 
**DOI**: [10.1007/s00338-023-02359-7](https://doi.org/10.1007/s00338-023-02359-7) in _Coral Reefs_ March 2023

**Authors**: Cher F Y Chow, Caitlin Bolton, Nader Boutros, Viviana Brambilla, Luisa Fontoura, Andrew S Hoey, Joshua S Madin, Oscar Pizarro, Damaris Torres-Pulliza, Rachael M Woods, Kyle J A Zawada, Miguel Barbosa, Maria Dornelas 

**Corresponding author**: Please direct any inquiries or questions to Cher via [Github](https://github.com/cherfychow) or [email](mailto:cher.fyc@gmail.com)    
  
This repository contains the data and code used for the analysis and figure generation for _Coral settlement and recruitment are negatively related to reef fish trait diversity_. The aim of this study was to investigate whether variation in foraging trait diversity of fish assemblages correlates with variation in coral settlement and subsequent recruitment to juvenile cohorts. We used various metrics to describe fish assemblage foraging impact on corals, both directly (i.e. observed foraging rates) and indirectly via foraging traits (trait richness, trait evenness, trait divergence, herbivore abundance, and benthic biter abundance).
  
**Disclaimer**: This repository has been archived in Zenodo on its initial release at publication, however the code may undergo revisions and changes following releases.

## Requirements
All files were written in R 4.0.0. We recommend executing the repository using this version along with the relevant versions of the packages below.  
**Package dependencies:**  
The package dependencies are declared in each file as `require(package name here)`
- `tidyverse` (v 1.3.1)
- `MuMIn`
- `lubridate`
- `worms`
- `rfishbase`
- `FD` (v 1.0-12)
- `patchwork` (v 1.1.1) * only for replicating visualisations
- `geometry` (v 0.4.5)
- `performance` (v 0.7.2)

## Repository structure
**`exec.R`** : This is the execution file which handles the entire analysis pipeline.
**`analysis_code`** : Contains analysis scripts separated by analysis component/primary data used. Not all scripts are directly used in the pipeline, such as reference scripts used for cleaning, diagnostic checks or extracting tables from FishBase. 
**`figure_code`** : Contains figure generation code. Dependent on object outputs from analysis code.
**`src`** : This directory contains all relevant data for fish and coral assemblages.

## Data
- **`coral_recruit.csv`** :  data on coral recruit counts
- **`coral_settlement.csv`** :  data on coral spat counts
- **`fish_assemblage_vid.csv`** :  fish abundance by site with video file information 
- **`fish_assemblage.csv`** :  fish abundance by site
- **`fish_bites.csv`** :  fish bite rate data with length class and screen duration fields
- **`fish_camera2.csv`** :  subset of secondary video data for sensitivity analyses
- **`fish_sptaxonomy.csv`** :  species abundances across all sites
- **`fish_traits.csv`** :  trait table for all observed fish species
- **`fish_traits2.csv`** :  trait table for sensitivity analyses (species in secondary videos)
- **`fish_traitsinfluence.csv`**:  adapted trait table according to bite influence for trait-weighted foraging rates

## Analysis
A large part of the analysis pipeline is in calculating predictor variables relating to fish diversity and foraging. The final steps of the pipeline fit the coral settlement and recruitment models as well as performing diagnostics and model selection.

- **`A01_traitspace.R`** Trait space construction and calculation of trait diversity indices from fish assemblage abundance and trait tables.
- **`A01b_pcoachecks.R`** Diagnostic checks of the principal coordinates analysis that the trait space is based on.
- **`A02_foragingrates.R`**: Calculation of trait-weighted foraging rates that takes into account length and trait differences in benthic bite impact
- **`A03_model_sett.R`**: Model construction and selection for coral settlement
- **`A04_model_rec.R`**: Model construction and selection for coral recruitment.
- **`A05_replicatesensitivity.R`**: Sensitivity analyses that incorporates secondary video data for a subset of sites to validate spatial sampling effort
- **`A06_SACtest.R`**: Species accumulation curve sensitivity analysis that tests for sufficient temporal sampling effort


## Figure generation
Figure code is largely self explanatory and roughly follows the structure of the analyses code files listed above. Because we used GLMMs, note that figures 4 and 5 will require longer computational time for bootstrapped confidence intervals.


