# Ecological scales of effect vary across space and time

### [Brent S. Pease](https://peaselab.com)

### Data/code DOI: [![DOI](https://zenodo.org/badge/794152577.svg)](https://zenodo.org/doi/10.5281/zenodo.11094099)

#### Please contact Brent Pease (bpease1@siu.edu) for questions about the code or data.
__________________________________________________________________________________________________________________________________________

## Abstract:  
The spatial scale at which an environmental variable is summarized can have considerable impacts on ecological inference of species distribution and abundance. While several analytical approaches have emerged to determine biologically relevant spatial scales— the spatial scale that most strongly influences the ecological patterns observed —identifying key ecological drivers of scale of effect is still underway. Additionally, several predicted ecological drivers of scale of effect can vary across space and time, but little research on spatiotemporal patterns has occurred. Here, I assessed spatial and temporal variation in scales of effect across 32 North American bird species using 18 years of empirical data from the North American Breeding Bird Survey. Scale estimation was then coupled with trait-based analyses and hypotheses testing of underlying processes of spatial and temporal variation in scales of effect. All 32 species tested exhibited varied scales of effect across years (average annual scales of effect ranging from 0.2 to 4.97 km) and Bird Conservation Regions, with spatial variability being the most pronounced. Trait-based analyses revealed a contrary relationship between hand-wing index, body size, and scale of effect, though the strength of this relationship was contingent on migratory status. Temporal variation in scales of effect was best explained by changes in human development over time, indicating that avian space use can be explained by an increasing human footprint. Additionally, relative population size, consistent with theoretical predictions stemming from density-dependent population dynamics, explained notable variation in spatial and temporal scales of effect. These findings contribute to the growing landscape ecology literature by providing empirical evidence for hypothesized drivers of scales of effect. By delineating species-specific scales of effect and elucidating their ecological drivers, this study enhances our understanding of spatial and temporal scales in ecological processes, aiding conservation efforts in a rapidly changing world.


## Repository Directory

### [Scripts](./Scripts): Contains code for cleaning, processing, and analyzing survey and supplemental data
* 101 - Data Prep - Get Species Observations.R - script used to process raw BBS data into format for analysis.
* 102 - Data Prep - Get MODIS Data.R - script used for downloading and process MODIS data for landcover analysis.
* 201 - Data Analysis - Est Scale for Species by Year.R - 'master' analysis script. This was used to calculate species- and year-specific scales of effect using single-species, single-season occupancy models in package `unmarked`. Manuscript methods subsection `Scales of Effect`.
* 202 - Data Analysis - Est Scale for Species by BCR.R - 'master' analysis script. This was used to calculate species- and BCR-specific scales of effect using single-species, single-season occupancy models in package `unmarked`. Manuscript methods subsection `Scales of Effect`.
* 203 - Data Analysis - Test SoE Temporal Driver Hypotheses.R - This script tested competing hypotheses of scale of effect drivers. Manuscript methods subsection `Ecological Drivers of Scale of Effect`.
* 204 - Data Analysis - Test SoE Spatial Driver Hypotheses.R - This script tested competing hypotheses of scale of effect drivers. Manuscript methods subsection `Ecological Drivers of Scale of Effect`.
* 205 - Data Analysis - Fragmentation Hypothesis.R - Intermediate processing script for spatiotemporal drivers of scales of effect.
* 206 - Data Analysis - Spp Rich Hypothesis.R - Intermediate processing script for spatiotemporal drivers of scales of effect.
* 207 - Data Analysis - Forest Hypothesis.R - Intermediate processing script for spatiotemporal drivers of scales of effect.
* 208 - Data Analysis - Development Hypothesis.R - Intermediate processing script for spatiotemporal drivers of scales of effect.
* 208 - Data Analysis - Development Hypothesis.R - Intermediate processing script for spatiotemporal drivers of scales of effect.
* 209 - Trait-based Analysis.R - This script was for assessing trait-based relationships with scales of effect. Manuscript methods subsection `Trait-based Relationships`.
* 301 - Figures - SoE over Space.R - Plotting spatial variation in scales of effect.
* 302 - Figures - SoE over Time.R - Plotting temporal variation in scales of effect.
* 401 - Summaries - Mean SoE.R - Table-form for estimated scales of effect.


### [Data](./Data): Contains raw and processed data used in scripts
BBS - Primary observation folder used in analyses. Uploaded directly to Dryad and not on github due to size
AVONET - trait database used in trait-based analysis
Global-HWI-v1.1 - hand-wing index database used in trait-based analysis

### [Results](./Results): Contains figures and output - will be created upon execution of code