# MS_iSTAY

**This repository contains three datasets and two sets of R source code associated with the paper "A continuum of information-based temporal stability measures and their decomposition across hierarchical levels" by Chao et al. (2026).**

## Datasets (based on the Jena Experiment)

### Data_Jena_20_metacommunities

Biomass time-series data from 2003 to 2024 (22 years) for the constituent communities within each of the 20 metacommunities in the Jena Experiment. The dataset contains a list of 20 data frames, each with plots (communities) as rows and time points as columns.

All 76 plots (communities) are grouped into 20 metacommunities according to block identity and species-richness level. Each metacommunity consists of three or four plots with the same species richness but different species compositions.

### Data_Jena_462_populations

Biomass time-series data for 462 populations across 76 plots from 2003 to 2024 (21 years, excluding 2004). The dataset is a population-by-time data frame with 462 rows and 21 columns. Each row represents a population, and each column represents one year of biomass data.

### Data_Jena_76_community_populations

Biomass time-series data for 76 communities in the Jena Experiment. Each community dataset contains the 21-year biomass time-series data for all species sown within each plot. Because the number of species sown per plot ranges from 1 to 16, the number of populations within a community also ranges from 1 to 16.

This dataset contains the same data as `Data_Jena_462_populations`, but is organized for analytical convenience into a list of 76 communities, each containing a population-by-time biomass data frame for its constituent populations.

All three datasets are provided in the `data/` folder.

## R Source Code

### Main_Document_R code.R

R code for performing the analyses and generating all figures (Figures 1–4) in the main text of Chao et al. (2026).

### Supporting_Information_R_code.R

R code for performing the analyses and generating all figures in Appendices S1–S5 of Chao et al. (2026).

The R scripts contain labeled sections for locating the code corresponding to each figure. Generated figures are saved in the `figures/` folder.

## Access Information

The datasets and associated analytical functions are also available through the R package **iSTAY** (information-based stability and synchrony measures), available from CRAN:

https://cran.r-project.org/web/packages/iSTAY/index.html

An online version is available at:

https://chao.shinyapps.io/iSTAY/

The biomass data from 2003 to 2023 were retrieved from:

https://jexis.idiv.de/ddm/Data/ShowData/624

The 2024 biomass data were obtained from:

https://jexis.idiv.de/ddm/Data/ShowData/695

## Reproducing the Figures

Open `MS_iSTAY.Rproj` in RStudio and run `Main_Document_R code.R` or `Supporting_Information_R_code.R` to reproduce the corresponding figures.

The datasets can be loaded directly from the `data/` folder.