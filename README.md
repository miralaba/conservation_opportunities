# Achieving net-positive outcomes in Amazonia

## Introduction

This repository hosts pre-processed data and R scripts for the analysis presented in the manuscript *"Combined interventions required to reverse biodiversity and carbon losses in tropical forest frontiers"*. The project compares the biodiversity and carbon benefits and costs of three major conservation strategies—avoiding deforestation, avoiding degradation, and promoting restoration—across two regions of the Brazilian Amazon using remote sensing, ecological field data, and socioeconomic assessments.

## Repository Structure

- `data/`: Pre-processed datasets derived from publicly available sources
  - `carbon.csv`: Field-based carbon estimates
  - `presence_records.csv`: Presence-only species records (birds and plants)
  - `property.csv`: Information on rural property boundaries
  - `species_summary.csv`: Summary data used in biodiversity modelling
  - `transect_harvest_value.csv`: Timber value estimates from field transects

- `scripts/`: R scripts used to prepare data layers, model biodiversity and carbon, and run cost-benefit analyses
  - `auxiliar.R`: Functions used across scripts
  - `cost_benefit_analysis.R`: Main script for analysing cost-effectiveness
  - `layer_buide_pgm.R`: Spatial processing for Paragominas region
  - `layer_buide_stm.R`: Spatial processing for Santarém region
  - `regional_model.R`: Models for biodiversity and carbon predictions
  - `variables_adjustments.R`: Variable cleaning and standardisation

- `shapes/`: Spatial masks and shapefiles
  - Files include spatial boundaries for `Paragominas` and `Santarém` study regions
  - File types: `.shp`, `.shx`, `.dbf`, `.prj`, `.cpg`, `.qix`, `.xml`

## Prerequisites

To run the scripts, you will need:

- **R** version 4.3.0
- R packages:
  - `tidyverse`
  - `raster` version 3.6.20
  - `sf` version 1.0.13
  - `sp` version 1.6.1
  - `virtualspecies` version 1.5.1
  - `usdm` version 1.1.18
  - `randomForest` version 4.7.1.1
  - `caret` version 6.0.94
  - `spThin` version 0.2.0
  - `ggplot2`

## Install all required packages using:

install.packages(c("tidyverse", "raster", "sp", "sf", "virtualspecies", "usdm", "randomForest", "caret", "spThin", "ggplot2"))


## Usage

1. Clone this repository: 
  git clone https://github.com/miralaba/conservation_opportunities.git

2. Open RStudio and load the necessary packages.

3. Run the scripts in the following suggested order:

  - layer_buide_pgm.R and layer_buide_stm.R: Process spatial layers for each region
  - variables_adjustments.R: Prepare input variables
  - regional_model.R: Run biodiversity and carbon predictions
  - cost_benefit_analysis.R: Calculate benefits, costs, and cost-effectiveness
  - Use auxiliar.R functions as needed in other scripts

## Raw Data Sources

This repository uses public datasets processed into standardised formats. Original raw data can be accessed via:

    MapBiomas Collection 8 – Land Use and Land Cover: https://brasil.mapbiomas.org/en/colecoes-mapbiomas/

    Silva Jr. et al. 2020 – Secondary forest age: https://github.com/celsohlsj/gee_brazil_sv

    INPE DEGRAD and DETER-B – Degradation monitoring:

        DEGRAD: http://www.obt.inpe.br/OBT/assuntos/programas/amazonia/degrad/acesso-ao-dados-do-degrad

        DETER-B: http://terrabrasilis.dpi.inpe.br/file-delivery/download/deter-amz/shape

    MapBiomas Fire Collection 2 – Fire scars: https://brasil.mapbiomas.org/en/colecoes-mapbiomas/

    NASA Earth Observations – Climate: https://neo.gsfc.nasa.gov/

    SRTM (USGS) – Elevation: https://earthexplorer.usgs.gov/

    HydroSHEDS – Rivers: https://www.hydrosheds.org/

    Imazon – Road network (available on request)

    SISCAR – Rural property data: https://www.car.gov.br/  

## Contact

For questions, please contact:

    Leonardo Miranda: miralaba@gmail.com