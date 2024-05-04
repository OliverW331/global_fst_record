# Script Repository README

## Overview
This repository contains all scripts used in the analysis, modeling, and simulation of the Global Species First Record Database. Each script has a specific role, from data processing and model fitting to sensitivity analysis and graphical representation.

## Script Descriptions

### Model Construction and Fitting
- **process_model_logistic_I_sigmoid_pd.R**: Constructs logistic growth models for invasive species detection and fits them to real data from the `data` repository.

### Graphical Analysis
- **global_fst_record.R**: Performs graphical analysis based on ecoregion, taxa, and country levels. Generates and saves graphs in the `graph` repository.

### Database Comparison and Analysis
- **delta_F.R**: Examines differences between two versions of the database from 1995 to 2016, aiming to identify trends suggesting a decline in actual invasions.
- **sensitivity_analysis.R**: Analyzes the correlation between the slopes of the differences in the two database versions and the slopes of actual invasions, expecting a positive correlation or at least a demonstration that an increase in invasions cannot explain the flat trend observed.
- **deltaF_I_cor.R**: Demonstrates that there is a positive correlation between the slopes of actual invasions (I) and the differences (delta F) if the invasion rate is linear.

### Sensitivity and Parameter Analysis
- **sensitivity_analysis_params.R**: Explores how variations in model parameters affect the modeled curves. Results are visualized in `model_graph/sensitivity_analysis`.

### Simulation and Model Validation
- **fit_simulations.R** & **fit_simulations.Rmd**: These scripts assess model fit by generating numerous simulations with random parameter combinations, comparing the fit of detection and reporting probability parameters versus actual invasion parameters. Results help evaluate the robustness of the model fitting process.
- **likelihood_check.R**: Investigates the likelihood consistency within the global species invasion model when detection and reporting processes are integrated, particularly focusing on the likelihood of species being recorded (invaded, detected and reported) or not.

## Improvements Needed
Some scripts require optimization and enhancement:
- **Parallelization**: Certain scripts are not yet parallelized, which could improve their efficiency and reduce runtime.
- **Sensitivity Analysis Range**: The selection ranges for sensitivity analysis need refinement to better demonstrate correlations and to enhance the interpretability of the findings. For example, although the correlation between slopes of difference and slopes of invasion from 1995 to 2016 is not strong, a negative or zero slope of difference occurs much more frequently when the slope of invasion is negative. The ranges for each parameter used to conduct sensitivity analysis and the function for invasion_growth need to be examined since we need to make sure the number of species recorded is always positive (which means that the rate of decline cannot be too large). Before conducting the sensitivity analysis, I should check if there is any anomaly with the combinations of parameters.

## Usage
To run these scripts, ensure you have R installed along with necessary libraries such as `ggplot2`, `dplyr`, `tidyr`, `gridExtra`, `foreach`, and `doParallel`. Adjust parameters and paths as needed for your specific setup.

## Contributions
Contributions to improve the scripts are welcome. Please fork this repository, make your changes, and submit a pull request.

## Contact
For more information or assistance with these scripts, please contact [Muchen Wang].

