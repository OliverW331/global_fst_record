# Detection and Reporting Model of Global Species First Record Database

## Introduction
The ability to detect and report invasive species accurately and promptly is critical in managing biodiversity and preventing the adverse effects of invasions. The detection and reporting model outlined here is designed to estimate the number of new invasions of species detected and reported over time, incorporating the probabilities associated with detection delays and reporting lags. This model is particularly valuable for ecological researchers and policymakers who need to understand the dynamics of species invasions and their impacts on ecosystems.

## Formula Description
The formula for estimating the number of new invasions detected in a particular year (\( T_d \)) and reported by a subsequent year (\( T_r \)) is given by:

\[ F(T_d) = \sum_{T_i = t_0}^{T_d} \left( I(T_i) \cdot \left( \prod_{j=0}^{T_d-T_i-1} (1-P_d(j)) \right) \cdot P_d(T_d-T_i) \cdot \left(1- \prod_{m=0}^{T_r-T_d} (1- P_r(m)) \right) \right) \]

### Process
- **Initial Invasion**: Accounting for invasions starting from year \( T_i \).
- **Non-Detection Period**: The period from the year of invasion up to just before detection.
- **Detection Year**: The specific year where detection occurs.
- **Reporting Delay**: The delay from detection to reporting.

## Variables Explanation
- **\( T_i \)**: Year of invasion.
- **\( T_d \)**: Year of detection.
- **\( T_r \)**: Year of reporting.
- **\( I(T_i) \)**: Number of invasions occurring in year \( T_i \).
- **\( P_d(j) \)**: Probability of detecting an invasion \( j \) years after invasion.
- **\( P_r(m) \)**: Probability of not reporting a detected invasion \( m \) years after detection.

These variables collectively contribute to modeling the detection and reporting processes, enabling us to simulate and predict the dynamics of invasive species reports over time. Understanding these dynamics helps in assessing the effectiveness of monitoring systems and in improving strategies for managing invasive species.

## Repositories Description

### 1. Data
The `data` repository houses three versions of the first record database for invasive species, available in both CSV and XLSX formats. These datasets are fundamental for all analysis and modeling tasks in the project.

### 2. Graph
The `graph` repository includes preliminary graphical analyses of the database. These graphs provide initial insights and visualizations of the data, helping to understand trends, distributions, and potential outliers in the invasive species records.

### 3. Model
Located in the `model` repository are the preliminary models fitted to the data. These models are constructed using the script `process_model_logistic_I_sig_pd_sig_pr.R` found in the `script` repository. They are crucial for predicting and understanding the dynamics of invasive species spread.

### 4. Model Graph
The `model_graph` repository contains analytical graphs generated to evaluate how well the models fit the data. These graphs are essential for visualizing the effectiveness of the models and for identifying areas where model adjustments may be necessary.

### 5. Model Result
This repository, `model_result`, includes the results of sensitivity analyses, simulations, and other outputs from model evaluations. It serves as a repository for storing detailed outcomes of various analyses that test the robustness and reliability of the models.

### 6. Script
The `script` repository encompasses all the scripts used to perform simulations, sensitivity analyses, and fitting of models to the database. This includes scripts for data preprocessing, analysis, and post-analysis evaluations.

## Usage
To use this project, clone the repositories and follow the setup instructions in each repository's specific README file (if available). Ensure you have the necessary software and libraries installed, as listed below:

- R software
- Required R packages: `ggplot2`, `dplyr`, `tidyr`, `gridExtra`, `foreach`, `doParallel`

## Contributing
Contributions to this project are welcome. To contribute, please fork the repository, make your changes, and submit a pull request. For substantial changes, please open an issue first to discuss what you would like to change.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments
Thanks to all contributors and researchers who have provided insights and feedback on the models and analyses.

## Contact
For any inquiries, please contact [Your Name or Project Email].

---

Please ensure to keep the documentation up-to-date as the project evolves.


