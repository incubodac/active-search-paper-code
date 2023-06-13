# Concurrent EEG and Eyetracker Data Analysis for Active Search Experiment

This repository contains the code for preprocessing and analyzing concurrent EEG and eyetracker data from a visual search experiment. The code generates the outputs presented in the Active Search paper.

## Repository Structure

The repository is organized into the following folders:

- **preprocessing**: Contains the code for preprocessing the raw EEG and eyetracker data. The preprocessing steps include artifact removal, filtering, referencing, and segmentation.
- **analysis**: Contains the code for analyzing the preprocessed data. The analysis code includes predictor extraction, deconvolution analyses and statistical analyses, and visualization.

## Preprocessing

The `preprocessing` folder contains the following files:

- `dep/preAnalysis_dac.m`: Matlab script for preprocessing the raw EEG data.
- `res/FP_basePreAnalysis.m`: Matlab class for preprocessing the raw eyetracker data.
- `res/FP_dac.m`:  Bundle of accesory functions.

## Analysis

The `analysis` folder contains the following subfolders:

- `data`: Results from the models and csv with the metadata from fixations.
- `matlab`: Main matlab scripts to perform the analysis.
- `python`: Scripts for extrating predictors and filtering conditions to analyze.
- `rawcode`: 
## Usage

To use this code, follow these steps:

1. Clone the repository:

   ```bash
   git clone https://github.com/your-username/repository-name.git
