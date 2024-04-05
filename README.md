# Comparative Analysis of SARS-CoV-2 Neutralization Titers across Assays and Human and Animal Sera

This repository contains the raw data and code to reproduce the analyses presented in `Comparative Analysis of SARS-CoV-2 Neutralization Titers across Assays and Human and Animal Sera`.

## Directory structure

`main_figures`: Subdirectories for each main figure containing the code to produce each figure, as well as the main figures themselves.

`som_figures`: Subdirectories for SOM figures containing the code to produce each figure, as well as the SOM figures themselves.

`data`: Data and code to produce the data.

`code`: Code and helper functions.


## Dependencies

All dependencies are listed in the `DEPENDENCIES` file. Apart from the dependencies available through cran, the following packages specific to antigenic cartography and titer analysis are required:

##### Racmacs
Available from [github.com/shwilks/Racmacs](github.com/shwilks/Racmacs).

##### titertools
Available from [github.com/shwilks/titertools](github.com/shwilks/titertools).

## Generating data and figures

To re-generate data associated with figures, run `code/generate_data.R`.

To re-generate the figures, run `code/generate_figures.R`.