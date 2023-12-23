# Interval-valued Data Analysis

## Interval-Valued Principal Component Analysis and Anomaly Detection

This repository holds code that implements the method proposed by Lin et al. (2022) (https://doi.org/10.1007/s11634-022-00527-1).

It also holds the code for the appplied work "Air Quality Data Analysis with Symbolic Principal Components" which was presented at the XXVI Congresso da Sociedade Portuguesa de Estatística (https://w3.math.uminho.pt/SPE2023/) and at the Symbolic Data Analysis Workshop 2023 (https://sda2018.wixsite.com/sda2023paris).

## Code

 - `setup.r` contains the implementation of several functions that compose the method.
 - `example.r` contains a simple example.
 - `air_quality.r` contains the application of the method on an air quality dataset.

## Data

Files containing the air quality dataset obtained from a monitoring station in Entrecampos, Lisbon. Each file corresponds to 9 pollutants' concentration measures during the years 2019, 2020, and 2021. These files were retrieved from the Portuguese Environment Agency database available at https://qualar.apambiente.pt/.
