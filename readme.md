# GBM-SIMUL
Basic computational framework for recreating glioblastoma progression in one dimension within microfluidic devices

This repository contains the main scripts for the reproduction of the results presented in the paper [Ayensa-Jimenez et al. (2020)](https://www.nature.com/articles/s41598-020-78215-3).

## Table of Contents
- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)
- [Contact](#contact)

## Description

The repository includes the following folders and files:

- `data`: This folder contains the experimental data:
    - `core.mat`: Pre-processing scripts.
    - `palisade.mat`: Scripts for the evaluation of new images.
    - `double_palisade.mat`: Post-processing scripts.
- `main.m`: Main script reproducing the simulation for the three experiments: necrotic core, pseudopalisade and double pseudopalisade.
- `simulationGBM.m`: Function for GBM progression simulation.

## Installation

1. Clone or download this repository to your local machine (you can use: `git clone https://github.com/Ayensa-Jimenez/GBM-SIMUL.git`).
2. Make sure you have Matlab R2020 or later installed.
3. Open Matlab and run the scripts-

## Usage 

### Reproducing the paper results
1. Run the script `main.m`

### Run the model with different model parameters
1. Modify the value of the parameters at the `main.m` script.
2. Run the `main.m` script.

## Contributing

If you want to contribute to this project, please follow these steps:

1. Fork the repository.
2. Create a new branch (`git checkout -b feature/new-feature`).
3. Make your changes and commit (`git commit -am 'Add a new feature'`).
4. Push your changes to your GitHub repository (`git push origin feature/new-feature`).
5. Create a new pull request.

## Contact

If you have any questions or suggestions, feel free to contact the author:

- Name: [Jacobo Ayensa Jiménez]
- Email: [jacoboaj@unizar.es]
