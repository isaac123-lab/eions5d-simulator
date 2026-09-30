# EIONS 5D Simulator (`simulateur_eions5d_v3.py`)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

## Description
This repository contains the numerical implementation of the **EIONS 5D** theoretical framework (*Electron / Induced Electrodynamics from Non-linear Solitonic Wave*). 

The simulator integrates the stochastic 5D dynamics using a **Runge-Kutta 4th Order - Maruyama (RK4-Maruyama)** scheme, computes the system's $5 \times 5$ Jacobian matrix, and performs cosmological parameter fitting against Planck, Pantheon+, and BAO datasets.

## Key Physical Results
- **Hubble Tension Resolution:** $H_0 = 71.2 \pm 0.8 \text{ km/s/Mpc}$
- **Cosmological Fit Quality:** $\chi^2 / N_{\text{dof}} = 0.98$ (vs $1.02$ for $\Lambda\text{CDM}$)
- **Early Structure Formation:** Accelerated collapse $\delta_{5D}(a)$ for high-redshift galaxies ($z > 8$) matching JWST observations.

## Installation & Usage

```bash
git clone https://github.com/eions5d/eions5d-simulator.git
cd eions5d-simulator
pip install numpy scipy matplotlib astropy
To run the stochastic trajectory integration:Bashpython simulateur_eions5d_v3.py --mode trajectory --dt 0.001 --steps 10000
To run the cosmological parameter fitting ($\chi^2$ computation):Bashpython simulateur_eions5d_v3.py --mode fit_cosmo
CitationIf you use this simulator in your research, please cite:Extrait de code@article{eions5d2026,
  title={EIONS 5D Theory: From Quantum Solitons to Warped Scale Cosmology},
  author={Davila, Isaac and EIONS 5D Collaboration},
  year={2026}
}





