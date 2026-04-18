# SdS-EMRI-dynamics
Contains infrastructure in Mathematica generating GW waveforms, and orbital trajectories of EMRIs under radiation reaction in SdS spacetime

# Extreme-Mass Ratio Inspirals in Schwarzschild-de Sitter Spacetimes (Weak-field)

**Author:** John Adrian Villanueva (University of the Philippines Diliman)  
**Associated Paper:** [arXiv:2602.17154](https://arxiv.org/abs/2602.17154) (Submitted to PRD)

## Overview
This repository contains the Mathematica infrastructure used to compute weak-field orbital dynamics and extract gravitational waveforms for extreme-mass ratio inspirals (EMRIs) in Schwarzschild-de Sitter backgrounds. 

By modeling the effects of the cosmological constant on EMRI trajectories, this code provides a framework for analyzing environmental perturbations relevant to the upcoming LISA mission.

## Key Features
* **Metric Formulation:** Exact analytical implementation of SdS metric perturbations.
* **Trajectory Integration:** Highly optimized `NDSolve` routines capable of tracking orbital phase evolution.
* **Waveform Extraction:** Modules to compute and visualize the resulting gravitational wave strain.

## Repository Structure
* `/notebooks`: The core Mathematica (`.nb`) files used for the derivations and numerical integration.
* `/readable_code`: **PDF exports of all notebooks.** (Recommended for quick browser viewing of equations, code, and plots without needing Mathematica).
* `/plots`: High-resolution figures generated for the manuscript.

## How to Run
1. Clone this repository to your local machine.
2. Open the desired notebook in `/notebooks` using Mathematica 12.0 or later.
3. The code is structured to be self-contained. Simply evaluate the initialization cells to load the spacetime parameters and run the integration sequence.
