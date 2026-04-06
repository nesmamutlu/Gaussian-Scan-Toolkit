# Gaussian-Scan-Toolkit
Automated Gaussian 1D Bond Scan Generator
This script automates the creation of Gaussian input files (.com) for 1D relaxed scans by identifying chemical bonds within an optimized structure. 
It is designed to streamline the transition from a completed geometry optimization to a potential energy surface (PES) scan.

Intelligent Bond Detection: Uses Covalent Radii logic to automatically detect valid chemical bonds between atoms (C, H, O, N, P, S, F, Cl) based on their spatial distance.

Interactive Selection: Lists all detected bonded pairs and allows the user to select a specific bond index for the scan.

ModRedundant Automation: Automatically generates the Opt=ModRedundant syntax and the required scan command (B atom_i atom_j S steps step_size).

HPC Optimized: Pre-configured with high-performance computing settings, including large memory allocation (224GB) and shared processor usage (112 cores).

Usage
Run the script by providing a Gaussian .log file from a previous optimization:
python auto_gaussian_1d_scan_from_log_bonded.py optimized.log

Workflow Details
Read Geometry: Extracts the final optimized "Standard Orientation" from the Gaussian log.

Identify Bonds: Calculates pairwise distances and compares them against a scaled covalent radii sum (default scale = 1.2).

Input Generation: Prompts for scan parameters (steps, step size, and QM method) to produce a ready-to-submit .com file.

## 📊 Gaussian Log to EXTXYZ Converter

This script is a post-processing tool designed to extract optimized geometries and their corresponding energies from Gaussian 1D scan or optimization log files. 
It outputs the data in the **Extended XYZ (.extxyz)** format, which is the standard for training Machine Learning Interatomic Potentials (MLIPs).

### Key Features
* **Energy Extraction**: Automatically locates `SCF Done` values for each optimized step and records them in Hartree.
* **Unit Conversion**: Automatically converts Hartree energies to **electron-volts (eV)** for standard ML compatibility.
* **Stationary Point Detection**: Specifically searches for `-- Stationary point found.` to ensure only fully converged geometries are extracted.
* **ASE Integration**: Uses the **Atomic Simulation Environment (ASE)** library to handle molecular structures and metadata.
* **Metadata Mapping**: Each frame in the output file contains the total energy, atomic symbols, and 3D positions.

### Usage
```bash
python extract_g16_scan_coordinates_to_extxyz.py optimized_scan.log -o output_data.extxyz
