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
