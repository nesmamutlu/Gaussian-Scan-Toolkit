#!/usr/bin/env python3

import argparse
import re
import os
from ase import Atoms
from ase.io import write

# Conversion factor from Hartree (Gaussian unit) to electron-volts (eV)
HARTREE_TO_EV = 27.211386245988

def parse_gaussian_log(logfile_path):
    """
    Parses a Gaussian .log file to:
    - Extract the most recent SCF Done energy.
    - Extract the most recent 'Standard orientation' geometry.
    - Commit the geometry and energy whenever '-- Stationary point found.' is encountered.
    """
    
    with open(logfile_path, 'r') as f:
        log_lines = f.readlines()

    # Regex to find the SCF Done energy line
    energy_pattern = re.compile(r"SCF Done:\s+E\([^)]+\)\s+=\s+([-]?\d+\.\d+)")
    
    # Mapping of atomic numbers to element symbols
    atom_symbols = {
        1: "H", 6: "C", 7: "N", 8: "O",
        9: "F", 15: "P", 16: "S", 17: "Cl"
    }

    atoms_list = []
    last_energy = None
    last_positions = None
    last_symbols = None

    capture_geom = False
    skip = 0
    tmp_positions = []
    tmp_symbols = []
    
    # Get the base filename (without extension) for unique labeling
    base_name = os.path.splitext(os.path.basename(logfile_path))[0]

    for line in log_lines:

        # ---- ENERGY EXTRACTION ----
        m = energy_pattern.search(line)
        if m:
            last_energy = float(m.group(1))

        # ---- GEOMETRY EXTRACTION (Standard Orientation) ----
        if "Standard orientation:" in line:
            capture_geom = True
            skip = 4 # Skip the table headers in Gaussian output
            tmp_positions = []
            tmp_symbols = []
            continue

        if capture_geom:
            if skip > 0:
                skip -= 1
                continue

            if "-----" in line: # End of the geometry table
                last_positions = tmp_positions
                last_symbols = tmp_symbols
                capture_geom = False
            else:
                parts = line.split()
                Z = int(parts[1]) # Atomic number
                x, y, z = map(float, parts[3:6]) # Coordinates
                tmp_symbols.append(atom_symbols.get(Z, "X"))
                tmp_positions.append((x, y, z))

        # ---- COMMIT STEP (Only when a stationary point is confirmed) ----
        if "-- Stationary point found." in line:
            if last_positions is None or last_energy is None:
                continue

            atoms = Atoms(
                symbols=last_symbols,
                positions=last_positions
            )

            # Assign a unique label using the step index to keep the header organized
            step_idx = len(atoms_list)
            atoms.info["label"] = f"{base_name}_step{step_idx}"
            
            # Store energy and step metadata for Machine Learning training
            atoms.info["energy"] = last_energy
            atoms.info["energy_ev"] = last_energy * HARTREE_TO_EV
            atoms.info["scan_step"] = step_idx

            atoms_list.append(atoms)

            # Reset variables for the next optimization step
            last_positions = None
            last_symbols = None
            last_energy = None

    return atoms_list, base_name


def main():
    parser = argparse.ArgumentParser(
        description="Gaussian Log -> EXTXYZ (Auto-naming & Unique Step Labeling)"
    )
    parser.add_argument("logfile", help="Gaussian .log file path")
    parser.add_argument(
        "-o", "--output",
        help="Custom output filename (Default: log_name.extxyz)"
    )

    args = parser.parse_args()

    if not os.path.exists(args.logfile):
        print(f"Error: {args.logfile} not found.")
        return

    atoms_list, base_name = parse_gaussian_log(args.logfile)

    # Set default output name if none is provided via -o
    output_filename = args.output if args.output else f"{base_name}.extxyz"

    if atoms_list:
        write(output_filename, atoms_list, format="extxyz")
        print(f"✔ Success: {len(atoms_list)} structures extracted.")
        print(f"✔ Output File: {output_filename}")
        print(f"✔ Label Format: {base_name}_stepX")
    else:
        print("✘ Failed: No stationary points detected in the log file.")


if __name__ == "__main__":
    main()
