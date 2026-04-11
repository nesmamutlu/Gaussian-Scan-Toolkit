import argparse
import re
import os
from ase import Atoms
from ase.io import write

# Conversion factor from Hartree to electron-Volts
HARTREE_TO_EV = 27.211386245988

def parse_gaussian_log(log_lines):
    """
    Parses Gaussian OPT / SCAN log files:
    - Tracks the latest SCF energy.
    - Captures the most recent Standard or Input orientation.
    - When '-- Stationary point found.' is encountered, saves the geometry and energy.
    """
    energy_pattern = re.compile(r"SCF Done:\s+E\([^)]+\)\s+=\s+([-]?\d+\.\d+)")
    
    # Dictionary to map atomic numbers to symbols
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

    for line in log_lines:
        # ---- ENERGY EXTRACTION ----
        m_energy = energy_pattern.search(line)
        if m_energy:
            last_energy = float(m_energy.group(1))

        # ---- GEOMETRY EXTRACTION ----
        # Checks for both Standard and Input orientations
        if "Standard orientation:" in line or "Input orientation:" in line:
            capture_geom = True
            skip = 4
            tmp_positions = []
            tmp_symbols = []
            continue

        if capture_geom:
            if skip > 0:
                skip -= 1
                continue
            if "-----" in line:
                if tmp_positions:
                    last_positions = tmp_positions
                    last_symbols = tmp_symbols
                capture_geom = False
            else:
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        Z = int(parts[1])
                        x, y, z = map(float, parts[3:6])
                        tmp_symbols.append(atom_symbols.get(Z, "X"))
                        tmp_positions.append((x, y, z))
                    except: continue

        # ---- FINALIZE OPTIMIZED STEP ----
        if "-- Stationary point found." in line:
            if last_positions and last_energy is not None:
                atoms = Atoms(symbols=last_symbols, positions=last_positions)
                atoms.info["energy_hartree"] = last_energy
                atoms.info["energy_ev"] = last_energy * HARTREE_TO_EV
                atoms_list.append(atoms)
                
                # Reset energy for the next step, but keep positions 
                # in case the next step doesn't print a new orientation immediately
                last_energy = None

    return atoms_list

def main():
    parser = argparse.ArgumentParser(
        description="Extract optimized geometries from Gaussian SCAN/OPT logs to separate PDB files."
    )
    parser.add_argument("logfile", help="Path to the Gaussian .log file")
    args = parser.parse_args()

    if not os.path.exists(args.logfile):
        print(f"✘ Error: {args.logfile} not found.")
        return

    # Get the base filename without extension for naming outputs
    base_name = os.path.splitext(os.path.basename(args.logfile))[0]

    with open(args.logfile, "r") as f:
        lines = f.readlines()

    all_structures = parse_gaussian_log(lines)

    if not all_structures:
        print("✘ No optimized structures found. Check if '-- Stationary point found.' exists in the log.")
        return

    print(f"✔ Found {len(all_structures)} structures. Saving with prefix '{base_name}'...")

    # Save each structure as logfilename_n.pdb
    for i, atoms in enumerate(all_structures, start=1):
        filename = f"{base_name}_{i}.pdb"
        
        try:
            # Attempt to write in Protein Data Bank format
            write(filename, atoms, format="proteindatabank")
            print(f"✔ Saved: {filename}")
        except Exception as e:
            # Fallback to XYZ if PDB writing fails
            xyz_name = filename.replace(".pdb", ".xyz")
            write(xyz_name, atoms)
            print(f"⚠ PDB Error! Saved as {xyz_name} instead.")

if __name__ == "__main__":
    main()
