#!/usr/bin/env python3
import numpy as np
import sys

# ============================================================
# 1. READ OPTIMIZED GEOMETRY FROM GAUSSIAN LOG
# ============================================================

def read_geometry_from_log(logfile):
    with open(logfile, "r") as f:
        lines = f.readlines()

    geometry = []
    capture = False
    skip = 0

    for line in lines:
        if "Standard orientation:" in line:
            geometry = []
            capture = True
            skip = 4
            continue

        if capture:
            if skip > 0:
                skip -= 1
                continue

            if "-----" in line:
                break

            parts = line.split()
            Z = int(parts[1])
            x, y, z = map(float, parts[3:6])
            geometry.append((Z, x, y, z))

    if not geometry:
        raise RuntimeError("No optimized geometry found in log file")

    return geometry


# ============================================================
# 2. COVALENT RADII (Å)
# ============================================================

COVALENT_RADII = {
    1: 0.31,  # H
    6: 0.76,  # C
    7: 0.71,  # N
    8: 0.66,  # O
    9: 0.57,  # F
    15: 1.07, # P
    16: 1.05, # S
    17: 1.02  # Cl
}

Z2SYM = {
    1: "H", 6: "C", 7: "N", 8: "O",
    9: "F", 15: "P", 16: "S", 17: "Cl"
}


# ============================================================
# 3. COMPUTE BONDED DISTANCES ONLY
# ============================================================

def compute_bonded_distances(geometry, scale=1.2):
    coords = np.array([[x, y, z] for _, x, y, z in geometry])
    Zs = [Z for Z, _, _, _ in geometry]

    bonded = []
    n = len(coords)

    for i in range(n):
        for j in range(i + 1, n):
            Zi, Zj = Zs[i], Zs[j]
            if Zi not in COVALENT_RADII or Zj not in COVALENT_RADII:
                continue

            rij = np.linalg.norm(coords[i] - coords[j])
            cutoff = scale * (COVALENT_RADII[Zi] + COVALENT_RADII[Zj])

            if rij <= cutoff:
                bonded.append((i + 1, j + 1, Zi, Zj, rij))

    return bonded


# ============================================================
# 4. WRITE GAUSSIAN RELAXED SCAN INPUT
# ============================================================

def write_gaussian_scan(
    filename,
    geometry,
    atom_i,
    atom_j,
    steps,
    step_size,
    method,
    charge,
    multiplicity
):
    with open(filename, "w") as f:
        f.write("%chk=scan.chk\n")
        f.write("%mem=224GB\n")
        f.write("%nprocshared=112\n")
        f.write(f"#p {method} Opt=(ModRedundant,loose)\n\n")
        f.write("Automated bonded 1D relaxed scan\n\n")
        f.write(f"{charge} {multiplicity}\n")

        for Z, x, y, z in geometry:
            f.write(f"{Z2SYM.get(Z,'X'):2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")

        f.write("\n")
        f.write(f"B {atom_i} {atom_j} S {steps} {step_size:.4f}\n\n\n")


# ============================================================
# 5. MAIN
# ============================================================

def main():
    if len(sys.argv) != 2:
        print("Usage: python auto_gaussian_1d_scan_from_log_bonded.py optimized.log")
        sys.exit(1)

    logfile = sys.argv[1]

    geometry = read_geometry_from_log(logfile)
    bonded = compute_bonded_distances(geometry)

    if not bonded:
        raise RuntimeError("No bonded atom pairs detected")

    print("\nBonded distances available for scan:\n")
    for idx, (i, j, Zi, Zj, d) in enumerate(bonded):
        print(f"[{idx}] {Z2SYM[Zi]}{i} – {Z2SYM[Zj]}{j} : {d:.3f} Å")

    choice = int(input("\nSelect bond index to scan: "))
    atom_i, atom_j, Zi, Zj, d0 = bonded[choice]

    steps = int(input("Number of scan steps: "))
    step_size = float(input("Step size (Å): "))
    method = input("QM method [wB97XD/6-311++G(3df,3pd)]: ") \
             or "wB97XD/6-311++G(3df,3pd)"

    outname = f"scan_{Z2SYM[Zi]}{atom_i}_{Z2SYM[Zj]}{atom_j}.com"

    write_gaussian_scan(
        outname,
        geometry,
        atom_i,
        atom_j,
        steps,
        step_size,
        method,
        charge=0,
        multiplicity=1
    )

    print(f"\nGaussian scan input written to: {outname}")


if __name__ == "__main__":
    main()
