#!/usr/bin/env python3
"""Cap one model of PDB 2RVD as ACE-YYDPETGTWY-NH2.

Takes the heavy atoms of the chosen NMR model, adds an acetyl group (ACE,
residue 0) on the N-terminus of Tyr1 and an amide nitrogen (NH2, residue 11)
on the C-terminus of Tyr10, and writes a PDB file for gmx pdb2gmx (run it
with -ignh). The caps are built with standard peptide geometry; energy
minimisation removes any small strain.

Usage:
    python build_caps.py 2RVD.pdb MODEL OUTPUT.pdb [--gromos]

--gromos names the acetyl methyl CA (GROMOS naming) instead of CH3 (AMBER).
"""
import argparse

import numpy as np


def read_model(path, model):
    """Heavy atoms of one MODEL of a multi-model PDB file."""
    atoms, inside = [], False
    for line in open(path):
        if line.startswith("MODEL"):
            inside = int(line.split()[1]) == model
        elif line.startswith("ENDMDL") and inside:
            break
        elif inside and line.startswith("ATOM") and line[76:78].strip() != "H":
            atoms.append((line[12:16].strip(), line[17:20], int(line[22:26]),
                          np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    if not atoms:
        raise SystemExit(f"model {model} not found in {path}")
    return atoms


def place(a, b, c, bond, angle, torsion):
    """Point d with |cd| = bond, angle bcd = angle and dihedral abcd = torsion (degrees)."""
    angle, torsion = np.radians(angle), np.radians(torsion)
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    return c + bond * (-np.cos(angle) * bc + np.sin(angle) * (np.cos(torsion) * m + np.sin(torsion) * n))


def dihedral(a, b, c, d):
    b1 = (c - b) / np.linalg.norm(c - b)
    v = (a - b) - np.dot(a - b, b1) * b1
    w = (d - c) - np.dot(d - c, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("pdb")
    parser.add_argument("model", type=int)
    parser.add_argument("output")
    parser.add_argument("--gromos", action="store_true", help="GROMOS atom names for the acetyl cap")
    args = parser.parse_args()

    atoms = read_model(args.pdb, args.model)
    xyz = {(resi, name): r for name, _, resi, r in atoms}
    N1, CA1, C1 = xyz[(1, "N")], xyz[(1, "CA")], xyz[(1, "C")]
    N10, CA10, C10, O10 = xyz[(10, "N")], xyz[(10, "CA")], xyz[(10, "C")], xyz[(10, "O")]

    # Acetyl: carbonyl C bonded to Tyr1 N, trans peptide bond, Tyr1 phi = -120 (beta strand)
    ace_c = place(C1, CA1, N1, 1.335, 121.7, -120.0)
    ace_o = place(CA1, N1, ace_c, 1.229, 122.7, 0.0)
    ace_me = place(CA1, N1, ace_c, 1.522, 116.2, 180.0)
    # Amide N bonded to Tyr10 C, in the carbonyl plane, opposite O
    nh2_n = place(N10, CA10, C10, 1.335, 116.2, dihedral(N10, CA10, C10, O10) + 180.0)

    methyl = "CA" if args.gromos else "CH3"
    records = [(methyl, "ACE", 0, ace_me), ("C", "ACE", 0, ace_c), ("O", "ACE", 0, ace_o)]
    records += atoms
    records += [("N", "NH2", 11, nh2_n)]

    with open(args.output, "w") as out:
        out.write(f"TITLE     CLN025 (ACE-YYDPETGTWY-NH2), capped from PDB 2RVD model {args.model}\n")
        for serial, (name, resn, resi, r) in enumerate(records, 1):
            field = name if len(name) == 4 else " " + name
            record = "HETATM" if resn in ("ACE", "NH2") else "ATOM  "
            out.write(f"{record}{serial:5d} {field:<4s} {resn:3s} A{resi:4d}    "
                      f"{r[0]:8.3f}{r[1]:8.3f}{r[2]:8.3f}  1.00  0.00          {name[0]:>2s}\n")
        out.write("TER\nEND\n")
    print(f"wrote {args.output}: model {args.model}, {len(records)} heavy atoms")


if __name__ == "__main__":
    main()
