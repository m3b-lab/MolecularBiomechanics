"""The course's fixed definitions for CLN025, shared by LABs 4, 5 and 6.

Every lab measures the state of the hairpin with the same functions, so that
results can be compared between labs, runs and force fields:

- rmsd_nmr: backbone (N, CA, C) RMSD of residues 1-10 from NMR model 1 of 2RVD
- q_value: fraction of native contacts (Best, Hummer & Eaton, PNAS 2013,
  doi:10.1073/pnas.1311599110), native contacts taken from NMR model 1
- end_to_end: distance between the CA atoms of Tyr1 and Tyr10
- folded: Q >= 0.7 and RMSD < 0.15 nm

Only residues 1-10 are used, so the caps (ACE, NH2) and the force field's naming
of hydrogens do not matter: the same code works for amber99sb-ildn and gromos54a7.

Use it from a notebook in a lab folder with:
    import sys; sys.path.append("../common/scripts")
    import cln025
"""
import os

import mdtraj as md
import numpy as np

KT = 0.008314 * 300          # kJ/mol at 300 K
Q_FOLDED = 0.7               # folded: Q at least this...
RMSD_FOLDED = 0.15           # ...and backbone RMSD from NMR model 1 below this (nm)

_HERE = os.path.dirname(os.path.abspath(__file__))
NMR = md.load(os.path.join(_HERE, "..", "structures", "cln025_capped.pdb"))


def xvg(path):
    """Columns of a GROMACS .xvg file as a 2D array."""
    return np.loadtxt(path, comments=["#", "@"])


def atoms(top, names):
    """Indices of the atoms of residues 1-10 with the given names, ordered by residue and name."""
    found = [(a.residue.resSeq, a.name, a.index) for a in top.atoms
             if 1 <= a.residue.resSeq <= 10 and a.name in names]
    return [i for _, _, i in sorted(found)]


def rmsd_nmr(traj):
    """Backbone RMSD (nm) of residues 1-10 from NMR model 1, after superposition."""
    backbone = ("N", "CA", "C")
    return md.rmsd(traj, NMR, atom_indices=atoms(traj.topology, backbone),
                   ref_atom_indices=atoms(NMR.topology, backbone))


def native_contacts():
    """Heavy-atom pairs of NMR model 1, in residues 1-10 more than 3 apart, closer than 0.45 nm.

    Returns a list of ((resSeq, name), (resSeq, name), native distance in nm)."""
    heavy = [a for a in NMR.topology.atoms if 1 <= a.residue.resSeq <= 10 and a.element.symbol != "H"]
    pairs = [(a, b) for i, a in enumerate(heavy) for b in heavy[i + 1:]
             if abs(a.residue.resSeq - b.residue.resSeq) > 3]
    r0 = md.compute_distances(NMR, [(a.index, b.index) for a, b in pairs])[0]
    return [((a.residue.resSeq, a.name), (b.residue.resSeq, b.name), d)
            for (a, b), d in zip(pairs, r0) if d < 0.45]


def q_value(traj, beta=50.0, lam=1.8):
    """Fraction of native contacts, frame by frame (beta in 1/nm, lambda dimensionless)."""
    contacts = native_contacts()
    index = {(a.residue.resSeq, a.name): a.index for a in traj.topology.atoms}
    pairs = [(index[a], index[b]) for a, b, _ in contacts]
    r0 = np.array([d for _, _, d in contacts])
    r = md.compute_distances(traj, pairs)
    x = np.clip(beta * (r - lam * r0), -50, 50)      # beyond ±50 a contact is fully formed or fully broken
    return np.mean(1.0 / (1.0 + np.exp(x)), axis=1)


def end_to_end(traj):
    """Distance (nm) between the CA atoms of Tyr1 and Tyr10, frame by frame."""
    ca = {a.residue.resSeq: a.index for a in traj.topology.atoms if a.name == "CA"}
    return md.compute_distances(traj, [(ca[1], ca[10])])[:, 0]


def observables(traj):
    """End-to-end distance, RMSD from NMR model 1 and Q, frame by frame."""
    return {"e2e": end_to_end(traj), "rmsd": rmsd_nmr(traj), "q": q_value(traj)}


def folded(obs):
    """True for the frames that count as folded."""
    return (obs["q"] >= Q_FOLDED) & (obs["rmsd"] < RMSD_FOLDED)


def read_xpm(path):
    """Values of a GROMACS .xpm file: bin edges x and y, and G[y, x] in the legend's units."""
    text = open(path).read().splitlines()
    header = next(i for i, line in enumerate(text) if line.startswith('"') and len(line.split()) >= 4)
    ncol, nrow, ncolors, cpp = map(int, text[header].strip('",').split())
    values, xs, ys, rows = {}, [], [], []
    for line in text[header + 1:header + 1 + ncolors]:
        values[line[1:1 + cpp]] = float(line.split('/* "')[1].split('"')[0])
    for line in text[header + 1 + ncolors:]:
        if line.startswith("/* x-axis:"):
            xs += [float(v) for v in line[10:].replace("*/", "").split()]
        elif line.startswith("/* y-axis:"):
            ys += [float(v) for v in line[10:].replace("*/", "").split()]
        elif line.startswith('"'):
            row = line.strip().strip('",')
            rows.append([values[row[i:i + cpp]] for i in range(0, len(row), cpp)])
    return np.array(xs), np.array(ys), np.array(rows[::-1])     # the file lists the top row first
