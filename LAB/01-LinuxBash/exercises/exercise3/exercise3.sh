#!/bin/bash
# Exercise 3: convert the models of a PDB file into an XYZ trajectory.
# Usage: ./exercise3.sh out-2RVD.pdb    (writes out-2RVD.xyz)
#
# An XYZ file has one block per frame:
#   number of atoms
#   a title line
#   name x y z          (one line per atom)
# Each MODEL of the PDB file becomes one frame.

input=$1
output="$(basename "$input" .pdb).xyz"

# number of atoms in the first model
n_atoms=$(awk '/^MODEL/ {m++} m == 1 && /^ATOM/ {n++} END {print n}' "$input")

# at every MODEL line start a new frame, then print name and x y z of each atom
awk -v n="$n_atoms" '/^MODEL/ {print n; print "frame", $2} /^ATOM/ {print $3, $7, $8, $9}' "$input" > "$output"
