#!/bin/bash
# Exercise 1: keep only the coordinate records of a PDB file.
# Usage: ./exercise1.sh 2RVD.pdb        (writes out-2RVD.pdb)
#
# grep -E accepts several patterns separated by |, and ^ anchors each
# pattern to the start of the line. We keep MODEL, ATOM, TER and ENDMDL,
# so the models stay separated, and drop everything else (HEADER, REMARK, ...).

input=$1
output="out-$(basename "$input")"

grep -E "^(MODEL|ATOM|TER|ENDMDL)" "$input" > "$output"
