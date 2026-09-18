#!/bin/bash
# Solution of exercise 2.2.3 (LAB 1 part 2).
# Usage: bash solution_modelstats.sh 2RVD.pdb
# Splits a multi-model PDB file into one file per model, then writes the
# Tyr1 CA - Tyr10 CA distance of every model to stat.csv and the model with
# the largest distance to maxdist.stat.

# Check that exactly one argument was given and that it is a file
if [ $# -ne 1 ]; then
    echo "Usage: $0 PDBFILE"
    exit 1
fi
input=$1
if [ ! -f "$input" ]; then
    echo "File not found: $input"
    exit 1
fi

# One file per model: at every MODEL line, switch the output file
mkdir -p models
awk '/^MODEL/ {file = sprintf("models/model_%02d.pdb", $2)} /^ATOM/ {print > file}' "$input"

# Header line, then one row per model
echo "#file,distance_A" > stat.csv
for model in models/model_*.pdb; do
    # $6 is the residue number, $7-$9 are x y z (in Angstrom)
    distance=$(awk '$3 == "CA" && $6 == 1  {x1 = $7; y1 = $8; z1 = $9}
                    $3 == "CA" && $6 == 10 {x2 = $7; y2 = $8; z2 = $9}
                    END {printf "%.2f", sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)}' "$model")
    echo "$(basename "$model"),$distance" >> stat.csv
done

# The largest distance: sort the rows (without the header) by the 2nd comma-separated field
best=$(grep -v '^#' stat.csv | sort -t, -k2 -g | tail -n 1)
echo "${best%,*} has the largest end-to-end distance (${best#*,} A)." > maxdist.stat
cat maxdist.stat
