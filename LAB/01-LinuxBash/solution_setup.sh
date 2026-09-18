#!/bin/bash
# Solution of the setup.sh exercise (LAB 1 part 2, section 3).
# Usage: bash solution_setup.sh MODEL          (MODEL = 1 ... 20)
# Prepares a working folder for the LAB 3 simulation of one 2RVD model.
# LAB 3 adds the GROMACS commands.

# Where the shared course files are
COMMON=${COMMON:-$HOME/molbiomech/labs/common}

# Accept exactly one argument, a whole number from 1 to 20
if [ $# -ne 1 ] || ! [[ $1 =~ ^[0-9]+$ ]] || [ "$1" -lt 1 ] || [ "$1" -gt 20 ]; then
    echo "Usage: $0 MODEL   (MODEL is a number from 1 to 20)"
    exit 1
fi
model=$1
folder=cln025_model$model

# Working folder with one subfolder per stage
mkdir -p "$folder"/{00-build,01-em,02-nvt,03-npt,04-md}
echo "Created $folder with 00-build 01-em 02-nvt 03-npt 04-md"

# The course parameter files
cp "$COMMON"/mdp/*.mdp "$folder"/
echo "Copied $(cd "$folder" && ls *.mdp | tr '\n' ' ')"

# The chosen model, capped as ACE-YYDPETGTWY-NH2
python3 "$COMMON"/scripts/build_caps.py "$COMMON"/structures/2RVD.pdb "$model" "$folder"/00-build/cln025_capped.pdb > /dev/null
echo "Capped model $model into $folder/00-build/cln025_capped.pdb"
