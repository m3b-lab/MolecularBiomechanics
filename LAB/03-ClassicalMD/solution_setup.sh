#!/bin/bash
# Solution of the setup.sh exercise (LAB 3, section 8): the LAB 1 script
# completed with the GROMACS commands of LAB 3.
# Usage: bash solution_setup.sh MODEL [LENGTH_NS]
#   MODEL      2RVD model to simulate, a number from 1 to 20
#   LENGTH_NS  length of the 300 K production run in ns (default 10)
# Builds capped CLN025 in water (amber99sb-ildn, TIP3P, 0.15 M NaCl) and runs
# energy minimisation, NVT and NPT equilibration and the production run.
# It takes hours: start it in the background from a terminal, e.g.
#   nohup bash solution_setup.sh 7 > setup_model7.out 2>&1 &

set -e      # stop at the first command that fails

# Where the shared course files are
COMMON=${COMMON:-$HOME/molbiomech/labs/common}

# Accept a model number from 1 to 20 and an optional run length
if [ $# -lt 1 ] || [ $# -gt 2 ] || ! [[ $1 =~ ^[0-9]+$ ]] || [ "$1" -lt 1 ] || [ "$1" -gt 20 ]; then
    echo "Usage: $0 MODEL [LENGTH_NS]   (MODEL is a number from 1 to 20)"
    exit 1
fi
model=$1
length_ns=${2:-10}
folder=cln025_model$model

# LAB 1: working folder, parameter files, capped model
mkdir -p "$folder"/{00-build,01-em,02-nvt,03-npt,04-md}
cp "$COMMON"/mdp/*.mdp "$folder"/
python3 "$COMMON"/scripts/build_caps.py "$COMMON"/structures/2RVD.pdb "$model" "$folder"/00-build/cln025_capped.pdb > /dev/null
echo "Capped model $model into $folder/00-build/cln025_capped.pdb"

# Topology, box, water and ions
cd "$folder"/00-build
gmx pdb2gmx -f cln025_capped.pdb -o cln025.gro -p topol.top -i posre.itp \
            -ff amber99sb-ildn -water tip3p -ignh > pdb2gmx.out 2>&1
gmx editconf -f cln025.gro -o box.gro -c -d 1.0 -bt dodecahedron > editconf.out 2>&1
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top > solvate.out 2>&1
gmx grompp -f ../em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1 > grompp_ions.out 2>&1
echo SOL | gmx genion -s ions.tpr -o system.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 > genion.out 2>&1
echo "Built the system: $(sed -n 2p system.gro | tr -d ' ') atoms"

# Energy minimisation
cd ../01-em
gmx grompp -f ../em.mdp -c ../00-build/system.gro -p ../00-build/topol.top -o em.tpr > grompp.out 2>&1
gmx mdrun -deffnm em -ntmpi 1 > mdrun.out 2>&1
echo "Minimised: $(grep 'Maximum force' em.log)"

# NVT equilibration, protein restrained
cd ../02-nvt
gmx grompp -f ../nvt.mdp -c ../01-em/em.gro -r ../01-em/em.gro -p ../00-build/topol.top -o nvt.tpr > grompp.out 2>&1
gmx mdrun -deffnm nvt -ntmpi 1 > mdrun.out 2>&1
echo "NVT equilibration done"

# NPT equilibration, protein restrained
cd ../03-npt
gmx grompp -f ../npt.mdp -c ../02-nvt/nvt.gro -r ../02-nvt/nvt.gro -t ../02-nvt/nvt.cpt -p ../00-build/topol.top -o npt.tpr > grompp.out 2>&1
gmx mdrun -deffnm npt -ntmpi 1 > mdrun.out 2>&1
echo "NPT equilibration done"

# Production at 300 K and 1 bar
cd ../04-md
nsteps=$(python3 -c "print(round($length_ns * 1000 / 0.002))")
sed "s/^nsteps .*/nsteps          = $nsteps/" ../md.mdp > md_300K.mdp
gmx grompp -f md_300K.mdp -c ../03-npt/npt.gro -t ../03-npt/npt.cpt -p ../00-build/topol.top -o md_300K.tpr > grompp.out 2>&1
echo "Production: $length_ns ns at 300 K, started $(date)"
gmx mdrun -deffnm md_300K -ntmpi 1 > mdrun.out 2>&1
echo "Production done $(date): $(grep 'Performance:' md_300K.log)"
