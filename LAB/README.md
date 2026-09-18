# Overview

These are the practical labs of the Molecular Biomechanics course. Every lab from LAB 1 onward works on the same molecule: **CLN025**, a 10-residue variant of chignolin (sequence YYDPETGTWY) that folds into a β-hairpin. It is the smallest system that behaves like a protein rather than a peptide, which makes it a good place to ask what holds a fold together and what it takes to pull it apart. The reference structure is **PDB 2RVD**.

All labs run inside the course virtual machine, which ships GROMACS, VMD, Python and the helper commands used below. Nothing has to be installed on your own machine apart from the virtualization software needed to import the VM (see LAB 0).

The labs build on each other:

```
LAB 0 → LAB 1 → LAB 2 → LAB 3 → LAB 4 ─┬→ LAB 5
                                       └→ LAB 6
```

LAB 6 closes the loop: its free-energy profile from pulling is compared against the equilibrium profile from LAB 4.

### LAB 0 - VM Setup and Benchmark

`00-VMSetup/`

Import the course VM into VirtualBox, give it a sensible share of your computer's resources, and check that everything works. You will also measure how fast GROMACS runs on your computer (**your own ns/day**): every later lab sizes its simulations against that number.

### LAB 1 - Linux, Bash and the Terminal

`01-LinuxBash/`

Working in the Linux shell: navigating the file system, handling files, `grep`/`awk`/`sed`, loops and scripts, monitoring and backgrounding long jobs with `top`, `htop` and `nohup`. All exercises use CLN025 structures and a small set of PDB entries. You will write the skeleton of a `setup.sh` script that LAB 3 completes, and extract statistics from PDB files with the command line.

### LAB 2 - VMD and the PDB Format

`02-VMD-PDB/`

The PDB format first (records, fixed columns, chains, occupancy and B-factor), then VMD: representations, coloring, the selection language, measuring distances, angles and dihedrals, and rendering. You will measure the interactions that stabilize the hairpin and produce a render of CLN025 that you will reuse for the rest of the course.

### LAB 3 - Classical MD

`03-ClassicalMD/`

The full GROMACS workflow on CLN025: topology, box and solvation, ions, energy minimization, NVT and NPT equilibration, production. Includes the box-shape comparison and a run matrix of folded and extended starting structures at two temperatures. The main force field is `amber99sb-ildn` with TIP3P water; `gromos54a7` with SPC water is run as a comparison, to see how much the result depends on the force field. Each student's replicas are pooled into a shared class dataset.

### LAB 4 - MD Analysis

`04-Analysis/`

RMSD, RMSF, radius of gyration, Ramachandran plots, end-to-end distance, hydrogen bonds and secondary structure, plus CLN025-specific observables such as the native contact fraction. Then principal component analysis, 2D free-energy surfaces, clustering of the conformational states, and the **equilibrium free-energy profile along the end-to-end distance** that LAB 6 is compared against.

### LAB 5 - Simulated Annealing

`05-SimulatedAnnealing/`

Heating and cooling schedules with the GROMACS `annealing` options, and a ladder of constant-temperature runs around the melting temperature of CLN025. You will build a melting curve, fit the melting temperature and the unfolding enthalpy, compare them with experiment, and look at hysteresis and at the unfolding pathway in VMD.

### LAB 6 - Steered MD

`06-SteeredMD/`

Unfolding CLN025 with the GROMACS pull code: force-extension curves, work, and free energies from many non-equilibrium pulls with Jarzynski's equality. Includes a sweep of pulling speeds to show how the result depends on the loading rate. The final step overlays the Jarzynski profile on the equilibrium profile from LAB 4: two independent routes to the same free energy.

## Reference data

Simulation outputs (trajectories, energies, logs) are not stored in this repository. Precomputed reference datasets are downloaded inside the VM with:

```bash
course-fetch <dataset>
course-fetch --all
```

This way every lab can be analysed on the same data, even if your own runs did not finish in time.

## Other folders

- `common/mdp/` holds the single copy of every GROMACS parameter file used in the labs.
- `Additional_Material/` holds optional material: `python-intro/` and `probability/` (the former introductory labs), `linux-bash/` (extra Linux and bash tutorials) and `amyloid-2BEG/` (the former steered-MD lab on the amyloid fibril 2BEG).

For a second walkthrough of a complete GROMACS run on a different protein, see the official [GROMACS introductory MD tutorial](https://gitlab.com/gromacs/online-tutorials/md-intro-tutorial).
