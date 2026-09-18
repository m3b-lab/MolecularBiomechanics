# Bash exercises on CLN025

These exercises work on `2RVD.pdb`, the NMR structure of CLN025. The entry contains 20 models of the same molecule, one after the other; each model starts with a `MODEL` line and ends with `ENDMDL`. You can treat the models like the frames of a short trajectory and look at them in VMD:

```bash
vmd 2RVD.pdb
```

For each exercise, write a bash script that does the job. Work in one folder, because each exercise uses the output of the previous one. Start by copying the structure:

```bash
mkdir -p ~/molbiomech/bash-exercises && cd ~/molbiomech/bash-exercises
cp ~/molbiomech/labs/common/structures/2RVD.pdb .
```

Each `exerciseN/` folder has a solution script and the expected output, so you can check your work. Try on your own first.

1. **Coordinates only.** From the PDB file, make a file with only the lines holding atom coordinates, dropping all comments. Keep the models separated by their `MODEL` and `ENDMDL` lines. This file is the input of the next exercises. (`grep`)
   Check: the result has 3,260 `ATOM` lines and 20 `MODEL` lines.
2. **Fields.** Find how many text fields an `ATOM` line has, then extract the number and name of every atom into a text file. (`grep`, `awk`)
3. **PDB to XYZ.** Convert the file to the XYZ format, one frame per model:
   ```
   number of atoms
   title
   name x y z
   ...
   ```
   Open the result in VMD to check it.
4. **Backbone.** Extract the backbone atoms (N, CA, C) of every frame and save them as an XYZ trajectory.
5. **Average structure.** Compute the average position of every backbone atom over the 20 frames and write it as an XYZ file. Compare it with the trajectory in VMD. The 20 models are already superimposed on each other, which is what makes a plain average meaningful here; for an MD trajectory you would first fit every frame onto a reference (LAB 4).
6. **RMSD.** For every frame, compute the RMSD of the CA atoms from the average structure of exercise 5:

   RMSD = sqrt( Σ (squared distance of each CA from its average position) / number of CA atoms )

   The 20 NMR models differ from their average by less than 1 Å: the ensemble is tight.
