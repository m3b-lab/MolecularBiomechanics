# Shared course files

Files used by several labs. Every lab from LAB 3 on simulates the same molecule, CLN025, with the same protocol.

## `structures/`

| File | Content |
|---|---|
| `2RVD.pdb` | PDB entry 2RVD: 20 NMR models of CLN025 (sequence YYDPETGTWY) |
| `cln025_capped.pdb` | Model 1, capped as ACE-YYDPETGTWY-NH2, AMBER atom names |
| `cln025_capped_gromos.pdb` | The same structure with GROMOS atom names |

The course caps both ends of the peptide: an acetyl group (ACE) on Tyr1 and an amide group (NH2) on Tyr10. Free ends carry a +1 and a −1 charge that attract each other and distort pulling experiments (Zhao & Cheng, *Amino Acids* 2012); the caps remove those charges. The experimental reference values used in the labs (melting temperature, folding times) were measured on uncapped CLN025, so keep that difference in mind when you compare.

To cap a different NMR model, for example model 5:

```bash
python scripts/build_caps.py structures/2RVD.pdb 5 cln025_model5.pdb            # AMBER names
python scripts/build_caps.py structures/2RVD.pdb 5 cln025_model5.pdb --gromos   # GROMOS names
```

## Topology

AMBER, the course's main force field:

```bash
gmx pdb2gmx -f cln025_capped.pdb -o cln025.gro -p topol.top -i posre.itp \
            -ff amber99sb-ildn -water tip3p -ignh
```

GROMOS, the comparison force field. When asked, choose `None` for both termini: the caps already are the ends of the chain.

```bash
gmx pdb2gmx -f cln025_capped_gromos.pdb -o cln025.gro -p topol.top -i posre.itp \
            -ff gromos54a7 -water spc -ignh -ter
```

Both give a total charge of −2 (Asp3 and Glu5).

## `mdp/`

| File | Stage |
|---|---|
| `em.mdp` | Energy minimisation (steepest descent, until the largest force is below 100 kJ mol⁻¹ nm⁻¹) |
| `nvt.mdp` | 100 ps at constant volume, 300 K, protein restrained |
| `npt.mdp` | 500 ps at constant pressure, 300 K and 1 bar, protein restrained |
| `md.mdp` | Production at 300 K and 1 bar, no restraints |

All four use a 2 fs time step with bonds to hydrogen constrained, PME electrostatics with 1.0 nm cut-offs, the v-rescale thermostat on the whole system and the C-rescale barostat. They are written for amber99sb-ildn with TIP3P water. Copy them into your working folder and change only what a lab asks you to change, usually `nsteps` and `ref-t`.

## `scripts/`

| File | Use |
|---|---|
| `build_caps.py` | Caps one 2RVD model as ACE-YYDPETGTWY-NH2 (see above) |
