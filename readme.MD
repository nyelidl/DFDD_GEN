# DFDD Input Generator

A program for preparing initial files for **[DFDD](https://github.com/nyelidl/DFDD)** to perform **LB-PaCS-MD (Ligand-Binding PaCS-MD) via PaCS-Q** — a dynamic docking method that drives a ligand from outside the binding pocket toward its binding site using [PaCS-Q](https://github.com/nyelidl/PaCS-Q).

---

## What this does

```
Co-crystal PDB or protein PDB
        ↓
Receptor preparation (PDB2PQR, protonation, binding site detection)
        ↓
Ligand preparation (SMILES / prep+frcmod / mol2 / lib)
        ↓
BSV-guided ligand placement (30+ Å outside binding pocket)
        ↓
Topology generation (AMBER ff19SB + GAFF2, octahedral box)
        ↓
Short minimization + heating
        ↓
Download ready-to-run PaCS-MD package (.zip)
```

---

## Two versions

| Version | File | Use case |
|---|---|---|
| Pure Colab | `LB_input.ipynb` | Run entirely in Google Colab |
| Colab + Streamlit | `LB_input_streamlit.ipynb` | Nicer UI via Streamlit, shareable URL |

---

## Quick start

### Pure Colab
1. Open `LB_input.ipynb` in Google Colab
2. Runtime → Change runtime type → **GPU** (recommended)
3. Run cells top to bottom
4. Download `LB_PaCS_MD_inputs.zip` from Cell 15
5. Extract and run `./run_pacs_md.sh`

### Colab + Streamlit
1. Open `LB_input_streamlit.ipynb` in Google Colab
2. Runtime → Change runtime type → **GPU**
3. Run the launch cell
4. Open the ngrok URL in your browser
5. Follow the step-by-step UI

---

## Requirements

All dependencies are installed automatically by the notebook. No local installation needed.

### Software installed automatically
- AmberTools (antechamber, tleap, parmchk2)
- OpenBabel
- RDKit
- OpenFF Toolkit
- PDB2PQR
- OpenMM
- PaCS-Q
- NGLView / py3Dmol

### For HPC runs
- AMBER ≥ 22 (if using AMBER engine)
- PaCS-Q (`pip install PaCS-Q`)
- GROMACS (if using GROMACS conversion output)

---

## Supported systems

| System type | Supported |
|---|---|
| Soluble proteins | ✅ |
| Single chain | ✅ |
| Multiple chains (dimer, heterodimer) | ✅ |
| Disulfide bonds (intra + inter-chain) | ✅ |
| AMBER modified residues (PTR, SEP, TPO, HYP...) | ✅ |
| Known cofactors (NAD+, FAD, ATP, CoA...) | ✅ |
| Simple metal ions (Zn²⁺, Mg²⁺, Ca²⁺) | ✅ (nonbonded) |
| Small molecule ligands | ✅ |
| DNA / RNA containing systems | ❌ |
| Membrane proteins / lipid bilayers | ❌ |
| Heme / Fe-S cluster proteins | ❌ |
| Covalently bound metals | ❌ |

---

## Ligand input modes

| Mode | Input | Notes |
|---|---|---|
| SMILES | Type SMILES string | Auto-parameterized with GAFF2 + AM1-BCC |
| prep + frcmod | Upload both files | From previous AMBER parameterization |
| mol2 | Upload .mol2 file | Charges detected automatically |
| lib | Upload .lib file | frcmod requested if non-standard atom types |

All modes force resname → **MOL** for consistency with PaCS-Q.

---

## Binding Site Vectors (BSV)

Ligand placement uses the BSV method from:

> Kuvek et al., *J. Chem. Inf. Model.* **2026**, 66, 1604–1616  
> DOI: [10.1021/acs.jcim.5c02705](https://doi.org/10.1021/acs.jcim.5c02705)

An icosphere is anchored at the binding site center. Vectors are traced to the nearest protein atom surface, determining the pocket exit direction. The ligand is placed along this exit vector at ≥ 30 Å — the minimum distance recommended for LB-PaCS-MD.

---

## Output files

```
LB_PaCS_MD_inputs.zip
├── system/
│   ├── complex.top          AMBER topology
│   ├── complex.inpcrd       Initial coordinates
│   └── last_frame.rst7      Minimized restart file
├── components/
│   ├── protein_prepared.pdb Protonated protein
│   ├── ligand_placed.pdb    BSV-positioned ligand
│   ├── guest.prep           GAFF2 prep file
│   └── guest.frcmod         GAFF2 frcmod file
├── config.json              OpenMM engine config
│   OR md.mdin               AMBER engine config
│   OR md.mdp                GROMACS MD parameters*
├── run_pacs_md.sh           Ready-to-run script
├── summary.txt              Full run parameters
└── README.md                Instructions
```

*GROMACS output is for manual MD workflow only — PaCS-Q does not support GROMACS engine.

---

## Running PaCS-MD

### Local machine
```bash
unzip LB_PaCS_MD_inputs.zip
cd LB_PaCS_MD_inputs
chmod +x run_pacs_md.sh
./run_pacs_md.sh
```

### HPC — SLURM
```bash
# Transfer files
scp -r LB_PaCS_MD_inputs/ user@cluster.hpc.ac.th:~/runs/

# Submit
cd ~/runs/LB_PaCS_MD_inputs
sbatch run_pacs_md.sh
```

### HPC — PBS
```bash
qsub run_pacs_md.sh
```

---

## PaCS-Q parameters

| Parameter | Default | Range | Notes |
|---|---|---|---|
| Cycles (-cy) | 20 | 10–100 | More = more complete unbinding pathway |
| Candidates (-cd) | 5 | 3–20 | More = better sampling, higher cost |
| Sim time/cycle | 20 ps | 10–100 ps | Time per MD cycle |
| Mode | -b | fixed | Binding only (ligand → pocket) |
| Engine | OpenMM or AMBER | — | GROMACS = conversion only |

---

## Simulation settings (fixed defaults)

| Setting | Value | Reason |
|---|---|---|
| Force field (protein) | ff19SB | Current AMBER standard |
| Force field (ligand) | GAFF2 + AM1-BCC | Standard small molecule |
| Water model | TIP3P | Consistent with ff19SB |
| Box type | Octahedral | ~30% fewer waters vs rectangular |
| Box padding | 15 Å | Accounts for 30+ Å ligand displacement |
| Timestep | 2 fs | Safe with SHAKE |
| SHAKE | H-bonds | Standard |
| Cutoff | 10 Å | Standard AMBER |
| Electrostatics | PME | Standard periodic |
| Thermostat | Langevin (γ = 2.0) | Standard |
| Barostat | Monte Carlo | Standard |
| Pressure | 1 bar | Standard |
| Neutralization | Na⁺/Cl⁻ (neutralize only) | No extra salt |
| Restraints | None | Protein moves freely |
| Placement distance | ≥ 30 Å | LB-PaCS-MD minimum |

---

## After PaCS-MD

Analysis is the user's responsibility. Suggested tools:

```bash
# Trajectory imaging
cpptraj -p system/complex.top \
        -y MDrun/cycle_001/traj.nc \
        -x traj_imaged.nc \
        autoimage

# Concatenate all cycles
# Binding free energy: MBAR, WHAM, umbrella sampling
# Visualization: VMD or PyMOL
```

**Reference:** Harada & Kitao, *J. Chem. Theory Comput.* **2015**, 11, 5032–5042

---

## Known limitations

- Soluble proteins only (no membrane, no DNA/RNA)
- Binding mode only (ligand placed outside pocket)
- No heme or metal cluster support
- No HMR (hydrogen mass repartitioning)
- No covalent ligands
- GROMACS output requires manual PaCS-MD cycling

---

## Citation

If you use this workflow, please cite:

**LB-PaCS-MD method:**
> Harada R, Kitao A. *Biased Molecular Simulations for Open-Boundary Systems.* J. Chem. Theory Comput. 2015, 11, 5032–5042.

**Binding Site Vectors:**
> Kuvek T, Jandová Z, Schleifer KJ, Oostenbrink C. *Binding Site Vectors Enable Mapping of Cytochrome P450 Functional Landscapes.* J. Chem. Inf. Model. 2026, 66, 1604–1616.

**PaCS-Q:**
> [https://github.com/nyelidl/PaCS-Q](https://github.com/nyelidl/PaCS-Q)

**AMBER force fields:**
> He X et al. *ff19SB: Amino-Acid-Specific Protein Backbone Parameters Trained against Quantum Mechanics Energy Surfaces in Solution.* J. Chem. Theory Comput. 2020, 16, 528–552.

> Wang J et al. *Development and Testing of a General Amber Force Field.* J. Comput. Chem. 2004, 25, 1157–1174. (GAFF/GAFF2)

---

## Contact

For issues and questions, please open a GitHub issue.

---

## License

MIT License — see [LICENSE](LICENSE) for details.
