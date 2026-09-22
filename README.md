# Structural Impact of Missense Mutation in Ubiquitin-Associated (UBA) Domain of Human Rad23 Protein: Molecular Dynamics Simulations Using Different Force-Fields

[![GROMACS](https://img.shields.io/badge/GROMACS-2019%20(GPU%2FMPI)-blue.svg)](https://www.gromacs.org/)
[![Python](https://img.shields.io/badge/Python-3.x-brightgreen.svg)](https://www.python.org/)
[![HPC](https://img.shields.io/badge/HPC-SLURM%20Cluster-orange.svg)](https://slurm.schedmd.com/)
[![License](https://img.shields.io/badge/License-Academic-lightgrey.svg)](#)

> **M.Tech. Thesis Project**  
> **Department:** Department of Biological Sciences and Bioengineering (BSBE)  
> **Institution:** Indian Institute of Technology (IIT) Kanpur  
> **Compute Facilities:** HPC2013 & Param Sanganak Supercomputing Facility, IIT Kanpur  

---

## Table of Contents
- [Overview & Biological Context](#overview--biological-context)
- [Motivation & The AlphaFold2 Dilemma](#motivation--the-alphafold2-dilemma)
- [Objectives](#objectives)
- [Simulation Systems & Protocols](#simulation-systems--protocols)
  - [Systems Overview](#systems-overview)
  - [Equilibration & Production Workflow](#equilibration--production-workflow)
- [Key Scientific Findings](#key-scientific-findings)
  - [1. Secondary Structure Dynamics (DSSP)](#1-secondary-structure-dynamics-dssp)
  - [2. Structural Deviations & Flexibility (RMSD & RMSF)](#2-structural-deviations--flexibility-rmsd--rmsf)
  - [3. Domain Compactness (Radius of Gyration)](#3-domain-compactness-radius-of-gyration)
  - [4. Hydrophobic Core Hydration](#4-hydrophobic-core-hydration)
- [Conclusions](#conclusions)
- [Repository Structure](#repository-structure)
- [Analysis Notebooks & HPC Scripts](#analysis-notebooks--hpc-scripts)
- [Reproducibility & Requirements](#reproducibility--requirements)
- [Acknowledgments](#acknowledgments)

---

## Overview & Biological Context

Human Rad23 (hHR23a) is a multifunctional protein vital to **Nucleotide Excision Repair (NER)** and **Ubiquitin-Proteasome System (UPS)** pathways. Acting as an escort / shuttle factor, Rad23 binds polyubiquitinated substrate proteins and delivers them directly to the 26S proteasome for degradation.

Its modular architecture comprises:
- **UBL domain** (Ubiquitin-Like): Coordinates interaction with the 26S proteasome subunit.
- **XPCB domain**: Interacts with the Xeroderma Pigmentosum group C protein during DNA repair.
- **UBA1 & UBA2 domains** (Ubiquitin-Associated): Conserved three-helix bundles ($\alpha_1, \alpha_2, \alpha_3$) presenting hydrophobic surface patches that recognize and bind mono- and poly-ubiquitin chains.

The target system in this study is the **UBA1 domain** (residues 161–200 from PDB ID: **1OQY**).

```
   Rad23 Architecture:
   [ N-term ]──( UBL )──[ Linker ]──( UBA1: 161-200 )──[ Linker ]──( XPCB )──[ Linker ]──( UBA2 )──[ C-term ]
                                             │
                                      L198A Point Mutation
                                (Hydrophobic Core Destabilization)
```

---

## Motivation & The AlphaFold2 Dilemma

While deep-learning-based structural prediction tools such as **AlphaFold2 (AF2)** have revolutionized structural biology, they exhibit documented deficiencies in predicting the destabilizing effects of single-point missense mutations. Because AF2 relies heavily on sequence homology in Multiple Sequence Alignments (MSAs) and static PDB structures:

- **AF2 Prediction for L198A Mutation:** Superimposing the predicted Wild-Type (WT) and L198A mutant UBA1 structures produces a C$\alpha$ RMSD of only **0.1 Å** (nearly identical, fully folded helical bundles).
- **Physical Reality:** In experimental contexts, mutating the buried, highly conserved leucine residue (**Leu198**) at the hydrophobic core of the three-helix bundle to alanine (**Ala198**) disrupts core packing, creates internal voids, and triggers conformational destabilization and disorder.

This project employs **All-Atom Molecular Dynamics (MD) Simulations** across multiple classical and modern force fields to:
1. Capture real-time conformational destabilization, unfolding transitions, and solvent penetration caused by L198A.
2. Evaluate which biomolecular force fields are best equipped to model mutational effects and transitions toward disordered or partially folded states.

---

## Objectives

1. **Elucidate Mechanism of Destabilization:** Determine how replacing Leu198 with Ala198 alters hydrophobic core packing, tertiary stability, and intramolecular hydrogen bonding.
2. **Comparative Force Field Benchmark:** Simulate the Wild-Type and L198A mutant across four prominent force fields:
   - **AMBER ff99SB-ILDN**
   - **CHARMM36**
   - **CHARMM36m** (specifically reparameterized for Intrinsically Disordered Proteins / Regions with refined protein-water interactions)
   - **OPLS-AA/L**
3. **Assess Force Field Sensitivity:** Determine which force fields permit structural destabilization versus those exhibiting excessive helical rigidity that prevent sampling of unfolded/disordered states.

---

## Simulation Systems & Protocols

### Systems Overview

A total of **8 simulation systems** (4 force fields $\times$ 2 protein states) were constructed and simulated for **1000 ns (1.0 µs)** each, yielding an aggregate of **8.0 µs** of production MD trajectories.

| System | State | Force Field | Water Model | Box Shape | Ions | Temp / Press | Production Time |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **WT-AMBER** | Wild-Type | AMBER ff99SB-ILDN | TIP3P | Cubic ($d = 1.0\text{ nm}$) | Na⁺ / Cl⁻ (Neutral) | 300 K / 1 bar | 1000 ns (1.0 µs) |
| **MT-AMBER** | L198A Mutant | AMBER ff99SB-ILDN | TIP3P | Cubic ($d = 1.0\text{ nm}$) | Na⁺ / Cl⁻ (Neutral) | 300 K / 1 bar | 1000 ns (1.0 µs) |
| **WT-CHARMM36** | Wild-Type | CHARMM36 | TIP3P | Cubic ($d = 1.0\text{ nm}$) | Na⁺ / Cl⁻ (Neutral) | 300 K / 1 bar | 1000 ns (1.0 µs) |
| **MT-CHARMM36** | L198A Mutant | CHARMM36 | TIP3P | Cubic ($d = 1.0\text{ nm}$) | Na⁺ / Cl⁻ (Neutral) | 300 K / 1 bar | 1000 ns (1.0 µs) |
| **WT-C36m** | Wild-Type | CHARMM36m | TIP3P | Cubic ($d = 1.0\text{ nm}$) | Na⁺ / Cl⁻ (Neutral) | 300 K / 1 bar | 1000 ns (1.0 µs) |
| **MT-C36m** | L198A Mutant | CHARMM36m | TIP3P | Cubic ($d = 1.0\text{ nm}$) | Na⁺ / Cl⁻ (Neutral) | 300 K / 1 bar | 1000 ns (1.0 µs) |
| **WT-OPLS** | Wild-Type | OPLS-AA/L | SPC/E | Cubic ($d = 1.0\text{ nm}$) | Na⁺ / Cl⁻ (Neutral) | 300 K / 1 bar | 1000 ns (1.0 µs) |
| **MT-OPLS** | L198A Mutant | OPLS-AA/L | SPC/E | Cubic ($d = 1.0\text{ nm}$) | Na⁺ / Cl⁻ (Neutral) | 300 K / 1 bar | 1000 ns (1.0 µs) |

---

### Equilibration & Production Workflow

The simulations followed a rigorous multi-stage restraint-relaxation protocol executed with GROMACS:

```
  ┌────────────────────────────────────────────────────────┐
  │ 1. System Preparation                                  │
  │    PDB clean -> pdb2gmx -> cubic box -> solvate -> ions│
  └──────────────────────────┬─────────────────────────────┘
                             │
  ┌──────────────────────────▼─────────────────────────────┐
  │ 2. Energy Minimization                                 │
  │    Steepest descent (em.mdp) until Fmax < 1000 kJ/mol  │
  └──────────────────────────┬─────────────────────────────┘
                             │
  ┌──────────────────────────▼─────────────────────────────┐
  │ 3. NVT Equilibration (100 ps, 300 K)                   │
  │    Harmonic position restraints: 4000 kJ/mol/nm²       │
  └──────────────────────────┬─────────────────────────────┘
                             │
  ┌──────────────────────────▼─────────────────────────────┐
  │ 4. Multi-Stage NPT Restraint Release (300 K, 1 bar)    │
  │    4000 -> 2000 -> 1000 -> 500 -> 250 -> 100 -> 50     │
  │    -> 10 -> 0 kJ/mol/nm²                               │
  └──────────────────────────┬─────────────────────────────┘
                             │
  ┌──────────────────────────▼─────────────────────────────┐
  │ 5. Unrestrained Relaxation Run                         │
  │    100 ns equilibration without position restraints    │
  └──────────────────────────┬─────────────────────────────┘
                             │
  ┌──────────────────────────▼─────────────────────────────┐
  │ 6. Production MD Simulation                            │
  │    1000 ns (1.0 µs) production run per system          │
  └──────────────────────────┬─────────────────────────────┘
                             │
  ┌──────────────────────────▼─────────────────────────────┐
  │ 7. Trajectory Post-Processing & Analysis               │
  │    PBC removal (nojump, mol center) -> fitting ->      │
  │    RMSD, RMSF, DSSP, H-bonds, Rg, Water penetration    │
  └────────────────────────────────────────────────────────┘
```

---

## Key Scientific Findings

### 1. Secondary Structure Dynamics (DSSP)
- **Wild-Type Systems:** All four force fields maintained stable $\alpha$-helical content throughout the full 1000 ns.
- **L198A Mutant Systems:**
  - **AMBER ff99SB-ILDN:** Captured severe structural unraveling. A significant fraction of helical residues transitioned into turns and disordered coil conformations.
  - **OPLS-AA/L:** Showed noticeable partial loosening of helices.
  - **CHARMM36 & CHARMM36m:** Exhibited strong helical propensity, retaining predominantly intact secondary structures despite the core mutation.

### 2. Structural Deviations & Flexibility (RMSD & RMSF)
- **C$\alpha$ Backbone RMSD:**
  - Wild-type systems remained stable with low fluctuations ($\approx 1.5 - 2.5\text{ \AA}$).
  - Mutant AMBER ff99SB-ILDN exhibited large deviations ($> 6\text{ \AA}$), indicating extensive conformational departure from the native starting structure.
- **Per-Residue RMSF:**
  - WT systems exhibited localized, minimal fluctuations primarily at terminal regions.
  - Mutant AMBER ff99SB-ILDN displayed pronounced fluctuations spanning all residues, reflecting generalized loss of rigidity.

### 3. Domain Compactness (Radius of Gyration)
- **AMBER ff99SB-ILDN Mutant:** Marked increase in $R_g$, reflecting domain expansion and unwinding of the core.
- **OPLS-AA/L Mutant:** Maintained lower $R_g$ than its WT counterpart, sampling collapsed/compact non-native conformations.
- **CHARMM36 / CHARMM36m:** Maintained uniform $R_g$ profiles corresponding to compact helical bundles.

### 4. Hydrophobic Core Hydration
Water molecules within a **4.0 Å cutoff sphere** around residue 198 were quantified:
- **Wild-Type (Leu198):** Solvent-excluded, deeply buried hydrophobic core; water count remained consistently near zero.
- **Mutant (Ala198):** The smaller alanine side chain created a packing defect, allowing water molecules to penetrate the core:
  - **AMBER ff99SB-ILDN:** Rapid and sustained water intrusion with high fluctuations.
  - **OPLS-AA/L & CHARMM36:** Exhibited sudden bursts of water entry after $\sim 550\text{ ns}$ as conformational fluctuations opened the pocket.

---

## Conclusions

1. **Mechanism of Destabilization:** Truncating Leu198 to Ala198 destabilizes the UBA1 domain by disrupting close hydrophobic packing. The resulting cavity permits water entry, weakening tertiary inter-helical contacts and inducing local unwinding.
2. **Force Field Sensitivity:**
   - **AMBER ff99SB-ILDN** demonstrated the highest sensitivity to the destabilizing point mutation, successfully capturing the transition toward disordered conformations.
   - **OPLS-AA/L** also identified structural loosening, making it a viable secondary candidate.
   - **CHARMM36 and CHARMM36m** showed higher-than-expected helical stabilization, resisting unfolding under identical simulation conditions.
3. **Implications for Computational Biology:** When modeling missense mutations that induce disorder or destabilization, force fields with strong helical bias may obscure critical destabilization phenotypes.

---

## Repository Structure

```
thesis-work/
│
├── README.md                           # Comprehensive project documentation
├── DefensePresentation.pdf             # Official M.Tech. defense presentation slides
│
└── Code/
    ├── hpcscripts/                     # HPC job submission & trajectory processing scripts
    │   ├── Jobsubmission/              # SLURM submission scripts for production MD
    │   │   ├── script_mam.sh           # Mutant AMBER (1 µs MD + equilibration pipeline)
    │   │   ├── script_mc36m.sh         # Mutant CHARMM36m (1 µs MD + equilibration pipeline)
    │   │   └── script_mch.sh           # Mutant CHARMM36 (1 µs MD + equilibration pipeline)
    │   │
    │   └── preliminaryanalysis/        # GROMACS post-processing & extraction scripts
    │       ├── script_analyze.sh       # Trajectory PBC correction (nojump/fit) & RMSD/RMSF
    │       ├── script_analyz1.sh       # WT Amber index generation, helical RMSD & H-bond extraction
    │       ├── script_analyz2.sh       # WT CHARMM36m index generation, helical RMSD & H-bond extraction
    │       ├── analyz_mch.sh           # MT CHARMM36 analysis pipeline
    │       ├── hbond_helical.sh        # Helical H-bond extraction across all force fields
    │       ├── rmsf_500.sh             # 500-1000 ns RMSF window calculations
    │       └── struct.sh               # Frame extraction at 50 ns intervals for structural rendering
    │
    └── Analysis/                       # Jupyter notebooks for statistical analysis & visualization
        ├── rmsd.ipynb                  # C-alpha and backbone RMSD trajectory plots
        ├── rmsf.ipynb                  # Per-residue fluctuation profiles (RMSF)
        ├── ss.ipynb                    # DSSP secondary structure evolution heatmaps
        ├── hbond.ipynb                 # Intramolecular and helical hydrogen bonding
        ├── Rg.ipynb                    # Radius of gyration (compactness) analysis
        ├── water.ipynb                 # Hydration analysis (4 Å sphere around residue 198)
        └── water3.ipynb                # Late-trajectory (700-1000 ns) water penetration dynamics
```

---

## Analysis Notebooks & HPC Scripts

| Notebook / Script | Description | Primary Libraries / Tools |
| :--- | :--- | :--- |
| [`rmsd.ipynb`](Code/Analysis/rmsd.ipynb) | Trajectory RMSD calculation, moving average smoothing, WT vs MT comparisons | `MDAnalysis`, `pandas`, `matplotlib` |
| [`rmsf.ipynb`](Code/Analysis/rmsf.ipynb) | Per-residue flexibility curves highlighting dynamic regions | `pandas`, `matplotlib`, `seaborn` |
| [`ss.ipynb`](Code/Analysis/ss.ipynb) | Matrix processing of DSSP `.xpm` matrices to visualize $\alpha$-helix, coil, turn transitions | `pandas`, `plotly`, `matplotlib` |
| [`hbond.ipynb`](Code/Analysis/hbond.ipynb) | Quantifies total and intra-helical hydrogen bonds over 1 µs | `pandas`, `matplotlib`, `numpy` |
| [`Rg.ipynb`](Code/Analysis/Rg.ipynb) | Radius of gyration about center of mass evaluating domain expansion | `pandas`, `matplotlib`, `seaborn` |
| [`water.ipynb`](Code/Analysis/water.ipynb) | Hydration shell counts within 4 Å of residue 198 across the trajectory | `pandas`, `matplotlib` |
| [`script_mam.sh`](Code/hpcscripts/Jobsubmission/script_mam.sh) | Full SLURM MPI submission script for GROMACS simulation pipeline | GROMACS 2019, MPI, CUDA |

> **Note on Raw Data:**  
> Raw coordinate trajectories (`.xtc`), energy files (`.edr`), and checkpoints (`.cpt`) occupy hundreds of gigabytes and are stored on the institutional HPC cluster storage. The notebooks in `Code/Analysis/` process the derived `.xvg` data tables and `.xpm` matrices generated by GROMACS.

---

## Reproducibility & Requirements

To view or rerun the analytical notebooks:

### Python Environment
```bash
pip install pandas numpy matplotlib seaborn MDAnalysis plotly openpyxl
```

### HPC & Simulation Environment
- **GROMACS**: Version 2019 or later (configured with MPI and GPU support)
- **DSSP**: Version 3.x+ (for secondary structure assignment)
- **Workload Manager**: SLURM (or PBS/Torque with script directive adaptation)

---

## Acknowledgments

- **HPC2013 and Param Sanganak Supercomputing Facility** at IIT Kanpur for computational grants and technical support.
- Faculty, staff, and computational biology lab members at the **Department of Biological Sciences and Bioengineering (BSBE), IIT Kanpur**.
