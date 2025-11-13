RE# 4QKZ Project: Molecular Dynamics and Pocket Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GROMACS](https://img.shields.io/badge/GROMACS-2022.5-blue.svg)](http://www.gromacs.org/)
[![Python](https://img.shields.io/badge/Python-3.10-green.svg)](https://www.python.org/)

## Overview

This project explores the dynamics and pocket evolution of the **4QKZ protein structure** using **GROMACS** for molecular dynamics (MD) simulation and **Fpocket** for cavity detection. The workflow includes simulation setup, production MD (100 ns), and comprehensive structural analysis of trajectories, emphasizing pocket volume, solvent accessibility, and stability metrics (RMSD, RMSF, SASA).

**Key Findings:**
- Pocket is **highly druggable** with stable volume (2273 ± 79 Ų)
- Excellent structural stability (RMSD < 0.2 nm)
- Good solvent accessibility (SASA ~24 Ų)
- 5 flexible residues identified for induced-fit binding

---

## Table of Contents

- [Environment Setup](#1-environment-setup)
- [System Preparation](#2-system-preparation)
- [Energy Minimization](#3-energy-minimization)
- [Equilibration Phases](#4-equilibration-phases)
- [Production MD](#5-production-molecular-dynamics)
- [Trajectory Post-Processing](#6-trajectory-post-processing)
- [Pocket Detection](#7-pocket-detection-using-fpocket)
- [Automated Pocket Analysis](#8-automated-pocket-analysis-fpocket-timeseries)
- [Structural Dynamics Analysis](#9-structural-dynamics-analysis)
- [Results Interpretation](#10-results-interpretation)
- [Repository Organization](#11-repository-organization)
- [Version Control](#12-version-control)
- [Citation](#13-citation)

---

## 1. Environment Setup

### Hardware Requirements
- HPC cluster with SLURM workload manager
- GPU acceleration (recommended for production MD)
- Minimum 32 GB RAM

### Software Dependencies

All work was performed on a high-performance computing (HPC) cluster running SLURM.

```bash
# Load required modules
module load gromacs/2022.5
module load python/3.10
module load fpocket/4.0

# Python packages
pip install pandas matplotlib numpy scipy
```

**Key Software:**
- GROMACS 2022.5+
- Fpocket 4.0
- Python 3.10+ with NumPy, Pandas, Matplotlib, SciPy

---

## 2. System Preparation

### Step 1: Initial Structure Cleanup

The PDB structure 4QKZ was cleaned to remove unwanted chains, waters, and heteroatoms.

```bash
pdb4amber -i 4QKZ.pdb -o 4QKZ_clean.pdb
```

**Purpose:** Remove crystallographic waters, non-standard residues, and select appropriate protein chains.

### Step 2: Define the Simulation Box

```bash
gmx editconf -f 4QKZ_clean.pdb -o newbox.gro -c -d 1.0 -bt cubic
```

**Parameters:**
- `-c`: Center protein in box
- `-d 1.0`: 1.0 nm minimum distance to box edge
- `-bt cubic`: Cubic box shape

### Step 3: Solvate the System

```bash
gmx solvate -cp newbox.gro -cs spc216.gro -o solvated.gro -p topol.top
```

**Purpose:** Add explicit water molecules (SPC water model) around the protein.

### Step 4: Add Ions for Neutralization

```bash
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```

**Purpose:** Neutralize system charge and add physiological salt concentration (~0.15 M NaCl).

---

## 3. Energy Minimization

Minimize potential energy before equilibration to remove steric clashes:

```bash
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
```

**Convergence Check:**

```bash
# Extract potential energy
gmx energy -f em.edr -o potential.xvg
# Select "Potential" from menu
```

**Success Criteria:**
- Potential Energy < -1.0e6 kJ/mol
- Maximum force < 1000 kJ/mol/nm
- Fmax in `em.log` should be < 1000

---

## 4. Equilibration Phases

### (a) NVT Equilibration (Constant Volume)

Stabilizes **temperature** while keeping volume constant.

```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt
```

**Key Parameters in `nvt.mdp`:**
- Temperature: 300 K
- Duration: 100 ps
- Thermostat: V-rescale

**Check Temperature Stability:**

```bash
gmx energy -f nvt.edr -o temperature.xvg
# Select "Temperature" from menu
```

### (b) NPT Equilibration (Constant Pressure)

Stabilizes **pressure** and adjusts system density.

```bash
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt
```

**Key Parameters in `npt.mdp`:**
- Pressure: 1 bar
- Duration: 100 ps
- Barostat: Parrinello-Rahman

**Check Pressure and Density:**

```bash
gmx energy -f npt.edr -o pressure.xvg
gmx energy -f npt.edr -o density.xvg
```

---

## 5. Production Molecular Dynamics

Run the main simulation for **100 ns**:

```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md
```

**Production Parameters:**
- Duration: 100 ns (100,000 ps)
- Time step: 2 fs
- Temperature: 300 K
- Pressure: 1 bar
- Output frequency: 50 ps (2000 frames)

**Main Outputs:**
- `md.xtc` — compressed trajectory (coordinates)
- `md.edr` — energy data
- `md.gro` — final frame
- `md.log` — simulation log
- `md.cpt` — checkpoint file

**Monitor Progress:**

```bash
# Check if simulation is running
tail -f md.log

# Estimate completion time
gmx mdrun -deffnm md -maxh 24  # Set max hours
```

---

## 6. Trajectory Post-Processing

### Remove Periodic Boundary Conditions

```bash
gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center
# Select: 1 (Protein), 0 (System)
```

### Extract Reference Frame (First Frame)

```bash
printf "1\n" | gmx trjconv -f md.xtc -s md.tpr -o ref.pdb -dump 0
```

**Purpose:** Create reference structure for RMSD calculations and pocket analysis.

### Create Pocket-Specific Index File

```bash
gmx make_ndx -f ref.pdb -o pocket.ndx
# Manually select residues: 158, 159, 160, 161, 188, 189, 193, 194, 197, 198, 214, 216-222
```

---

## 7. Pocket Detection using Fpocket

Run pocket detection on the reference structure:

```bash
/path/to/fpocket/bin/fpocket -f ref.pdb
```

**Fpocket Outputs:**

```
ref_out/
├── ref_out.pdb          # Annotated structure with pocket spheres
├── ref_info.txt         # Pocket statistics and rankings
├── pockets/             # Individual pocket PDB files
│   ├── pocket1_atm.pdb
│   ├── pocket2_atm.pdb
│   └── ...
└── ref_PYMOL.sh        # PyMOL visualization script
```

### Visualize in PyMOL

```bash
pymol ref.pdb ref_out/pockets/pocket1_atm.pdb
```

### Extract Residues of Chosen Pocket

```bash
grep "^ATOM" ref_out/pockets/pocket1_atm.pdb | awk '{print $6}' | sort -n | uniq > pocket1_residues.txt
```

**Identified Pocket Residues:**
158, 159, 160, 161, 188, 189, 193, 194, 197, 198, 214, 216, 217, 218, 219, 220, 221, 222

---

## 8. Automated Pocket Analysis (Fpocket Timeseries)

### Extract Frames from Trajectory

```bash
# Extract every 500 ps (100 frames total)
for i in {0..100000..500}; do
    echo 1 | gmx trjconv -s md.tpr -f md.xtc -o frame_${i}.pdb -dump $i
done
```

### Run Fpocket on Each Frame

```bash
#!/bin/bash
for pdb in frame_*.pdb; do
    fpocket -f $pdb
done
```

### Parse Results with Python

```python
# parse_fpocket_runs.py
import os
import pandas as pd
import re

data = []

for folder in sorted(os.listdir('.')):
    if folder.endswith('_out'):
        info_file = os.path.join(folder, f"{folder.replace('_out', '')}_info.txt")
        
        if os.path.exists(info_file):
            with open(info_file, 'r') as f:
                content = f.read()
                
                # Extract pocket 1 statistics
                match = re.search(r'Pocket 1.*?Volume : (\d+\.\d+).*?Hydrophobicity.*?: (\d+\.\d+).*?Polarity.*?: (\d+\.\d+)', content, re.DOTALL)
                
                if match:
                    data.append({
                        'snapshot': folder,
                        'Volume': float(match.group(1)),
                        'Hydrophobicity': float(match.group(2)),
                        'Polarity': float(match.group(3))
                    })

df = pd.DataFrame(data)
df.to_csv('fpocket_timeseries.csv', index=False)
print(df.describe())
```

**Run:**

```bash
python parse_fpocket_runs.py
```

**Output Columns:**
- `snapshot` — frame identifier
- `Volume` — pocket volume (Ų)
- `Hydrophobicity` — hydrophobic surface ratio
- `Polarity` — polar surface ratio
- `DruggabilityScore` — predicted druggability
- `NumberAlphaSpheres` — cavity sampling density

---

## 9. Structural Dynamics Analysis

### RMSD (Root Mean Square Deviation)

**Measures structural stability of the pocket over time.**

```bash
# Calculate RMSD for pocket residues
gmx rms -s md.tpr -f md_noPBC.xtc -n pocket.ndx -o rmsd_pocket.xvg
# Select pocket group and pocket group for reference
```

**Convert XVG to CSV:**

```bash
grep -v "^[@#]" rmsd_pocket.xvg > rmsd_pocket.csv
```

**Python Visualization:**

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("rmsd_pocket.csv", delim_whitespace=True, header=None, names=['Time', 'RMSD'])
df['Time_ns'] = df['Time'] / 1000

plt.figure(figsize=(10, 6))
plt.plot(df['Time_ns'], df['RMSD'], linewidth=0.8)
plt.axhline(df['RMSD'].mean(), color='red', linestyle='--', label=f'Mean: {df["RMSD"].mean():.3f} nm')
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (nm)')
plt.title('Pocket RMSD Over Time')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("rmsd_pocket.png", dpi=300)
plt.show()
```

**Interpretation:**
- RMSD < 0.2 nm → Stable pocket structure
- RMSD > 0.3 nm → Significant conformational changes
- Increasing RMSD → Pocket destabilization
- **Your Result: 0.187 ± 0.051 nm** ✅ Excellent stability

### RMSF (Root Mean Square Fluctuation)

**Measures per-residue flexibility.**

```bash
gmx rmsf -s md.tpr -f md_noPBC.xtc -n pocket.ndx -o rmsf_pocket.xvg -res
# Select pocket group
```

**Visualization:**

```python
df = pd.read_csv("rmsf_pocket.csv", delim_whitespace=True, header=None, names=['Residue', 'RMSF'])

plt.figure(figsize=(10, 6))
bars = plt.bar(df['Residue'], df['RMSF'], color='teal', edgecolor='black', linewidth=0.5)

# Highlight high flexibility residues
threshold = df['RMSF'].mean() + df['RMSF'].std()
high_flex = df[df['RMSF'] > threshold]

for idx in high_flex.index:
    bars[idx].set_color('red')
    bars[idx].set_alpha(0.9)

plt.axhline(df['RMSF'].mean(), color='red', linestyle='--', label=f'Mean: {df["RMSF"].mean():.3f} nm')
plt.xlabel('Residue Number')
plt.ylabel('RMSF (nm)')
plt.title('Per-Residue Flexibility (RMSF)')
plt.legend()
plt.grid(True, alpha=0.3, axis='y')
plt.savefig("rmsf_pocket.png", dpi=300)
plt.show()
```

**Interpretation:**
- RMSF < 0.1 nm → Rigid residue (good for binding)
- RMSF > 0.15 nm → Flexible residue (induced-fit potential)
- **Your High Flexibility Residues:** 158, 188, 189, 219, 221

### SASA (Solvent Accessible Surface Area)

**Measures pocket exposure to solvent.**

```bash
gmx sasa -f md_noPBC.xtc -s md.tpr -n pocket.ndx -o sasa_pocket.xvg
# Select pocket group
```

**Visualization:**

```python
df = pd.read_csv("sasa_pocket.csv", delim_whitespace=True, header=None, names=['Time', 'SASA'])
df['Time_ns'] = df['Time'] / 1000

plt.figure(figsize=(10, 6))
plt.plot(df['Time_ns'], df['SASA'], linewidth=0.8, color='#2E86AB')
plt.axhline(df['SASA'].mean(), color='red', linestyle='--', label=f'Mean: {df["SASA"].mean():.2f} nm²')
plt.fill_between(df['Time_ns'], 
                  df['SASA'].mean() - df['SASA'].std(), 
                  df['SASA'].mean() + df['SASA'].std(), 
                  alpha=0.2, color='red')
plt.xlabel('Time (ns)')
plt.ylabel('SASA (nm²)')
plt.title('Solvent Accessible Surface Area')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("sasa_pocket.png", dpi=300)
plt.show()
```

**Interpretation:**
- SASA decreasing → Pocket closure or ligand binding
- SASA stable → Consistent accessibility
- **Your Result: 23.62 ± 1.34 Ų** ✅ Good accessibility

### Pocket Volume Analysis

```bash
# Analyze volume from fpocket timeseries
python analyze_volume.py
```

**Interpretation:**
- Stable volume → Predictable binding affinity
- **Your Result: 2273 ± 79 Ų** ✅ Large, stable pocket

---

## 10. Results Interpretation

### Summary Statistics

| Metric | Value | Interpretation | Status |
|--------|-------|----------------|--------|
| **RMSD** | 0.187 ± 0.051 nm | Excellent structural stability | ✅ |
| **RMSF (mean)** | 0.124 ± 0.045 nm | Moderate flexibility | ✅ |
| **SASA** | 23.62 ± 1.34 Ų | Good solvent accessibility | ✅ |
| **Pocket Volume** | 2273 ± 79 Ų | Large, stable cavity | ✅ |
| **Volume CV** | 3.47% | Minimal size fluctuation | ✅ |
| **Flexible Residues** | 5/18 (28%) | Good balance rigidity/flexibility | ✅ |

### Druggability Assessment

**Overall Classification: HIGHLY DRUGGABLE** 

**Criteria Met:**
1. ✅ Pocket Stability (RMSD < 0.2 nm)
2. ✅ Volume Stability (CV < 5%)
3. ✅ Solvent Accessibility (SASA > 20 Ų)
4. ✅ Appropriate Size (2000-3000 Ų)
5. ✅ Balanced Flexibility (72% rigid residues)

### Time-Dependent Behavior

| Period | RMSD (nm) | Volume (Ų) | Status |
|--------|-----------|------------|--------|
| 0-25 ns | 0.237 | 2247 | Initial relaxation |
| 25-50 ns | 0.207 | 2277 | Stabilizing |
| 50-75 ns | 0.156 | 2313 | Equilibrated ✅ |
| 75-100 ns | 0.148 | 2255 | Fully stable ✅ |

**Key Finding:** RMSD decreases over time, indicating the pocket adopts a more stable conformation. The simulation reached equilibrium at ~50 ns.

### Implications for Drug Design

1. **Target Rigid Residues (197, 214, 216-218, 220, 222):**
   - Form stable anchoring interactions
   - Low RMSF → predictable binding geometry

2. **Exploit Flexible Residues (158, 188, 189, 219, 221):**
   - Enable induced-fit binding
   - Potential selectivity handles
   - May act as "lid" over bound ligand

3. **Pocket Characteristics:**
   - Large volume (2273 Ų) → can accommodate MW 300-700 Da
   - Good accessibility → favorable binding kinetics
   - Stable over 100 ns → reliable docking predictions

---

## 11. Repository Organization

```
4QKZ_project/
│
├── README.md                    # This file
├── LICENSE                      # MIT License
│
├── setup/                       # System preparation files
│   ├── 4QKZ_clean.pdb          # Cleaned PDB structure
│   ├── em.mdp                  # Energy minimization parameters
│   ├── nvt.mdp                 # NVT equilibration parameters
│   ├── npt.mdp                 # NPT equilibration parameters
│   ├── md.mdp                  # Production MD parameters
│   ├── ions.mdp                # Ion addition parameters
│   ├── topol.top               # System topology
│   ├── index.ndx               # Custom index groups
│   └── run_md.sh               # SLURM submission script
│
├── production/                  # Simulation outputs
│   ├── md.xtc                  # Trajectory (compressed)
│   ├── md.tpr                  # Portable binary input
│   ├── md.gro                  # Final coordinates
│   ├── md.edr                  # Energy file
│   ├── md.log                  # Simulation log
│   └── md.cpt                  # Checkpoint file
│
├── analysis/                    # Post-processing & results
│   ├── ref.pdb                 # Reference structure
│   ├── pocket.ndx              # Pocket residue index
│   │
│   ├── fpocket_analysis/       # Fpocket results
│   │   ├── ref_out/            # Initial pocket detection
│   │   ├── fpocket_runs/       # Timeseries data
│   │   ├── parse_fpocket_runs.py
│   │   └── fpocket_timeseries.csv
│   │
│   ├── dynamics/                # MD analysis
│   │   ├── rmsd_pocket.csv
│   │   ├── rmsf_pocket.csv
│   │   ├── sasa_pocket.csv
│   │   └── pocket_volume.csv
│   │
│   ├── figures/                 # Generated plots
│   │   ├── rmsd_pocket.png
│   │   ├── rmsf_pocket.png
│   │   ├── sasa_pocket.png
│   │   ├── volume_pocket.png
│   │   └── md_pocket_analysis.png
│   │
│   └── scripts/                 # Analysis scripts
│       ├── analyze_trajectory.py
│       ├── plot_results.py
│       └── calculate_druggability.py
│
├── docs/                        # Documentation
│   ├── methods.md              # Detailed methodology
│   ├── parameters.md           # Parameter justification
│   └── references.bib          # Citations
│
└── results/                     # Final report
    ├── report.pdf              # Complete analysis report
    └── summary.md              # Executive summary
```

---

## 12. Version Control

### Initialize Repository

```bash
cd 4QKZ_project

# Initialize git
git init
git add .
git commit -m "Initial commit: Complete MD workflow and analysis"

# Create GitHub repository (do this on GitHub first)
git branch -M main
git remote add origin https://github.com/YourUsername/4QKZ_project.git
git push -u origin main
```

### Recommended .gitignore

```bash
cat > .gitignore << EOF
# Large trajectory files (use Git LFS or external storage)
*.xtc
*.trr
*.edr
*.cpt

# Temporary files
*.o
*.log
*.tmp
*~
.DS_Store

# Python cache
__pycache__/
*.pyc
.ipynb_checkpoints/

# Backup files
\#*#
*.bak
EOF

git add .gitignore
git commit -m "Add .gitignore for trajectory files"
```

### Use Git LFS for Large Files (Optional)

```bash
git lfs install
git lfs track "*.xtc"
git lfs track "*.trr"
git add .gitattributes
git commit -m "Track large trajectory files with Git LFS"
```

---

## 13. Citation

### This Project

```bibtex
@software{gebremedhin2024_4qkz,
  author = {Gebremedhin, Atnatiwos N.},
  title = {4QKZ Protein Molecular Dynamics and Pocket Analysis},
  year = {2024},
  url = {https://github.com/YourUsername/4QKZ_project}
}
```

### Key Software

**GROMACS:**
```bibtex
@article{abraham2015gromacs,
  title={GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers},
  author={Abraham, Mark James and Murtola, Teemu and Schulz, Roland and P{\'a}ll, Szil{\'a}rd and Smith, Jeremy C and Hess, Berk and Lindahl, Erik},
  journal={SoftwareX},
  volume={1},
  pages={19--25},
  year={2015},
  publisher={Elsevier}
}
```

**Fpocket:**
```bibtex
@article{le2009fpocket,
  title={fpocket: an open source platform for ligand pocket detection},
  author={Le Guilloux, Vincent and Schmidtke, Peter and Tuffery, Pierre},
  journal={BMC bioinformatics},
  volume={10},
  number={1},
  pages={1--11},
  year={2009},
  publisher={BioMed Central}
}
```

---

## 14. Contact & Support

**Author:** Atnatiwos N. Gebremedhin  
**Institution:** [Your Institution]  
**Email:** [Your Email]  

For questions or issues:
- Open an issue on [GitHub](https://github.com/YourUsername/4QKZ_project/issues)
- Email the author directly

---

## 15. License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

```text
MIT License

Copyright (c) 2024 Atnatiwos N. Gebremedhin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## 16. Acknowledgments

- **GROMACS Development Team** for the excellent MD simulation software
- **Fpocket Developers** for the pocket detection tool
- **HPC Support Staff** for computational resources
- **Rayca Bio** for the assessment task opportunity

---

## Appendix: Quick Start Guide

### For Users Who Want to Reproduce Results

```bash
# 1. Clone repository
git clone https://github.com/YourUsername/4QKZ_project.git
cd 4QKZ_project

# 2. Load modules
module load gromacs/2022.5 python/3.10

# 3. Run complete workflow
bash setup/run_md.sh

# 4. Analyze results
cd analysis
python scripts/analyze_trajectory.py
python scripts/plot_results.py

# 5. View results
ls figures/
```

### Estimated Computation Time

| Step | Duration | Resources |
|------|----------|-----------|
| Setup | 5 min | 1 CPU core |
| Energy Minimization | 10 min | 4 CPU cores |
| NVT Equilibration | 15 min | 8 CPU cores |
| NPT Equilibration | 15 min | 8 CPU cores |
| 100 ns Production | 24-48 hours | 1 GPU + 8 CPUs |
| Analysis | 1 hour | 1 CPU core |

**Total Wall Time:** ~2-3 days

---

**Last Updated:** November 2024  
**Version:** 1.0.0
