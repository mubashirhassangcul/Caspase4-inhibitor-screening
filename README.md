# AI Accelerated Chemical Screening Integrates ChemBERTa to Identify Repositionable CASP4 Inhibitors via MD Simulations and MM/PBSA Analysis

## Abstract

This study employed an integrated ligand-based virtual screening pipeline to identify potential CASP4 inhibitors from the DrugBank database, leveraging docking-score prioritization, SMILES-derived ChemBERTa embeddings, and key physicochemical descriptors. Building on this foundation, the workflow incorporated virtual screening, cheminformatics modeling, PK–PD evaluation, molecular docking, molecular dynamics simulations, and MM/PBSA free-energy analysis to systematically prioritize repurposed DrugBank compounds for CASP4-targeted Alzheimer's disease therapy. A Random Forest classifier trained on the hybrid ChemBERTa–physicochemical feature set distinguished active from inactive compounds with ~95% accuracy (ROC-AUC = 0.73) and achieved a ~3.5-fold enrichment of active compounds among the top-ranked hits, while a companion Random Forest regressor, trained on the experimental pIC50 values of active compounds, ranked candidates by predicted potency. In addition, integrated cheminformatics modeling, PK-PD analysis, and molecular docking further narrowed the selection to the top five candidate compounds: DB00519, DB01068, DB06202, DB08882, and DB05316. Finally, MD simulations and MM/PBSA calculations established clear thermodynamic support for DB05316 and DB00519, whose binding free energies (−21.4 and −20.9 kcal/mol, respectively) exceeded even the reference compound donepezil, highlighting them as the most promising CASP4 inhibitors. Overall, this workflow provides an efficient and robust strategy for prioritizing repositioned DrugBank compounds as potential CASP4 inhibitors against AD.

---

## Pipeline Overview

```
DrugBank Drug-Lib compounds
       ↓
Structure-based virtual screening (VSTH/Tianhe-2 vs. CASP4, PDB 6NRY) — ~1,739 initial hits
       ↓
AI-based filtering of screening hits
       ↓
ChEMBL (fetch experimental pIC50 labels; active if pIC50 ≥ 8.2)
       ↓
Feature extraction:
  • ChemBERTa SMILES embeddings (768-dim)
  • RDKit physicochemical descriptors (6)
       ↓
Random Forest Classifier (active vs. inactive; 500 trees, 774-dim hybrid features)
       ↓
Random Forest Regressor (pIC50 prediction on actives)
       ↓
Top-10 candidate selection
       ↓
ADMET (ADMET Tab 3.0 / DMPNN) + SA score + drug-likeness
       ↓
PK-PD simulations + dose-response curves (SciPy / PySB)
       ↓
Binding-site ID (PrankWeb) + AutoDock Vina docking vs. CASP4
       ↓
Top-5 candidates + donepezil reference → 100-ns MD in triplicate (GROMACS, CHARMM36)
       ↓
MM/PBSA free-energy ranking (gmx_MMPBSA)
       ↓
Benchmarking vs. donepezil
```

---

## Repository Structure

```
├── Caspase4_github.ipynb                        # Main analysis notebook
├── Screening energy VINA_caspase-4.csv          # AutoDock Vina docking scores (input)
├── screening_with_drugbank_ids.csv              # Processed screening results
├── chemberta_input_features.csv / .xlsx         # ChemBERTa + RDKit hybrid feature set (774-dim)
├── qsar_data_with_descriptors.csv               # QSAR descriptor table for screened compounds
├── sascorer.py                                  # Synthetic accessibility scorer
├── fpscores.pkl.gz                               # Fragment scores for SA calculation
├── collect_scores_sdf.py                        # Script to collect Vina scores from SDF
├── make_ligands_pdbqt.py                        # Script to prepare ligands for docking
├── run_vina_sdf.sh                              # AutoDock Vina docking shell script
├── config.txt                                   # Vina docking configuration
├── top10_smiles.csv                             # SMILES of top-10 predicted compounds
├── top10_molecular_properties.csv               # Molecular properties of top-10 hits
├── top10_molecular_properties_with_PK_like_cols.csv  # Extended PK-like properties
├── Table2_structures/                           # 2D structures of the 10 screened DrugBank compounds
│   ├── DB00439_Cerivastatin.png
│   ├── DB00519_Trandolapril.png
│   ├── DB01068_Clonazepam.png
│   ├── DB01544_Flunitrazepam.png
│   ├── DB05316_Pimavanserin.png
│   ├── DB06202_Lasofoxifene.png
│   ├── DB06203_Alogliptin.png
│   ├── DB08882_Linagliptin.png
│   ├── DB08897_Aclidinium.png
│   └── DB09477_Enalaprilat.png
└── New graphs/                                  # MD simulation analysis figures
    ├── RMSD_triplicates.png
    ├── RMSF_triplicates.png
    ├── SASA_triplicates.png
    ├── Gyration_triplicates.png
    ├── ProtLigDist_triplicates.png
    ├── TotalEnergy_triplicates.png
    ├── plot_all.py
    └── plot_rmsd.py
```

---

## Requirements

```bash
pip install transformers torch scikit-learn pandas numpy rdkit-pypi matplotlib seaborn
pip install chembl_downloader chembl-webresource-client
```

> **Note:** DrugBank data (`drugbank.xml`) is required to reproduce the full screening but is not included due to licensing restrictions. Register and download it from [drugbank.ca](https://go.drugbank.com/).

The molecular dynamics and MM/PBSA stages additionally require external, non-pip tools that are not bundled in this repo: [GROMACS](https://www.gromacs.org) 2019.3 (with the CHARMM36 force field), [CHARMM-GUI](https://www.charmm-gui.org) for system setup, and [gmx_MMPBSA](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA) v1.6.3 for free-energy post-processing.

---

## Usage

1. Clone the repository:
```bash
git clone https://github.com/mubashirhassangcul/Caspase4-inhibitor-screening.git
cd Caspase4-inhibitor-screening
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Open and run the notebook:
```bash
jupyter notebook Caspase4_github.ipynb
```

Run cells sequentially. The notebook covers:

**Data & Feature Extraction**
- Loading and processing AutoDock Vina docking scores
- Extracting DrugBank IDs and molecular descriptors
- Fetching experimental pIC50 values from ChEMBL (active if pIC50 ≥ 8.2)
- Generating ChemBERTa SMILES embeddings (768-dim, frozen pretrained encoder)

**Machine Learning**
- Training Random Forest classifier (active vs. inactive compounds; hybrid 774-dim feature set)
- Training Random Forest regressor (pIC50 prediction on actives)
- Ranking and selecting top-10 candidate compounds

**Model Evaluation Plots**
- ROC curve (AUC = 0.73)
- Enrichment curve (top 10% captures ~3.5× more actives than random)
- Actual vs. predicted pIC50 error plot for top-10 hits
- t-SNE 2D and 3D chemical space visualization
- Feature distribution plots (MW, LogP, HBA, HBD, rotatable bonds, aromatic rings) for active vs. inactive compounds

**Drug-likeness & ADMET Analysis**
- Molecular property table (MW, LogP, HBD, HBA, TPSA, RotBonds)
- Synthetic Accessibility (SA) score for top-10 hits
- ADMET Tab 3.0 (DMPNN-based) pharmacokinetic and toxicity predictions

**Pharmacology Plots**
- Radar plots for individual compound PK profile
- Radar plots for all 10 compounds (5×2 grid)
- Simulated dose–response curves for top-10 hits
- Integrated PK-PD simulations with variable half-lives per compound
- Dual-axis concentration and effect vs. time plots

**Molecular Docking & MM/PBSA Free Energy**
- Binding-site identification on CASP4 (PDB 6NRY) via PrankWeb
- AutoDock Vina docking of all screened compounds against the CASP4 catalytic pocket
- 100-ns triplicate MD simulations (GROMACS, CHARMM36) of the top-5 complexes plus donepezil
- MM/PBSA binding free-energy ranking (gmx_MMPBSA) across the full production trajectories

**Benchmarking**
- Docking-energy and MM/PBSA comparison against the reference Alzheimer's drug donepezil

---

## Key Results

| Metric | Value |
|--------|-------|
| Classifier Accuracy | ~95% |
| ROC-AUC | 0.73 |
| Enrichment Factor (top 10%) | ~3.5× |
| Top candidates | DB00519, DB01068, DB06202, DB08882, DB05316 |
| Best docking energy | DB01068, −9.75 kcal/mol (vs. donepezil −7.9 kcal/mol) |
| Best MM/PBSA ΔG_bind | DB05316, −21.4 kcal/mol; DB00519, −20.9 kcal/mol (vs. donepezil −17.1 kcal/mol) |

MM/PBSA thermodynamic ranking across the full 100-ns trajectories: **DB05316 > DB00519 > DB06202 > donepezil > DB08882 > DB01068** — DB05316 and DB00519 are the two candidates that outperform the donepezil reference.

---

## MD Simulation Figures

Molecular dynamics simulations were performed in triplicate for the top protein–ligand complexes. Results are in `New graphs/`:

| Figure | Description |
|--------|-------------|
| `RMSD_triplicates.png` | Backbone RMSD over simulation time |
| `RMSF_triplicates.png` | Per-residue flexibility |
| `SASA_triplicates.png` | Solvent accessible surface area |
| `Gyration_triplicates.png` | Radius of gyration |
| `ProtLigDist_triplicates.png` | Protein–ligand distance |
| `TotalEnergy_triplicates.png` | System total energy |

---

## Citation

> Manuscript is accepted. Citation will be added upon publication.

---

## License

This project is for academic research purposes. DrugBank data usage is subject to [DrugBank's terms of use](https://go.drugbank.com/legal/terms_of_use).
