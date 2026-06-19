# LLM-Enhanced ChemBERTa Screening Reveals Repositionable Caspase-4 Inhibitors: Insights from Molecular Docking and MD Simulations

## Abstract

This study presents an integrated ligand-based virtual screening pipeline to identify potential caspase-4 inhibitors from the DrugBank database by combining docking scores, SMILES-derived ChemBERTa embeddings, and key physicochemical descriptors. A random forest model trained on experimental pIC50 values demonstrated strong predictive performance, achieving ~95% accuracy and ~3.5-fold early enrichment of active compounds. Top candidates included DB00519, DB01068, DB06202, DB08882, and DB05316, showing high predicted activity and favorable docking scores. Molecular dynamics simulations further confirmed the stability and consistent energetic behavior of the top protein–ligand complexes. Overall, this workflow provides an efficient and robust strategy for prioritizing repurposable DrugBank compounds as potential caspase-4 inhibitors.

---

## Pipeline Overview

```
DrugBank compounds
       ↓
AutoDock Vina (molecular docking)
       ↓
Feature extraction:
  • ChemBERTa SMILES embeddings
  • RDKit physicochemical descriptors
       ↓
Random Forest Classifier (active vs. inactive)
       ↓
Random Forest Regressor (pIC50 prediction)
       ↓
Top-10 candidate selection
       ↓
ADMET + SA score + PK-PD simulation
       ↓
MD simulation validation
```

---

## Repository Structure

```
├── Caspase4_github.ipynb                        # Main analysis notebook
├── Screening energy VINA_caspase-4.csv          # AutoDock Vina docking scores (input)
├── screening_with_drugbank_ids.csv              # Processed screening results
├── sascorer.py                                  # Synthetic accessibility scorer
├── fpscores.pkl.gz                              # Fragment scores for SA calculation
├── collect_scores_sdf.py                        # Script to collect Vina scores from SDF
├── make_ligands_pdbqt.py                        # Script to prepare ligands for docking
├── run_vina_sdf.sh                              # AutoDock Vina docking shell script
├── config.txt                                   # Vina docking configuration
├── top10_smiles.csv                             # SMILES of top-10 predicted compounds
├── top10_molecular_properties.csv               # Molecular properties of top-10 hits
├── top10_molecular_properties_with_PK_like_cols.csv  # Extended PK-like properties
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
- Loading and processing docking scores
- Fetching pIC50 values from ChEMBL
- Generating ChemBERTa embeddings
- Training the Random Forest classifier and regressor
- Selecting and analyzing top-10 hits
- PK-PD simulations and drug-likeness analysis

---

## Key Results

| Metric | Value |
|--------|-------|
| Classifier Accuracy | ~95% |
| ROC-AUC | ~0.72 |
| Enrichment Factor (top 10%) | ~3.5× |
| Top candidates | DB00519, DB01068, DB06202, DB08882, DB05316 |

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

> Manuscript under review. Citation will be added upon publication.

---

## License

This project is for academic research purposes. DrugBank data usage is subject to [DrugBank's terms of use](https://go.drugbank.com/legal/terms_of_use).
