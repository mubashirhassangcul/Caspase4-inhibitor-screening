import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Load SMILES (drugbank_id, smiles)
df = pd.read_csv("top10_smiles.csv")
df = df[["drugbank_id", "smiles"]].drop_duplicates()

os.makedirs("ligands_pdbqt", exist_ok=True)

def smiles_to_pdbqt(dbid, smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"Failed SMILES for {dbid}")
        return
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xf00d) != 0:
        print(f"Embed failed for {dbid}")
        return
    AllChem.UFFOptimizeMolecule(mol)

    pdb_block = Chem.MolToPDBBlock(mol)
    out_path = os.path.join("ligands_pdbqt", f"{dbid}.pdbqt")
    with open(out_path, "w") as f:
        for line in pdb_block.splitlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                f.write(line.replace("HETATM", "ATOM  ") + "\n")

for _, row in df.iterrows():
    smiles_to_pdbqt(str(row["drugbank_id"]), row["smiles"])
