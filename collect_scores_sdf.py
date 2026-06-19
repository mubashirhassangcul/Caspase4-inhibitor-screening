import os, glob, re
import pandas as pd

INPUT_CSV = "top10_smiles.csv"      # has drugbank_id column
LOG_DIR   = "docking_out_sdf"

df = pd.read_csv(INPUT_CSV)

def parse_vina_log(logfile):
    scores = []
    with open(logfile) as f:
        for line in f:
            if re.match(r"\s*\d+\s+[-0-9.]+\s", line):
                parts = line.split()
                mode = int(parts[0])
                score = float(parts[1])
                scores.append((mode, score))
    return scores[0][1] if scores else None

rows = []
for log in glob.glob(os.path.join(LOG_DIR, "*.log")):
    base = os.path.basename(log).replace(".log", "")
    dbid = base                      # because logs are named DBxxxx.log
    best = parse_vina_log(log)
    rows.append({"drugbank_id": dbid, "vina_score": best})

dock_df = pd.DataFrame(rows)
merged = df.merge(dock_df, on="drugbank_id", how="left")
merged = merged.sort_values("vina_score")   # more negative = better
merged.to_csv("top10_with_vina_sdf.csv", index=False)
print(merged[["drugbank_id", "vina_score"]])
