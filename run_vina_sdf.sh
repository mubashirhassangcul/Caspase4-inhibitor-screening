#!/bin/bash

VINA="./vina_1.2.7_mac_aarch64"

mkdir -p docking_out_sdf

for f in ligands_sdf/*.sdf
do
    base=$(basename "$f" .sdf)
    echo "Docking $base"
    "$VINA" --receptor caspase4_receptor.pdbqt \
            --ligand "$f" \
            --config config.txt \
            --out docking_out_sdf/"${base}_out.pdbqt"
done
