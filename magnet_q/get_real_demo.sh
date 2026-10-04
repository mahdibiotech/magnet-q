#!/usr/bin/env bash
# Jeu réel minuscule (tutoriel Galaxy, Zenodo 17661262) : contigs MEGAHIT + bins CONCOCT.
# Sert à valider la chaîne complète ; trop petit pour juger la qualité des bins.
set -euo pipefail
mkdir -p real/bins && cd real
Z=https://zenodo.org/record/17661262/files
wget -nc -q $Z/MEGAHIT_contigs.fasta
for i in 0 1 2 3 4 5 6 7 8 9; do wget -nc -q -P bins $Z/concoct_$i.fasta; done
ls -la . bins
