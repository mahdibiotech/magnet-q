#!/usr/bin/env bash
set -eu
mkdir -p coffee/bins && cd coffee
Z="https://zenodo.org/records/7845138/files"
curl -fSL -o bins.zip "$Z/26_%20MetaBAT2%20on%20data%20ERR2231567_%20Bins.zip?download=1"
curl -fSL -o checkm_stats.txt "$Z/CheckM_lineage_wf_on_data_ERR2231567__Bin_statistics.txt?download=1"
python3 -m zipfile -e bins.zip bins_raw
find bins_raw -type f | head -20
find bins_raw -type f | while read -r f; do cp "$f" "bins/$(basename "$f" | tr ' ' '_')"; done
cat bins/* > assembly.fa
echo "Bins : $(ls bins | wc -l) ; contigs : $(grep -c '^>' assembly.fa)"
head -c 300 assembly.fa | head -3
