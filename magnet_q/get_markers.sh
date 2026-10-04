#!/usr/bin/env bash
set -u
mkdir -p real && cd real
BASE=https://busco-data.ezlab.org/v5/data
V=$(curl -fsSL "$BASE/file_versions.tsv" 2>/dev/null | awk -F'\t' '$1=="bacteria_odb10"{print $2; exit}')
FILE=""
for d in $V 2024-01-08 2020-03-06; do
  if curl -fsL --range 0-99 -o /dev/null "$BASE/lineages/bacteria_odb10.$d.tar.gz"; then
    FILE="bacteria_odb10.$d.tar.gz"; break
  fi
done
[ -n "$FILE" ] || { echo "Aucun fichier bacteria_odb10 trouvé (réseau ou site BUSCO indisponible)."; exit 1; }
echo "Téléchargement de $FILE"
curl -fSL -o "$FILE" "$BASE/lineages/$FILE"
tar -xzf "$FILE"
H=$(find . -type d -name hmms | head -1)
[ -n "$H" ] || { echo "Dossier hmms introuvable dans l'archive :"; ls; exit 1; }
cat "$H"/*.hmm > bacteria_markers.hmm
echo "Marqueurs : $(grep -c '^NAME' bacteria_markers.hmm)"
