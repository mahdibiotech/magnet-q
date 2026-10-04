# MAGNET-Q

Prototype indépendant (Python/Dash) pour explorer la curation manuelle de MAGs : estimation rapide de la complétude et de la contamination d'un bin à partir de gènes marqueurs, effet estimé d'un retrait de contigs avant validation, journal annulable (SQLite), export.
Je ne connais pas le code de MAGNET (LABGeM/Genoscope) : ce dépôt n'est ni un remplacement ni une copie, seulement une exploration.

## Données
Jeu public réel : bins MetaBAT2 de l'échantillon ERR2231567 (fermentation du café), Zenodo, enregistrement 7845138 (licence et auteurs : voir la page). L'assemblage est reconstitué en concaténant les contigs des bins, donc les contigs non classés sont absents. Pas de fichier de couverture : seule la vue « Composition » (tétranucléotides, PCA) est utilisée.

## Reproduire
```bash
conda install -c conda-forge -c bioconda prodigal hmmer      # + pip install dash plotly pandas numpy
cd magnet_q
bash get_markers.sh                      # marqueurs BUSCO bacteria_odb10 (124)
bash get_coffee_demo.sh                  # jeu café + tableau CheckM de référence
python prepare_data.py --assembly coffee/assembly.fa --bins coffee/bins \
    --hmm real/bacteria_markers.hmm --cutoffs real/bacteria_odb10/scores_cutoff \
    --cutoff-factor 0.3 --threads 4
python app.py                            # http://127.0.0.1:8050 ; curer le bin 5
python report.py 5                       # bilan avant/après de l'outil
python validate_checkm2.py --checkm-stats coffee/checkm_stats.txt --tag avant_checkm
```
Pour CheckM2 : environnement séparé (`conda create -n checkm2 --override-channels -c conda-forge -c bioconda checkm2`), puis `checkm2 database --download`, `python validate_checkm2.py --tag apres` pour les bins curés, et `checkm2 predict -i coffee/bins -x fasta -o data/checkm2_avant` pour les bins d'origine.

## Résultats (7 bins, un seul échantillon)
**1. Estimation rapide contre CheckM, bins d'origine** (écart absolu moyen, en points) :

| Facteur appliqué aux seuils de score BUSCO | Complétude | Contamination |
|---|---|---|
| 0 (aucun filtre) | 4,1 | 4,7 |
| 0,3 (réglage retenu) | 3,7 | 3,4 |
| 0,5 | 3,1 | 3,9 |
| 0,7 | 3,7 | 4,6 |
| 1 (seuils BUSCO) | 4,5 | 5,1 |

Le facteur 0,3 a été choisi en regardant l'écart avec CheckM sur ces mêmes 7 bins : le résultat est optimiste et non validé sur un autre jeu. Deux bins (1 et 5) concentrent presque tout l'écart.

**2. Estimation rapide contre CheckM2** (bins après curation du bin 5) : écart moyen de 12,3 points (complétude) et 6,0 (contamination) sur les 7 bins ; 5,0 et 2,5 si l'on exclut le bin 5.

**3. Curation du bin 5** (4 142 contigs, 12,85 Mb → 3 195 contigs, 9,83 Mb ; 947 contigs retirés, 23 %) :

| | Complétude | Contamination |
|---|---|---|
| CheckM2 avant | 100 | 243,5 |
| CheckM2 après | 100 | 56,9 |
| Estimation de l'outil, avant → après | 56,5 → 43,5 | 53,2 → 29,8 |

Les six autres bins sont identiques avant et après dans CheckM2.

## Limites
- Un seul bin curé, un seul curateur (l'auteur de l'outil), pas de témoin sans le panneau de qualité : on ne peut pas dire que l'outil a aidé plus qu'une sélection à l'œil sur la composition.
- Le bin 5 reste très contaminé après curation (56,9 %) : il faudrait le séparer en plusieurs génomes.
- L'estimation de l'outil sous-estime la contamination du bin 5 (29,8 % contre 56,9 % selon CheckM2) ; il donne la direction du changement, pas son ampleur.
- Sur ce bin, CheckM (59,6 %) et CheckM2 (243,5 %) divergent fortement avant curation, et CheckM2 affiche 100 % de complétude avant et après : le gain de complétude ou son coût n'est pas mesurable avec CheckM2 ici.
- Un marqueur en double signale un conflit entre deux contigs, pas lequel est l'intrus.
- Pas de couverture, donc la vue GC × couverture est inutilisable sur ce jeu.

## Fichiers
`core.py` (logique), `app.py` (interface), `prepare_data.py` (import), `validate_checkm2.py` (comparaison), `report.py` (bilan avant/après), `benchmark.py` (jeu simulé de démonstration uniquement ; ses résultats ne sont pas ceux ci-dessus).
