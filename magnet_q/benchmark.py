"""Benchmark sur données simulées avec vérité terrain.

Mesure : (1) latence du retour qualité ; (2) valeur du classement des contigs suspects,
en simulant un curateur qui retire le contig suspect n°1 tant que contamination >= 5 %.
Limite assumée : l'estimateur et la simulation partagent le même jeu de marqueurs ;
la fidélité vis-à-vis de CheckM2 doit être mesurée sur données réelles (voir README).
"""
import os
import tempfile
import time

import pandas as pd

from core import Curator, label


def true_foreign(cur, bin_id):
    m = cur.members(bin_id)
    major = m.true_genome.mode()[0]
    return set(m.index[m.true_genome != major])


def run(db=None):
    cur = Curator(db or os.path.join(tempfile.mkdtemp(), "bench.db"))
    t = time.perf_counter()
    for b in cur.bins():
        cur.quality(b)
    latency_ms = 1000 * (time.perf_counter() - t) / len(cur.bins())

    rows = []
    for b in cur.bins():
        before = cur.quality(b)
        foreign, removed, steps = true_foreign(cur, b), [], 0
        while cur.quality(b)["contamination"] >= 5 and steps < 60:
            s = cur.suspects(b, 1)
            if s.empty:
                break
            cid = s.index[0]
            cur.move([cid], b, "unbinned", "retrait")
            removed.append(cid)
            steps += 1
        after = cur.quality(b)
        rows.append({"bin": b, "cont_avant": before["contamination"], "cont_apres": after["contamination"],
                     "comp_avant": before["completeness"], "comp_apres": after["completeness"],
                     "label_avant": before["label"], "label_apres": after["label"],
                     "retires": len(removed), "retires_reellement_etrangers": len(set(removed) & foreign)})
    df = pd.DataFrame(rows).round(1)
    os.makedirs("data", exist_ok=True)
    df.to_csv("data/benchmark.csv", index=False)
    print(df.to_string(index=False))
    n_hq0, n_hq1 = (df.label_avant == "HQ").sum(), (df.label_apres == "HQ").sum()
    print(f"\nLatence estimation qualité : {latency_ms:.1f} ms / bin")
    print(f"Bins HQ : {n_hq0} -> {n_hq1} sur {len(df)}")
    print(f"Perte moyenne de complétude : {(df.comp_avant - df.comp_apres).mean():.1f} points")
    r = df.retires.sum()
    print(f"Précision des retraits (contigs réellement étrangers) : "
          f"{df.retires_reellement_etrangers.sum()}/{r}" if r else "Aucun retrait")


if __name__ == "__main__":
    run()
