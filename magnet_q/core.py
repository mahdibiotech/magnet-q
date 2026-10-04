"""MAGNET-Q : logique métier (sans Dash) pour pouvoir la tester et la benchmarker."""
import json
import os
import sqlite3
import time

import numpy as np
import pandas as pd

N_MARKERS = 107  # taille du jeu de gènes marqueurs bactériens mono-copie (type CheckM)
UNBINNED = "unbinned"


def simulate(seed=7, n_genomes=10):
    """Jeu de démo synthétique avec vérité terrain (contigs, marqueurs, binning bruité)."""
    rng = np.random.default_rng(seed)
    gc0 = rng.uniform(0.38, 0.66, n_genomes)
    cov0 = rng.uniform(1.5, 3.5, n_genomes)  # ln(couverture)
    contigs, markers, cid = [], [], 0
    for g in range(n_genomes):
        n = int(rng.integers(60, 220))
        present = rng.random(N_MARKERS) < rng.uniform(0.93, 0.995)
        owner = rng.integers(0, n, N_MARKERS)
        for i in range(n):
            name = f"ctg_{cid:05d}"
            contigs.append((name, g, float(np.clip(rng.normal(gc0[g], 0.025), 0.2, 0.8)),
                            float(rng.normal(cov0[g], 0.3)), int(rng.lognormal(9.5, 0.8))))
            markers += [(name, int(m)) for m in np.where((owner == i) & present)[0]]
            cid += 1
    df = pd.DataFrame(contigs, columns=["contig_id", "true_genome", "gc", "cov", "length"])
    # binning bruité : 12 % des contigs vont dans le bin du génome le plus proche
    z = (df[["gc", "cov"]] - df[["gc", "cov"]].mean()) / df[["gc", "cov"]].std()
    cent = z.groupby(df.true_genome).mean().to_numpy()
    bins = []
    for zi, g in zip(z.to_numpy(), df.true_genome):
        if rng.random() < 0.12:
            d = np.linalg.norm(cent - zi, axis=1)
            d[g] = np.inf
            g = int(d.argmin())
        bins.append(f"bin_{g:03d}")
    df["bin_id"] = bins
    r2 = np.random.default_rng(seed + 1)
    ctr = r2.normal(0, 3, (n_genomes, 2))
    df[["pc1", "pc2"]] = ctr[df.true_genome.to_numpy()] + r2.normal(0, 0.8, (len(df), 2))
    return df, pd.DataFrame(markers, columns=["contig_id", "marker"])


def label(comp, cont):
    if comp > 90 and cont < 5:
        return "HQ"
    if comp >= 50 and cont < 10:
        return "MQ"
    return "LQ"


class Curator:
    def __init__(self, db="data/magnet_q.db"):
        os.makedirs(os.path.dirname(db) or ".", exist_ok=True)
        new = not os.path.exists(db)
        self.cx = sqlite3.connect(db, check_same_thread=False)
        if new:
            c, m = simulate()
            c.to_sql("contigs", self.cx, index=False)
            m.to_sql("markers", self.cx, index=False)
            self.cx.execute("CREATE TABLE actions(id INTEGER PRIMARY KEY, ts TEXT, kind TEXT, ids TEXT,"
                            " src TEXT, dst TEXT, undone INT DEFAULT 0)")
            self.cx.commit()
        self.n_markers, self.assembly = N_MARKERS, ""
        try:
            meta = dict(self.cx.execute("SELECT k, v FROM meta").fetchall())
            self.n_markers, self.assembly = int(meta.get("n_markers", N_MARKERS)), meta.get("assembly", "")
        except sqlite3.OperationalError:
            pass
        self.contigs = pd.read_sql("SELECT * FROM contigs", self.cx).set_index("contig_id")
        self.markers = pd.read_sql("SELECT * FROM markers", self.cx)
        for a in self.cx.execute("SELECT kind, ids, src, dst FROM actions WHERE undone=0 ORDER BY id"):
            self.contigs.loc[json.loads(a[1]), "bin_id"] = a[3]

    # ---------- lecture ----------
    def bins(self):
        return sorted(b for b in self.contigs.bin_id.unique() if b != UNBINNED)

    def members(self, bin_id):
        return self.contigs[self.contigs.bin_id == bin_id]

    def quality(self, bin_id, exclude=()):
        ids = self.members(bin_id).index.difference(list(exclude))
        counts = self.markers[self.markers.contig_id.isin(ids)].marker.value_counts()
        comp = 100 * len(counts) / self.n_markers
        cont = 100 * (counts.sum() - len(counts)) / self.n_markers
        return {"completeness": comp, "contamination": cont, "label": label(comp, cont),
                "size": int(self.contigs.loc[ids, "length"].sum()), "n": len(ids)}

    def suspects(self, bin_id, top=8):
        """Contigs classés par gain de contamination moins perte de complétude si retirés."""
        ids = self.members(bin_id).index
        m = self.markers[self.markers.contig_id.isin(ids)].copy()
        counts = m.marker.value_counts()
        m["dup"] = m.marker.map(counts) > 1
        g = m.groupby("contig_id").agg(dup=("dup", "sum"), tot=("marker", "size"))
        g["d_cont"] = -100 * g.dup / self.n_markers
        g["d_comp"] = -100 * (g.tot - g.dup) / self.n_markers
        g["score"] = g.dup * 1.0 - (g.tot - g.dup) * 1.5
        return g[g.dup > 0].sort_values("score", ascending=False).head(top)

    # ---------- écriture + journal ----------
    def _apply(self, ids, src, dst):
        self.contigs.loc[ids, "bin_id"] = dst

    def move(self, ids, src, dst, kind):
        ids = [i for i in ids if i in self.contigs.index and self.contigs.at[i, "bin_id"] == src]
        if not ids or src == dst:
            return 0
        self.cx.execute("DELETE FROM actions WHERE undone=1")
        self.cx.execute("INSERT INTO actions(ts,kind,ids,src,dst) VALUES (?,?,?,?,?)",
                        (time.strftime("%H:%M:%S"), kind, json.dumps(ids), src, dst))
        self.cx.commit()
        self._apply(ids, src, dst)
        return len(ids)

    def new_bin_name(self):
        nums = [int(b.split("_")[1]) for b in self.bins() if b.split("_")[1].isdigit()]
        return f"bin_{max(nums, default=0) + 1:03d}"

    def undo(self):
        r = self.cx.execute("SELECT id, ids, src, dst FROM actions WHERE undone=0 ORDER BY id DESC LIMIT 1").fetchone()
        if r:
            self._apply(json.loads(r[1]), r[3], r[2])
            self.cx.execute("UPDATE actions SET undone=1 WHERE id=?", (r[0],))
            self.cx.commit()
        return bool(r)

    def redo(self):
        r = self.cx.execute("SELECT id, ids, src, dst FROM actions WHERE undone=1 ORDER BY id LIMIT 1").fetchone()
        if r:
            self._apply(json.loads(r[1]), r[2], r[3])
            self.cx.execute("UPDATE actions SET undone=0 WHERE id=?", (r[0],))
            self.cx.commit()
        return bool(r)

    def journal(self, n=6):
        rows = self.cx.execute("SELECT ts, kind, ids, src, dst FROM actions WHERE undone=0 ORDER BY id DESC LIMIT ?", (n,))
        return [f"{t} · {k} {len(json.loads(i))} contigs : {s} → {d}" for t, k, i, s, d in rows]

    def export(self, path="data/curated_bins.tsv"):
        self.contigs[["bin_id"]].to_csv(path, sep="\t")
        return path
