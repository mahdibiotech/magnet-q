"""Construit la base MAGNET-Q à partir de vraies données.

Entrées :
  --assembly contigs.fa            assemblage (FASTA)
  --bins DIR | bins.tsv            dossier de FASTA (un fichier par bin) ou TSV contig<TAB>bin
  --depth depth.txt                sortie de jgi_summarize_bam_contig_depths (MetaBAT2) : contigName, contigLen, totalAvgDepth
  Marqueurs, au choix :
    --hmm markers.hmm              HMM de gènes mono-copie (ex. Bacteria_71 d'anvi'o) ; nécessite prodigal + hmmsearch
    --markers-tsv m.tsv            table contig<TAB>marker déjà calculée (avec --n-markers)
Sortie : data/magnet_q.db  (puis : python app.py)
"""
import argparse
import itertools
import os
import re
import sqlite3
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

COMP = str.maketrans("ACGT", "TGCA")
BASES = "ACGT"


def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq).upper()
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq).upper()


def canonical_kmers(k=4):
    idx = {}
    for t in map("".join, itertools.product(BASES, repeat=k)):
        rc = t.translate(COMP)[::-1]
        idx.setdefault(min(t, rc), len(idx))
    return idx


def tnf_vector(seq, idx, k=4):
    v = np.zeros(len(idx))
    for i in range(len(seq) - k + 1):
        t = seq[i:i + k]
        if set(t) <= set(BASES):
            rc = t.translate(COMP)[::-1]
            v[idx[min(t, rc)]] += 1
    return v


def composition(assembly, min_len):
    idx = canonical_kmers()
    rows, vecs = [], []
    for name, seq in read_fasta(assembly):
        if len(seq) < min_len:
            continue
        gc = (seq.count("G") + seq.count("C")) / len(seq)
        rows.append((name, len(seq), gc))
        v = tnf_vector(seq, idx) + 1.0
        f = v / v.sum()
        vecs.append(np.log(f) - np.log(f).mean())  # transformation log-ratio centrée
    X = np.array(vecs)
    X -= X.mean(0)
    u, s, _ = np.linalg.svd(X, full_matrices=False)
    pcs = u[:, :2] * s[:2]
    df = pd.DataFrame(rows, columns=["contig_id", "length", "gc"])
    df[["pc1", "pc2"]] = pcs
    return df


def header_coverage(assembly):
    """Couverture lue dans les en-têtes : MEGAHIT (multi=) ou SPAdes (_cov_)."""
    out = {}
    for line in open(assembly):
        if line.startswith(">"):
            m = re.search(r"multi=([0-9.]+)", line) or re.search(r"_cov_([0-9.]+)", line)
            if m:
                out[line[1:].split()[0]] = float(m.group(1))
    return out


def read_bins(path, pattern="*.f*a*"):
    p = Path(path)
    if p.is_dir():
        rows = [(n, f.stem) for f in sorted(p.glob(pattern)) for n, _ in read_fasta(f)]
        df = pd.DataFrame(rows, columns=["contig_id", "bin_id"])
    else:
        df = pd.read_csv(p, sep="\t", header=None, names=["contig_id", "bin_id"])
    dup = df.contig_id.duplicated().sum()
    if dup:
        print(f"Attention : {dup} contigs présents dans plusieurs bins (premier bin gardé). "
              "Utilisez --bins-glob pour ne garder qu'un seul binner.")
    return df.drop_duplicates("contig_id")


def read_depth(path):
    d = pd.read_csv(path, sep="\t")
    d = d.rename(columns={d.columns[0]: "contig_id", "totalAvgDepth": "depth"})
    return d[["contig_id", "depth"]]


def read_cutoffs(path, factor=1.0):
    """Seuils de score par marqueur (fichier BUSCO scores_cutoff : marqueur<TAB>score)."""
    return {l.split()[0]: factor * float(l.split()[1]) for l in open(path) if len(l.split()) >= 2}


def parse_tblout(tbl, cutoffs=None):
    """Meilleur marqueur par gène ; un hit sous le seuil de son marqueur est ignoré."""
    best = {}
    for line in open(tbl):
        if line.startswith("#"):
            continue
        f = line.split()
        gene, marker, ev, score = f[0], f[2], float(f[4]), float(f[5])
        if cutoffs and score < cutoffs.get(marker, 0):
            continue
        if gene not in best or ev < best[gene][1]:
            best[gene] = (marker, ev)
    return pd.DataFrame([(g.rsplit("_", 1)[0], v[0]) for g, v in best.items()], columns=["contig_id", "marker"])


def hmm_markers(assembly, hmm, threads, evalue, cutoffs=None, factor=1.0):
    names = re.findall(r"^NAME\s+(\S+)", open(hmm).read(), flags=re.M)
    tmp = tempfile.mkdtemp()
    faa, tbl = f"{tmp}/prot.faa", f"{tmp}/hits.tbl"
    subprocess.run(["prodigal", "-i", assembly, "-a", faa, "-p", "meta", "-q", "-o", "/dev/null"], check=True)
    subprocess.run(["hmmsearch", "--cpu", str(threads), "-E", str(evalue), "--tblout", tbl, "-o", "/dev/null",
                    hmm, faa], check=True)
    return parse_tblout(tbl, read_cutoffs(cutoffs, factor) if cutoffs else None), len(names)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--bins", required=True)
    ap.add_argument("--depth", help="facultatif ; sans lui, la couverture est ignorée")
    ap.add_argument("--cov-from-header", action="store_true", help="couverture lue dans les en-têtes MEGAHIT/SPAdes")
    ap.add_argument("--bins-glob", default="*.f*a*", help="motif des fichiers de bins, ex. 'concoct_*.fasta'")
    ap.add_argument("--hmm")
    ap.add_argument("--cutoffs", help="seuils de score par marqueur (BUSCO scores_cutoff)")
    ap.add_argument("--cutoff-factor", type=float, default=1.0, help="multiplie les seuils (0 = pas de filtre, 1 = seuils BUSCO)")
    ap.add_argument("--markers-tsv")
    ap.add_argument("--n-markers", type=int)
    ap.add_argument("--min-len", type=int, default=1500)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--evalue", type=float, default=1e-10)
    ap.add_argument("--out", default="data/magnet_q.db")
    a = ap.parse_args()

    df = composition(a.assembly, a.min_len)
    if a.depth:
        df = df.merge(read_depth(a.depth), on="contig_id", how="left")
        df["cov"] = np.log(df.depth.fillna(0) + 1)
        df = df.drop(columns="depth")
    elif a.cov_from_header:
        hc = header_coverage(a.assembly)
        print(f"Couverture trouvée dans les en-têtes pour {len(hc)} contigs sur {len(df)}")
        df["cov"] = np.log(df.contig_id.map(hc).fillna(0) + 1)
    else:
        print("Pas de fichier de couverture : la vue « Composition » sera utilisée.")
        df["cov"] = 0.0
    df = df.merge(read_bins(a.bins, a.bins_glob), on="contig_id", how="left")
    print(f"{df.bin_id.notna().sum()} contigs sur {len(df)} appartiennent à un bin")
    df["bin_id"] = df.bin_id.fillna("unbinned")

    if a.markers_tsv:
        mk = pd.read_csv(a.markers_tsv, sep="\t", header=None, names=["contig_id", "marker"])
        n_markers = a.n_markers or mk.marker.nunique()
    elif a.hmm:
        mk, n_markers = hmm_markers(a.assembly, a.hmm, a.threads, a.evalue, a.cutoffs, a.cutoff_factor)
    else:
        ap.error("fournir --hmm ou --markers-tsv")
    mk = mk[mk.contig_id.isin(df.contig_id)]

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    if os.path.exists(a.out):
        os.remove(a.out)
    cx = sqlite3.connect(a.out)
    df.to_sql("contigs", cx, index=False)
    mk.to_sql("markers", cx, index=False)
    cx.execute("CREATE TABLE actions(id INTEGER PRIMARY KEY, ts TEXT, kind TEXT, ids TEXT, src TEXT, dst TEXT,"
               " undone INT DEFAULT 0)")
    cx.execute("CREATE TABLE meta(k TEXT, v TEXT)")
    cx.executemany("INSERT INTO meta VALUES (?,?)", [("n_markers", str(n_markers)),
                                                     ("assembly", os.path.abspath(a.assembly))])
    cx.commit()
    print(f"{len(df)} contigs, {df.bin_id.nunique()} bins, {len(mk)} hits marqueurs ({n_markers} marqueurs) -> {a.out}")


if __name__ == "__main__":
    main()
