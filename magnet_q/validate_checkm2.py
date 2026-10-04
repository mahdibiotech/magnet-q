"""Compare l'estimation rapide de MAGNET-Q à CheckM2.

  python validate_checkm2.py            # exporte les bins en FASTA, lance checkm2, compare
  python validate_checkm2.py --report existing/quality_report.tsv   # réutilise un rapport CheckM2 déjà produit
À lancer avant et après curation (changer --tag) pour mesurer le gain réel.
"""
import argparse
import os
import re
import subprocess

import pandas as pd

from core import UNBINNED, Curator
from prepare_data import read_fasta


def export_bins(cur, outdir):
    os.makedirs(outdir, exist_ok=True)
    wanted = cur.contigs.bin_id.to_dict()
    handles = {}
    for name, seq in read_fasta(cur.assembly):
        b = wanted.get(name)
        if b and b != UNBINNED:
            if b not in handles:
                handles[b] = open(f"{outdir}/{b}.fa", "w")
            handles[b].write(f">{name}\n{seq}\n")
    for h in handles.values():
        h.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default="data/magnet_q.db")
    ap.add_argument("--tag", default="curated")
    ap.add_argument("--report")
    ap.add_argument("--checkm-stats", help="tableau Bin statistics de CheckM (lineage_wf)")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    cur = Curator(a.db)
    report = a.report
    if a.checkm_stats:
        key = lambda s: re.findall(r"\d+", str(s))[-1]
        ref = pd.read_csv(a.checkm_stats, sep="\t")
        c2 = ref.assign(k=ref["Bin Id"].map(key)).set_index("k")[["Completeness", "Contamination"]]
        est = pd.DataFrame({b: cur.quality(b) for b in cur.bins()}).T[["completeness", "contamination"]].astype(float)
        est.index = [key(b) for b in est.index]
        report = "stats"
    elif not report:
        if not cur.assembly or not os.path.exists(cur.assembly):
            raise SystemExit("Assemblage introuvable : relancer prepare_data.py avec --assembly.")
        d = f"data/bins_{a.tag}"
        export_bins(cur, d)
        subprocess.run(["checkm2", "predict", "-i", d, "-x", "fa", "-o", f"data/checkm2_{a.tag}",
                        "--threads", str(a.threads), "--force"], check=True)
        report = f"data/checkm2_{a.tag}/quality_report.tsv"
    if report != "stats":
        c2 = pd.read_csv(report, sep="\t", dtype={"Name": str}).set_index("Name")[["Completeness", "Contamination"]]
        est = pd.DataFrame({b: cur.quality(b) for b in cur.bins()}).T[["completeness", "contamination"]].astype(float)
    df = est.join(c2, how="inner")
    df.columns = ["est_compl", "est_cont", "checkm2_compl", "checkm2_cont"]
    os.makedirs("data", exist_ok=True)
    df.round(1).to_csv(f"data/validation_{a.tag}.csv")
    print(df.round(1).to_string())
    for k in ("compl", "cont"):
        e, c = df[f"est_{k}"], df[f"checkm2_{k}"]
        print(f"{k}: r = {e.corr(c):.2f}, écart absolu moyen = {(e - c).abs().mean():.1f} points (n = {len(df)})")


if __name__ == "__main__":
    main()
