"""Bilan avant/après curation d'un bin.   python report.py 5"""
import sys

import pandas as pd

from core import Curator, label


def main():
    b = sys.argv[1] if len(sys.argv) > 1 else "5"
    db = sys.argv[2] if len(sys.argv) > 2 else "data/magnet_q.db"
    cur = Curator(db)
    orig = pd.read_sql("SELECT contig_id, bin_id FROM contigs", cur.cx).set_index("contig_id").bin_id
    ids0 = orig[orig == b].index  # composition d'origine, avant toute action

    def q(ids):
        counts = cur.markers[cur.markers.contig_id.isin(ids)].marker.value_counts()
        comp = 100 * len(counts) / cur.n_markers
        cont = 100 * (counts.sum() - len(counts)) / cur.n_markers
        return comp, cont, label(comp, cont), len(ids), cur.contigs.loc[ids, "length"].sum() / 1e6

    a, z = q(ids0), q(cur.members(b).index)
    n_act = cur.cx.execute("SELECT COUNT(*) FROM actions WHERE undone=0").fetchone()[0]
    print(f"Bin {b} : {n_act} actions enregistrées")
    print(f"{'':8}{'compl.':>8}{'cont.':>8}{'classe':>8}{'contigs':>9}{'Mb':>7}")
    for name, r in (("avant", a), ("après", z)):
        print(f"{name:8}{r[0]:8.1f}{r[1]:8.1f}{r[2]:>8}{r[3]:9d}{r[4]:7.2f}")
    dc, dp = a[1] - z[1], a[0] - z[0]
    print(f"\nContamination : -{dc:.1f} points ; complétude : -{dp:.1f} points ; "
          f"contigs retirés : {a[3] - z[3]} ({100 * (a[3] - z[3]) / a[3]:.0f} %)")
    if dc > 0:
        print(f"Coût : {dp / dc:.2f} point de complétude perdu par point de contamination retiré")


if __name__ == "__main__":
    main()
