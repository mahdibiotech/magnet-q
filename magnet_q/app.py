"""MAGNET-Q : panneau qualité en temps réel pour la curation de MAGs.  Lancer : python app.py"""
import numpy as np
import plotly.graph_objects as go
from dash import Dash, Input, Output, State, ctx, dcc, html, no_update

from core import UNBINNED, Curator

cur = Curator()
START = max(cur.bins(), key=lambda b: cur.quality(b)["contamination"])  # bin le plus mélangé
START_VIEW = "cov" if cur.contigs["cov"].nunique() > 1 else "tnf"
app = Dash(__name__, title="MAGNET-Q")
server = app.server
COLORS = {"HQ": "#1B7F5C", "MQ": "#B07A00", "LQ": "#B3261E"}


def bar(value, maxv, marks, color):
    """Barre de jauge avec repères de seuils MIMAG."""
    pct = min(100, 100 * value / maxv)
    ticks = [html.Span(className="tick", style={"left": f"{100 * m / maxv}%"}) for m in marks]
    return html.Div([html.Div(className="fill", style={"width": f"{pct}%", "background": color}), *ticks],
                    className="gauge")


def quality_panel(q):
    c = COLORS[q["label"]]
    return [
        html.Div([html.Span("Complétude"), html.B(f"{q['completeness']:.1f} %")], className="row"),
        bar(q["completeness"], 100, [50, 90], c),
        html.Div([html.Span("Contamination"), html.B(f"{q['contamination']:.1f} %")], className="row"),
        bar(q["contamination"], 20, [5, 10], c),
        html.Div(f"Qualité {q['label']}", className="verdict", style={"color": c}),
        html.Div(f"{q['n']} contigs · {q['size'] / 1e6:.2f} Mb", className="muted"),
    ]


app.layout = html.Div([
    dcc.Store(id="version", data=0),
    html.Header([
        html.H1("MAGNET-Q"),
        dcc.Dropdown(id="bin", options=cur.bins(), value=START, clearable=False, className="dd"),
        html.Span(className="spacer"),
        html.Button("Annuler", id="undo"), html.Button("Rétablir", id="redo"),
        html.Button("Exporter", id="export"), html.Span(id="exported", className="muted"),
    ]),
    html.Details([
        html.Summary("Comment utiliser cette page (30 secondes)"),
        html.Ol([
            html.Li("Chaque bin est un génome de bactérie reconstruit à partir de morceaux d'ADN (les points). "
                    "Parfois, des morceaux d'une autre bactérie s'y sont glissés : c'est la contamination."),
            html.Li("Sur le graphique, les morceaux d'un même génome forment un groupe. "
                    "Les points isolés ou en bordure sont les intrus probables. Le tableau « Contigs suspects » "
                    "à droite les liste du plus au moins probable."),
            html.Li("Entourez des points à la souris (outil lasso, en haut du graphique). "
                    "La ligne sous le graphique annonce le résultat avant que vous ne validiez."),
            html.Li("Si la contamination baisse sans perdre trop de complétude, cliquez sur « Retirer du bin ». "
                    "Vous pouvez toujours annuler."),
        ]),
        html.Div([html.B("Complétude : "), "part des gènes attendus (une copie chacun) retrouvés dans le bin. ",
                  html.B("Contamination : "), "part de gènes présents en double, signe d'un mélange. ",
                  html.B("HQ : "), "complétude > 90 % et contamination < 5 %."], className="muted"),
    ], open=True, className="guide"),
    html.Main([
        html.Section([
            dcc.RadioItems(id="view", value=START_VIEW, inline=True, className="muted", options=[
                {"label": " GC × couverture", "value": "cov"},
                {"label": " Composition (tétranucléotides, PCA)", "value": "tnf"}]),
            dcc.Graph(id="scatter", config={"displaylogo": False}, style={"height": "58vh"}),
            html.Div(id="whatif", className="whatif"),
            html.Div([
                html.Button("Retirer du bin", id="remove", className="primary"),
                dcc.Dropdown(id="target", placeholder="Déplacer vers…", className="dd"),
                html.Button("Déplacer", id="move"),
                html.Button("Créer un nouveau bin", id="new"),
            ], className="actions"),
        ], className="plot"),
        html.Aside([
            html.H2("Qualité du bin"), html.Div(id="quality"),
            html.H2("Contigs suspects"),
            html.Div("Retirer ces contigs fait baisser la contamination (colonne 2) au prix d'une perte de complétude (colonne 3).", className="muted"),
            html.Div(id="suspects"),
            html.H2("Journal"), html.Div(id="journal", className="journal"),
        ]),
    ]),
])


def selected_ids(sel, bin_id):
    ids = [p.get("customdata") for p in (sel or {}).get("points", [])]
    members = set(cur.members(bin_id).index)
    return [i for i in ids if i in members]


@app.callback(Output("version", "data"),
              Input("remove", "n_clicks"), Input("move", "n_clicks"), Input("new", "n_clicks"),
              Input("undo", "n_clicks"), Input("redo", "n_clicks"),
              State("bin", "value"), State("scatter", "selectedData"), State("target", "value"),
              State("version", "data"), prevent_initial_call=True)
def act(_r, _m, _n, _u, _d, bin_id, sel, target, v):
    t = ctx.triggered_id
    ids = selected_ids(sel, bin_id)
    if t == "undo":
        cur.undo()
    elif t == "redo":
        cur.redo()
    elif t == "remove":
        cur.move(ids, bin_id, UNBINNED, "retrait")
    elif t == "move" and target:
        cur.move(ids, bin_id, target, "déplacement")
    elif t == "new":
        cur.move(ids, bin_id, cur.new_bin_name(), "nouveau bin")
    return v + 1


@app.callback(Output("scatter", "figure"), Output("quality", "children"), Output("suspects", "children"),
              Output("journal", "children"), Output("target", "options"), Output("bin", "options"),
              Input("bin", "value"), Input("version", "data"), Input("view", "value"))
def render(bin_id, _v, view):
    try:
        return _render(bin_id, view)
    except Exception:  # affiche l'erreur dans la page au lieu d'un écran vide
        import traceback
        return go.Figure(), html.Pre(traceback.format_exc(), className="muted"), "", "", [], no_update


def _render(bin_id, view):
    m = cur.members(bin_id)
    x, y, xt, yt = ((m.gc * 100, m["cov"], "GC (%)", "ln(couverture)") if view == "cov"
                    else (m.pc1, m.pc2, "PC1 composition", "PC2 composition"))
    fig = go.Figure(go.Scattergl(
        x=np.asarray(x, dtype=float), y=np.asarray(y, dtype=float), mode="markers", customdata=[str(i) for i in m.index],
        marker=dict(size=(m.length / m.length.max() * 14 + 4).to_numpy(dtype=float), color="#0B5C66", opacity=0.65),
        selected=dict(marker=dict(color="#C2185B", opacity=1)),
        unselected=dict(marker=dict(opacity=0.35)),
        hovertemplate="%{customdata}<br>%{x:.2f} · %{y:.2f}<extra></extra>"))
    fig.update_layout(dragmode="lasso", margin=dict(l=50, r=10, t=10, b=45), paper_bgcolor="rgba(0,0,0,0)",
                      plot_bgcolor="#F7F9FA", xaxis_title=xt, yaxis_title=yt,
                      uirevision=f"{bin_id}{view}", font=dict(family="Instrument Sans, system-ui, sans-serif"))
    s = cur.suspects(bin_id)
    sus = html.Table([html.Tr([html.Th("Contig"), html.Th("Contamination (points)"), html.Th("Complétude (points)")])] +
                     [html.Tr([html.Td(i), html.Td(f"{r.d_cont:.1f}"), html.Td(f"{r.d_comp:.1f}")])
                      for i, r in s.iterrows()]) if len(s) else html.Div("Aucun marqueur dupliqué.", className="muted")
    j = [html.Div(x) for x in cur.journal()] or [html.Div("Aucune action pour l'instant.", className="muted")]
    others = [b for b in cur.bins() if b != bin_id]
    return fig, quality_panel(cur.quality(bin_id)), sus, j, others, cur.bins()


@app.callback(Output("whatif", "children"), Input("scatter", "selectedData"),
              Input("bin", "value"), Input("version", "data"))
def whatif(sel, bin_id, _v):
    ids = selected_ids(sel, bin_id)
    if not ids:
        return html.Div("Entourez des points sur le graphique : l'effet du retrait s'affichera ici.", className="muted")
    a, b = cur.quality(bin_id), cur.quality(bin_id, exclude=ids)
    kb = cur.contigs.loc[ids, "length"].sum() / 1e3
    return [html.Div(f"{len(ids)} contigs · {kb:.0f} kb sélectionnés", className="muted"),
            html.Div([html.Span("Si retirés  "),
                      html.B(f"contamination {a['contamination']:.1f} → {b['contamination']:.1f} %"),
                      html.Span("  ·  "),
                      html.B(f"complétude {a['completeness']:.1f} → {b['completeness']:.1f} %"),
                      html.Span(f"  ·  {a['label']} → {b['label']}")], className="delta")]


@app.callback(Output("exported", "children"), Input("export", "n_clicks"), prevent_initial_call=True)
def export(_):
    return f"Écrit : {cur.export()}"


if __name__ == "__main__":
    app.run(debug=False, port=8050)
