#!/usr/bin/env python3

from typing import Dict, Any, Sequence, Tuple, Optional
import os
from datetime import datetime

__all__ = ["render_report"]

# ---- deps ----
import numpy as np
import plotly.graph_objs as go

# ---- styles ----
_INLINE_CSS = """
body { margin:0; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,'Noto Sans KR',Arial,sans-serif; background:#fff; color:#111827; }
.container { max-width:1400px; margin:0 auto; padding:24px; }
.header { background:#fff; border-bottom:1px solid #e5e7eb; position:sticky; top:0; z-index:10; }
.h1 { font-size:24px; font-weight:700; padding:16px 0; color:#111827; letter-spacing:.3px; }
.section { margin-top:24px; background:#fff; border:1px solid #e5e7eb; border-radius:16px; overflow:hidden; }
.section .title { padding:14px 16px; font-size:16px; font-weight:700; background:#f9fafb; border-bottom:1px solid #e5e7eb; color:#111827; }
.section .content { padding:16px; }
/* summary dark */
.summary-box { background:#0b0d12; }
.summary-box .title { background:#0b0d12; color:#fff; border-bottom:1px solid #1f2937; }
pre.summary { background:#0b0d12; color:#fff; border:1px solid #1f2937; padding:12px; border-radius:12px; max-height:320px; overflow:auto; white-space:pre-wrap; line-height:1.5; }
.grid { display:grid; grid-template-columns:1fr 1fr; grid-gap:14px; }
.card { background:#fff; border:1px solid #e5e7eb; border-radius:14px; padding:12px; }
.card .subtitle { font-size:14px; font-weight:600; margin-bottom:6px; color:#374151; }
.plot { width:100%; height:520px; }
.footer { color:#6b7280; font-size:12px; margin-top:24px; text-align:right; }
"""

_HTML_TMPL = """<!doctype html>
<html><head><meta charset="utf-8"/><title>{title}</title>
<style>{css}</style>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
</head>
<body>
<div class="container">
  <div class="header"><div class="h1">{title}</div></div>

  <div class="section summary-box">
    <div class="title">Summary</div>
    <div class="content">
      <pre class="summary">{summary_text}</pre>
    </div>
  </div>

  <div class="section">
    <div class="title">UMAPs (2x2)</div>
    <div class="content"><div class="grid">{umap_cards}</div></div>
  </div>

  <div class="section">
    <div class="title">Reference Markers Heatmap</div>
    <div class="content">{heatmap_content}</div>
  </div>

  <div class="footer">Generated: {generated_at}</div>
</div>
</body></html>"""


# ---- helpers: load UMAP from CSV (FIG_DIR/umap.csv) ----
def _load_umap_csv(base_dir: str) -> Dict[str, Sequence[Any]]:
    import os, csv
    path = os.path.join(base_dir, "Annot_results", "Plots", "umap.csv")
    if not os.path.exists(path):
        return {}
    cols = {}
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter=",")
        for row in reader:
            for k, v in row.items():
                cols.setdefault(k, []).append(v)
    # Convert numeric coords
    def _to_float_list(vals):
        out = []
        for v in vals or []:
            try:
                out.append(float(v))
            except Exception:
                out.append(float("nan"))
        return out
    x = _to_float_list(cols.get("umap_x"))
    y = _to_float_list(cols.get("umap_y"))
    # Depth candidates
    depth = cols.get("n_fragment") or cols.get("passed_filters") or None
    # Majority candidates
    majority = cols.get("majority") or cols.get("best_label") or None
    umap = {
        "x": x,
        "y": y,
        "sample": cols.get("sample"),
        "depth": depth,
        "leiden": cols.get("leiden"),
        "majority": majority,
    }
    return umap
# ---- helpers: UMAP normalization ----
def _pick_first(d: Dict[str, Any], keys) -> Optional[Sequence[Any]]:
    for k in keys:
        if k in d and d[k] is not None:
            return d[k]
    return None

def _to_1d_list(arr) -> Sequence[Any]:
    try:
        a = np.asarray(arr)
        if a.ndim == 1:
            return a.tolist()
        if a.ndim == 2 and a.shape[1] == 1:
            return a[:,0].tolist()
        if a.ndim == 2 and a.shape[0] == 1:
            return a[0,:].tolist()
        return a.reshape(-1).tolist()
    except Exception:
        try:
            return list(arr)
        except Exception:
            return []

def _extract_xy(ud: Dict[str, Any]) -> Tuple[Sequence[float], Sequence[float]]:
    x = ud.get("x"); y = ud.get("y")
    if x is not None and y is not None:
        return _to_1d_list(x), _to_1d_list(y)
    # try common embeddings
    for key in ["X","xy","coords","coord","umap","UMAP","obsm_X_umap","X_umap","embedding"]:
        val = ud.get(key)
        if val is None:
            continue
        A = np.asarray(val, dtype=object)
        if A.ndim == 2 and A.shape[1] >= 2:
            try:
                return np.asarray(A[:,0], dtype=float).tolist(), np.asarray(A[:,1], dtype=float).tolist()
            except Exception:
                return [row[0] for row in A], [row[1] for row in A]
        if A.ndim == 1 and len(A) and isinstance(A[0], (list, tuple)) and len(A[0]) >= 2:
            return [row[0] for row in A], [row[1] for row in A]
    return [], []

def _align3(x, y, c):
    n = min(len(x), len(y), len(c))
    if n < 3: return [], [], []
    return x[:n], y[:n], c[:n]

def _normalize_umap_payload(umap: Dict[str, Any]) -> Dict[str, Sequence[Any]]:
    U = dict(umap or {})
    x, y = _extract_xy(U)
    sample   = _pick_first(U, ["sample","Sample","batch","Batch","library","library_id","dataset","sample_id"])
    depth    = _pick_first(U, ["depth","n_fragment","nCount","nCount_ATAC","total_counts","coverage","fragments","reads"])
    leiden   = _pick_first(U, ["leiden","clusters","cluster","louvain","seurat_clusters","cluster_id"])
    majority = _pick_first(U, ["majority","best_label","predicted","prediction","annotation","annot","celltype","CellType","major_label","major"])
    sample   = _to_1d_list(sample) if sample is not None else []
    depth    = _to_1d_list(depth) if depth is not None else []
    leiden   = _to_1d_list(leiden) if leiden is not None else []
    majority = _to_1d_list(majority) if majority is not None else []
    sx, sy, ss = _align3(x, y, sample)   if sample   else ([],[],[])
    dx, dy, dd = _align3(x, y, depth)    if depth    else ([],[],[])
    lx, ly, ll = _align3(x, y, leiden)   if leiden   else ([],[],[])
    mx, my, mm = _align3(x, y, majority) if majority else ([],[],[])
    return {"x": x, "y": y, "sample": ss, "depth": dd, "leiden": ll, "majority": mm}

# ---- plotting ----

def _plotly_umap(x, y, color, name_hint: str) -> str:
    x = np.asarray(x); y = np.asarray(y)
    n = min(x.size, y.size)
    if n < 3:
        return "<div class='plot'>No data</div>"
    x = x[:n]; y = y[:n]
    layout = dict(
        margin=dict(l=10, r=10, t=30, b=10),
        height=520,
        template="plotly_white",
        paper_bgcolor="white",
        plot_bgcolor="white",
        xaxis=dict(showgrid=True, zeroline=False),
        yaxis=dict(showgrid=True, zeroline=False),
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        title=None,
    )
    try:
        c = np.asarray(color, dtype=object) if color is not None else None
        # If this is 'leiden', force categorical rendering (skip numeric casting)
        if (name_hint or '').lower() == 'leiden' and c is not None:
            c = c.astype(str)
            cats = list(dict.fromkeys([str(v) for v in c]))
            traces = []
            cstr = c.astype(str)
            for lab in cats:
                msk = (cstr == lab)
                cnt = int(msk.sum()) if hasattr(msk, 'sum') else sum(bool(v) for v in msk)
                traces.append(go.Scattergl(x=x[msk], y=y[msk], mode="markers",
                                           marker=dict(size=6, opacity=0.85), name=f"{lab} ({cnt})"))
            fig = go.Figure(data=traces)
            fig.update_layout(**layout)
            return fig.to_html(full_html=False, include_plotlyjs=False)
        ## LEIDEN EARLY CATEGORICAL ##



        if c is None or c.size != n:
            fig = go.Figure(data=go.Scattergl(x=x, y=y, mode="markers",
                                              marker=dict(size=6, opacity=0.85),
                                              name=(name_hint or "UMAP")))
            fig.update_layout(**layout)
            return fig.to_html(full_html=False, include_plotlyjs=False)
        # numeric
        try:
            cnum = c.astype(float)
            if (name_hint or '').lower() == 'depth':
                import numpy as _np
                clog = _np.log10(_np.clip(cnum, 1.0, None))
                _cd = _np.column_stack([cnum, clog])
                fig = go.Figure(data=go.Scattergl(
                    x=x, y=y, mode="markers",
                    marker=dict(size=6, color=clog, colorscale="Blues", showscale=True, opacity=0.85,
                                colorbar=dict(title="log10(depth)")),
                    name=name_hint,
                    customdata=_cd,
                    hovertemplate="z=%{customdata[0]:.3g} (log10=%{customdata[1]:.3f})<br>x=%{x:.2f}<br>y=%{y:.2f}<extra></extra>"
                ))
            else:
                fig = go.Figure(data=go.Scattergl(
                    x=x, y=y, mode="markers",
                    marker=dict(size=6, color=cnum, colorscale="Blues", showscale=True, opacity=0.85),
                    name=name_hint
                ))
            fig.update_layout(**layout)
            return fig.to_html(full_html=False, include_plotlyjs=False)
        except Exception:
            # categorical (force leiden to categorical)
            if (name_hint or '').lower() == 'leiden':
                c = c.astype(str)
            cats = list(dict.fromkeys([str(v) for v in c]))
            traces = []
            cstr = c.astype(str)
            for lab in cats:
                msk = (cstr == lab)
                cnt = int(msk.sum()) if hasattr(msk, 'sum') else sum(bool(v) for v in msk)
                traces.append(go.Scattergl(x=x[msk], y=y[msk], mode="markers",
                                           marker=dict(size=6, opacity=0.85), name=f"{lab} ({cnt})"))
            fig = go.Figure(data=traces)
            fig.update_layout(**layout)
            return fig.to_html(full_html=False, include_plotlyjs=False)
    except Exception:
        fig = go.Figure(data=go.Scattergl(x=x, y=y, mode="markers",
                                          marker=dict(size=6, opacity=0.85),
                                          name=(name_hint or "UMAP")))
        fig.update_layout(**layout)
        return fig.to_html(full_html=False, include_plotlyjs=False)


def _plotly_heatmap(z, x_labels=None, y_labels=None, height=480) -> str:
    z = np.asarray(z, dtype=float)
    if z.ndim != 2 or z.size == 0:
        return "<div class='plot'>No heatmap data</div>"
    x = [str(v) for v in (x_labels if x_labels is not None else list(range(z.shape[1])))]
    y = [str(v) for v in (y_labels if y_labels is not None else list(range(z.shape[0])))]
    fig = go.Figure(data=go.Heatmap(z=z, x=x, y=y, colorscale="RdBu", reversescale=True, zmid=0))
    fig.update_layout(margin=dict(l=80, r=10, t=30, b=80), height=height,
                      template="plotly_white", paper_bgcolor="white", plot_bgcolor="white",
                      xaxis=dict(showgrid=False), yaxis=dict(showgrid=False))
    return fig.to_html(full_html=False, include_plotlyjs=False)

# ---- public API ----
def render_report(summary_text: str,
                  umap: Dict[str, Sequence[Any]],
                  heatmap: Dict[str, Any],
                  base_dir: str = "Test",
                  out_name: str = "annot_summary.html",
                  title: str = "Annotation Summary") -> str:
    """Render interactive report and return absolute output path."""
    # Normalize provided UMAP; if missing, load from CSV
    U_in = _normalize_umap_payload(umap) if umap else {}
    if not U_in.get('x') or not U_in.get('y'):
        U_in = _normalize_umap_payload(_load_umap_csv(base_dir))
    umap = U_in
    x = umap.get('x') or [];
    y = umap.get('y') or []

    def _card(lbl, vec):
        if vec:
            html = _plotly_umap(x, y, vec, lbl)
        else:
            # still plot points if x/y exist
            html = _plotly_umap(x, y, None, lbl)
        return f"<div class='card'><div class='subtitle'>{lbl}</div><div class='plot'>{html}</div></div>"

    cards = [
        _card("sample", umap.get("sample")),
        _card("depth", umap.get("depth")),
        _card("leiden", umap.get("leiden")),
        _card("majority", umap.get("majority")),
    ]

    heat_html = "<div class='plot'>No heatmap data</div>"
    if heatmap and ("z" in heatmap):
        heat_html = _plotly_heatmap(heatmap["z"], heatmap.get("x_labels"), heatmap.get("y_labels"), height=480)

    html = _HTML_TMPL.format(
        title=title,
        css=_INLINE_CSS,
        summary_text=(summary_text or ""),
        umap_cards="".join(cards),
        heatmap_content=heat_html,
        generated_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    )

    out_dir = os.path.join(base_dir, "Annot_results")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, out_name)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    return os.path.abspath(out_path)

# Safety: ensure attribute exists even if import aliasing behaves oddly
# (Some loaders might shadow names; this guarantees the symbol.)
render_report = render_report
