#!/usr/bin/env python3

"""
Final single-file report generator (CSV-first), v4.0
  - Prevent gene labels from overlapping titles by lifting outside-title higher and increasing top margin.
  - Bottom type strip colors switched to a softer qualitative palette (Dark24 → Set3 → Glasbey fallback) with slight lightening/desaturation.
  - Type names rendered as centered annotations over contiguous color runs in the strip.
  - Leiden table scroller tightened; UMAPs remain square.

CLI:
  python run_finalReport.py -b <BASE_DIR> -o /custom/out -r /path/markers.csv
"""

import argparse
import re
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.colors as pc

# ---- Utils ----
def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def fmt_int(x) -> str:
    try:
        return f"{int(x):,d}"
    except Exception:
        return str(x)

def now_kst() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def read_text(p: Path) -> Optional[str]:
    for enc in ("utf-8", "utf-8-sig"):
        try:
            return p.read_text(encoding=enc)
        except Exception:
            pass
    return None

def resolve_paths(base_dir: Path) -> Dict[str, Path]:
    return {
        "base": base_dir,
        "filter_summary": base_dir / "Filter_results" / "summary.txt",
        "umap_csv": base_dir / "Annot_results" / "plots" / "umap.csv",
        "dot_dir": base_dir / "Annot_results" / "dotplots",
        "dot_frac": base_dir / "Annot_results" / "dotplots" / "dot_frac_leiden_dendro.csv",
        "dot_mean": base_dir / "Annot_results" / "dotplots" / "dot_mean_leiden_dendro.csv",
        "label_map": base_dir / "Annot_results" / "dotplots" / "label_map.csv",
        "update_summary": base_dir / "Annot_results" / "update_h5ad" / "update_summary.csv",
        "default_report_dir": base_dir / "Annot_results" / "report",
    }

# ---- HTML and CSS ----
BASE_CSS = """
:root { --bg:#ffffff; --fg:#101316; --muted:#60656f; --card:#f7f7f8; --border:#e5e7eb; }
*{box-sizing:border-box}
body{background:var(--bg);color:var(--fg);font:14px/1.55 system-ui,-apple-system,"Segoe UI",Roboto,"Noto Sans",Arial,"Apple SD Gothic Neo","Malgun Gothic","맑은 고딕",sans-serif;margin:0}
header{position:sticky;top:0;z-index:10;background:rgba(255,255,255,0.95);backdrop-filter:blur(6px);border-bottom:1px solid var(--border)}
.header-inner{display:flex;align-items:center;justify-content:space-between;gap:12px;padding:8px 14px}
.title{font-weight:800;font-size:18px}
.meta{color:var(--muted);font-size:12px}
.container{display:flex}
nav{width:184px;border-right:1px solid var(--border);padding:12px 10px;position:sticky;top:48px;height:calc(100vh - 48px);overflow:auto}
nav a{display:block;padding:6px 8px;color:var(--fg);text-decoration:none;border-radius:8px;font-size:13px}
nav a:hover{background:var(--card)}
main{flex:1;padding:16px 18px}
section{margin:18px 0}
h2{margin:8px 0 10px;font-size:20px;font-weight:800}
h3{margin:6px 0 8px;font-size:16px;color:#222;font-weight:700}
.subtle{color:var(--muted);font-size:12.5px;margin-top:-4px}

/* Cards */
.cards{display:grid;grid-template-columns:repeat(3,1fr);gap:8px;margin:0 0 8px}
.card{background:var(--card);border:1px solid var(--border);border-radius:8px;padding:8px 10px}
.card .k{color:var(--muted);font-size:12px}
.card .v{font-weight:800;font-size:18px}

/* Tables */
.table-wrap{margin:6px 0 10px;overflow:auto}
.table-wrap.narrow{max-width:920px}
.table-wrap.wide{max-width:none}
.table-wrap.tall{max-height:360px;overflow-y:auto}
.table-caption{color:var(--muted);font-size:12px;margin:2px 0 6px}
table.tbl{border-collapse:collapse;width:100%;min-width:320px;table-layout:fixed}
table.tbl th,table.tbl td{border:1px solid var(--border);padding:4px 6px;text-align:left;font-size:13px;word-break:break-word}
table.tbl th{background:#fafafa;font-weight:700}
table.tbl td.num, table.tbl th.num { text-align:center; }

/* Layouts */
.grid-1x3{display:grid;grid-template-columns:repeat(3,1fr);gap:10px;align-items:stretch}
.grid-1x2{display:grid;grid-template-columns:repeat(2,1fr);gap:12px;align-items:stretch}
.figure{background:#fff;border:1px solid var(--border);border-radius:10px;padding:6px}
.vcol{display:flex;flex-direction:column;gap:8px}
.block{background:#fff;border:1px solid var(--border);border-radius:10px;padding:10px}

code.small{background:#f3f4f6;padding:1px 6px;border-radius:6px;border:1px solid var(--border);font-size:12px}
.muted{color:var(--muted);font-style:italic;margin:6px 0}
.hr{height:1px;background:var(--border);margin:10px 0}
footer{color:var(--muted);font-size:12px;padding:12px;border-top:1px solid var(--border);margin-top:12px}
"""

def build_nav(items: List[Tuple[str, str]]) -> str:
    return "<nav>" + "\n".join([f'<a href="#{i}">{t}</a>' for i, t in items]) + "</nav>"

def wrap_section(section_id: str, title: str, subtitle: Optional[str], inner_html: str) -> str:
    sub = f'<div class="subtle">{subtitle}</div>' if subtitle else ""
    return f"""
<section id="{section_id}">
  <h2>{title}</h2>
  {sub}
  <div>{inner_html}</div>
</section>
"""

def assemble_html(project: str, base_dir: Path, sections: List[str]) -> str:
    nav = build_nav([("filter", "Filter & Clustering"), ("annotation", "Annotation"), ("update", "update_h5ad")])
    head = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>Final Report - {project}</title>
<style>{BASE_CSS}</style>
</head>
<body>
<header>
  <div class="header-inner">
    <div class="title">Final Report — {project}</div>
    <div class="meta">Generated: {now_kst()} · BASE_DIR: <code class="small">{base_dir}</code></div>
  </div>
</header>
<div class="container">
{nav}
<main>
"""
    tail = """
</main>
</div>
<footer>Fully-embedded FinalReport. © 2025 Computational Genomics & AI Lab</footer>
</body>
</html>
"""
    return head + "\n".join(sections) + tail

# ---- Filter section ----
def parse_filter_summary(txt: str) -> Tuple[Optional[float], Optional[Tuple[int, int]], Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    if not txt:
        return None, None, None, None

    tss_cut = None
    umi_pair = None

    m = re.search(r"TSS\s*cutoff\s*:\s*([0-9.]+)", txt)
    if m:
        try:
            tss_cut = float(m.group(1))
        except Exception:
            tss_cut = None

    m = re.search(r"UMI\s*min/max\s*:\s*([0-9,]+)\s*/\s*([0-9,]+)", txt)
    if m:
        try:
            umi_pair = (int(m.group(1).replace(",", "")), int(m.group(2).replace(",", "")))
        except Exception:
            umi_pair = None

    def parse_block(title_prefix: str) -> Optional[pd.DataFrame]:
        pat = rf"{re.escape(title_prefix)}.*?(?:(?:\n\n)|\Z)"
        mm = re.search(pat, txt, flags=re.DOTALL | re.IGNORECASE)
        if not mm:
            return None
        block = mm.group(0)
        rows = []
        for line in block.splitlines():
            line = line.strip()
            if line.startswith("- "):
                m2 = re.match(r"-\s*([^:]+):\s*([0-9,]+)\s*→\s*([0-9,]+)\s*\(retention=([0-9.]+)\)", line)
                if m2:
                    name = m2.group(1).strip()
                    pre = int(m2.group(2).replace(",", ""))
                    post = int(m2.group(3).replace(",", ""))
                    ret = float(m2.group(4))
                    rows.append([name, pre, post, ret])
            elif line.upper().startswith("TOTAL"):
                m3 = re.match(r"TOTAL:\s*([0-9,]+)\s*→\s*([0-9,]+)\s*\(retention=([0-9.]+)\)", line, flags=re.IGNORECASE)
                if m3:
                    pre = int(m3.group(1).replace(",", ""))
                    post = int(m3.group(2).replace(",", ""))
                    ret = float(m3.group(3))
                    rows.append(["TOTAL", pre, post, ret])
        if rows:
            return pd.DataFrame(rows, columns=["name", "pre", "post", "retention"])
        return None

    df_A = parse_block("A) Filtering")
    df_B = parse_block("B) Doublet removal")
    return tss_cut, umi_pair, df_A, df_B

def _colgroup_html(cols: List[str]) -> str:
    key = [c.lower() for c in cols]
    if key == ["name", "pre", "post", "retention"]:
        return '<colgroup><col style="width:38%"><col style="width:20%"><col style="width:20%"><col style="width:22%"></colgroup>'
    w = max(10.0, 100.0 / max(1, len(cols)))
    return '<colgroup>' + ''.join([f'<col style="width:{w:.1f}%">' for _ in cols]) + '</colgroup>'

def df_to_html_table(df: pd.DataFrame, caption: str, extra_class: str = "narrow") -> str:
    if df is None or df.empty:
        return f'<div class="muted">No data for {caption}</div>'
    fmt_df = df.copy()
    num_cols = set()
    for c in fmt_df.columns:
        if c in ("pre", "post"):
            fmt_df[c] = fmt_df[c].map(fmt_int)
            num_cols.add(c)
        elif c == "retention":
            fmt_df[c] = fmt_df[c].map(lambda x: f"{float(x):.4f}")
            num_cols.add(c)
        elif pd.api.types.is_numeric_dtype(df[c]):
            num_cols.add(c)

    head = "".join([f"<th class='num'>{c}</th>" if c in num_cols else f"<th>{c}</th>" for c in fmt_df.columns])
    rows = []
    for _, r in fmt_df.iterrows():
        tds = "".join([
            (f"<td class='num'>{r[c]}</td>" if c in num_cols else f"<td>{r[c]}</td>")
            for c in fmt_df.columns
        ])
        rows.append(f"<tr>{tds}</tr>")
    body = "\n".join(rows)
    colgroup = _colgroup_html(list(fmt_df.columns))
    return f"""
<div class="table-wrap {extra_class}">
  <div class="table-caption">{caption}</div>
  <table class="tbl">
    {colgroup}
    <thead><tr>{head}</tr></thead>
    <tbody>
      {body}
    </tbody>
  </table>
</div>
"""

# ---- Clustering (UMAP) ----

def load_umap_csv(p: Path) -> Tuple[pd.DataFrame, List[str]]:
    errs: List[str] = []
    if not p.exists():
        errs.append(f"Missing: {p}")
        return pd.DataFrame(), errs
    try:
        df = pd.read_csv(p)
    except Exception as e:
        errs.append(f"Failed to read {p.name}: {e}")
        return pd.DataFrame(), errs

    required = ["cell", "sample", "n_fragment", "leiden", "umap_x", "umap_y"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        errs.append(f"Required columns missing in umap.csv: {missing}")
        return pd.DataFrame(), errs

    df = df.copy()
    df["cell"] = df["cell"].astype(str)
    df["sample"] = df["sample"].astype(str).replace({"": "NA"})
    df["leiden"] = df["leiden"].astype(str).replace({"": "NA"})
    df["umap_x"] = pd.to_numeric(df["umap_x"], errors="coerce")
    df["umap_y"] = pd.to_numeric(df["umap_y"], errors="coerce")
    df["n_fragment"] = pd.to_numeric(df["n_fragment"], errors="coerce")
    df = df.dropna(subset=["umap_x", "umap_y"])
    return df, errs

def group_counts(series: pd.Series) -> Dict[str, int]:
    return series.value_counts(dropna=False).astype(int).to_dict()

def limit_points(df: pd.DataFrame, max_points: int = 50000, by: Optional[str] = None) -> pd.DataFrame:
    n = len(df)
    if n <= max_points:
        return df
    if by is None or by not in df.columns:
        return df.sample(n=max_points, random_state=1337).sort_index()
    vc = df[by].astype(str).value_counts()
    total = vc.sum()
    target = (vc / total * max_points).round().astype(int)
    parts = []
    for g, k in target.items():
        sub = df[df[by].astype(str) == g]
        parts.append(sub if len(sub) <= k else sub.sample(n=int(k), random_state=1337))
    return pd.concat(parts, axis=0)

class PlotlyOnce:
    def __init__(self):
        self.used = False
    def to_html(self, fig, full_html=False):
        html = pio.to_html(fig, full_html=full_html, include_plotlyjs=(not self.used))
        self.used = True
        return html

def title_outside(fig: go.Figure, text: str, top_margin: int = 180) -> None:
    fig.update_layout(title=None, margin=dict(t=top_margin))
    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.0, y=1.26, xanchor="left", yanchor="top",
        showarrow=False, text=text,
        font=dict(size=18, family="Arial, sans-serif")
    )

def _apply_grid_layout(fig: go.Figure) -> None:
    fig.update_xaxes(showgrid=True, zeroline=False)
    fig.update_yaxes(showgrid=True, zeroline=False, scaleanchor="x", scaleratio=1)

def glasbey_colors(n: int) -> List[str]:
    base = [
        "#000000","#ff0000","#00ff00","#0000ff","#ffff00","#ff00ff","#00ffff","#ffffff",
        "#7f0000","#007f00","#00007f","#7f7f00","#7f007f","#007f7f","#7f7f7f","#ff7f7f",
        "#7fff7f","#7f7fff","#ffff7f","#ff7fff","#7fffff","#bfbfbf","#400000","#004000",
        "#000040","#404000","#400040","#004040","#404040","#bf0000","#00bf00","#0000bf",
        "#bfbf00","#bf00bf","#00bfbf","#bfbfbf","#ffbf7f","#bfff7f","#bf7fff","#ffffbf",
        "#ffbfff","#bfffff","#3f3f3f","#9f0000","#009f00","#00009f","#9f9f00","#9f009f",
        "#009f9f","#9f9f9f","#ff9f7f","#9fff7f","#9f7fff","#ffff9f","#ff9fff","#9fffff",
        "#1f1f1f","#df0000","#00df00","#0000df","#dfdf00","#df00df","#00dfdf","#dfdfdf"
    ]
    if n <= len(base):
        return base[:n]
    out, i = [], 0
    while len(out) < n:
        out.append(base[i % len(base)])
        i += 1
    return out

def make_umap_by_depth(df: pd.DataFrame, title: str, p1: PlotlyOnce, fig_size: int = 560) -> str:
    if df is None or df.empty:
        return f'<div class="muted">No UMAP data for {title}</div>'
    show = limit_points(df, max_points=50000)
    logdepth = np.log10(show["n_fragment"].fillna(0).astype(float) + 1.0)
    fig = go.Figure(go.Scattergl(
        x=show["umap_x"], y=show["umap_y"],
        mode="markers",
        marker=dict(
            size=3, opacity=0.85,
            color=logdepth, colorscale="Magma",
            showscale=True, colorbar=dict(title="log10(n_fragment+1)")
        ),
        text=show["cell"],
        hovertemplate="<b>%{text}</b><br>UMAP1=%{x:.3f}<br>UMAP2=%{y:.3f}<br>log10(n_frag+1)=%{marker.color:.3f}<extra></extra>",
        showlegend=False
    ))
    fig.update_layout(
        xaxis=dict(title="UMAP1"),
        yaxis=dict(title="UMAP2"),
        margin=dict(l=10, r=10, t=180, b=10),
        template="plotly_white",
        hovermode="closest", dragmode="pan",
        showlegend=False, font=dict(size=13),
        height=fig_size
    )
    title_outside(fig, title, top_margin=180)
    _apply_grid_layout(fig)
    return p1.to_html(fig, full_html=False)

def make_umap_grouped(df: pd.DataFrame, color_col: str, title: str, p1: PlotlyOnce,
                      palette_map: Optional[Dict[str, str]] = None,
                      fig_size: int = 560) -> str:
    if df is None or df.empty:
        return f'<div class="muted">No UMAP data for {title}</div>'
    show = limit_points(df, max_points=50000, by=color_col if color_col in ("sample", "leiden", "type") else None)
    fig = go.Figure()
    counts = group_counts(show[color_col].astype(str)) if color_col in show.columns else {}
    levels = list(dict.fromkeys(show[color_col].astype(str).tolist()))
    if palette_map:
        palette = {lv: palette_map.get(lv) for lv in levels}
    elif color_col == "leiden":
        palette = dict(zip(levels, glasbey_colors(len(levels))))
    else:
        palette = {lv: None for lv in levels}

    for val, sub in show.groupby(color_col, sort=False):
        name = str(val)
        cnt = counts.get(name, len(sub))
        color = palette.get(name, None)
        fig.add_trace(go.Scattergl(
            x=sub["umap_x"], y=sub["umap_y"],
            mode="markers",
            marker=dict(size=3, opacity=0.78, color=color, line=dict(width=0.2, color="#ffffff")),
            name=f"{name} ({fmt_int(cnt)})",
            text=sub["cell"],
            hovertemplate="<b>%{text}</b><br>UMAP1=%{x:.3f}<br>UMAP2=%{y:.3f}<br>"+f"{color_col}={name}"+ "<extra></extra>",
        ))
    fig.update_layout(
        xaxis=dict(title="UMAP1"),
        yaxis=dict(title="UMAP2"),
        legend=dict(itemsizing="constant"),
        margin=dict(l=10, r=10, t=180, b=10),
        template="plotly_white",
        hovermode="closest", dragmode="pan",
        font=dict(size=13),
        height=fig_size
    )
    title_outside(fig, title, top_margin=180)
    _apply_grid_layout(fig)
    return p1.to_html(fig, full_html=False)


# ---- Marker/type helpers ----

def _hex_to_rgb(h: str) -> Tuple[int, int, int]:
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

def _rgb_to_hex(rgb: Tuple[int, int, int]) -> str:
    r, g, b = [max(0, min(255, int(round(v)))) for v in rgb]
    return f"#{r:02x}{g:02x}{b:02x}"

def _blend_hex(c1: str, c2: str, alpha: float) -> str:
    r1, g1, b1 = _hex_to_rgb(c1); r2, g2, b2 = _hex_to_rgb(c2)
    r = (1-alpha)*r1 + alpha*r2
    g = (1-alpha)*g1 + alpha*g2
    b = (1-alpha)*b1 + alpha*b2
    return _rgb_to_hex((r, g, b))

def soft_palette(colors: List[str]) -> List[str]:
    # Slightly desaturate (towards gray) and lighten (towards white)
    out = []
    for c in colors:
        c1 = _blend_hex(c, "#808080", 0.15)  # desaturate ~15%
        c2 = _blend_hex(c1, "#f5f5f5", 0.05) # lighten ~5%
        out.append(c2)
    return out

def choose_type_palette(n: int) -> List[str]:
    if n <= 24:
        base = pc.qualitative.Dark24[:n]
    elif n <= 32:
        need = n
        base = (pc.qualitative.Set3 + pc.qualitative.Dark24)[:need]
    else:
        base = glasbey_colors(n)
    return soft_palette([pc.convert_colors_to_same_type([c], colortype='rgb')[0] if c.startswith("rgb") else c for c in base])

def parse_marker_types(marker_csv: Optional[Path], gene_cols: List[str]) -> Optional[Dict[str, str]]:
    if marker_csv is None or (not marker_csv.exists()):
        return None
    try:
        mk = pd.read_csv(marker_csv)
    except Exception:
        return None
    mk.columns = [c.lower() for c in mk.columns]
    if "type" not in mk.columns:
        return None

    candidates = [c for c in ("gene_name", "gene", "gene_id") if c in mk.columns]
    if not candidates:
        return None

    gene_cols_set = set(str(g) for g in gene_cols)
    best_col, best_hit = None, -1
    for c in candidates:
        hits = set(str(x) for x in mk[c].astype(str)) & gene_cols_set
        if len(hits) > best_hit:
            best_hit, best_col = len(hits), c

    if best_col is None or best_hit == 0:
        return None

    mk = mk[[best_col, "type"]].dropna()
    mk[best_col] = mk[best_col].astype(str)
    mk["type"] = mk["type"].astype(str)
    return dict(zip(mk[best_col], mk["type"]))

def make_type_palette(type_levels: List[str]) -> Dict[str, str]:
    colors = choose_type_palette(len(type_levels)) if type_levels else []
    pal = dict(zip(type_levels, colors))
    pal["NA"] = "#9aa0a6"
    return pal


# ---- Dotplot layout ----

def _build_discrete_colorscale(colors: List[str]) -> List[Tuple[float, str]]:
    k = max(1, len(colors))
    if k == 1:
        return [(0.0, colors[0]), (1.0, colors[0])]
    stops = []
    for i, c in enumerate(colors):
        a = i / k
        b = (i + 1) / k
        stops.append((a, c)); stops.append((b, c))
    return stops

def _run_lengths(labels: List[str]) -> List[Tuple[int, int, str]]:
    runs = []
    if not labels:
        return runs
    s, cur = 0, labels[0]
    for i in range(1, len(labels)):
        if labels[i] != cur:
            runs.append((s, i-1, cur))
            s, cur = i, labels[i]
    runs.append((s, len(labels)-1, cur))
    return runs

def dotplot_with_typestrip_bottom(df_frac: pd.DataFrame, df_mean: pd.DataFrame, title: str,
                                  p1: PlotlyOnce,
                                  gene_type_map: Optional[Dict[str, str]] = None,
                                  type_palette: Optional[Dict[str, str]] = None) -> str:
    groups = list(df_frac.index.astype(str))
    genes  = list(df_frac.columns.astype(str))

    if gene_type_map is None:
        gene_type_map = {}
    gene_types = [gene_type_map.get(g, "NA") for g in genes]
    type_levels = list(dict.fromkeys(gene_types))
    if type_palette is None:
        type_palette = make_type_palette(type_levels)

    observed_types = list(dict.fromkeys(type_levels + (["NA"] if "NA" in type_palette else [])))
    palette_colors = [type_palette.get(t, "#9aa0a6") for t in observed_types]
    color_index = {t: i for i, t in enumerate(observed_types)}
    strip_vals = [color_index[t] for t in gene_types]

    xs, ys, colors, sizes, hover_txt = [], [], [], [], []
    fvals = df_frac.values.flatten()
    fmin, fmax = float(np.nanmin(fvals)), float(np.nanmax(fvals))
    if np.isnan(fmin) or np.isnan(fmax) or fmax <= 0:
        fmin, fmax = 0.0, 1.0
    size_min, size_max = 4.0, 16.0

    for gy, g in enumerate(groups):
        for gx, gene in enumerate(genes):
            xs.append(gene); ys.append(g)
            f = float(df_frac.iat[gy, gx]); m = float(df_mean.iat[gy, gx])
            frac_norm = 0.0 if fmax == fmin else (f - fmin) / (fmax - fmin)
            sizes.append(size_min + frac_norm * (size_max - size_min))
            colors.append(m)
            hover_txt.append(f"group={g}<br>gene={gene}")

    fig = make_subplots(rows=2, cols=1, shared_xaxes=False, vertical_spacing=0.02,
                        row_heights=[0.973, 0.027])

    fig.add_trace(
        go.Scatter(
            x=xs, y=ys, mode="markers",
            marker=dict(size=sizes, sizemode="diameter", color=colors, colorscale="Magma",
                        showscale=True, colorbar=dict(title="mean")),
            text=hover_txt, hovertemplate="%{text}<extra></extra>", showlegend=False
        ),
        row=1, col=1
    )
    for v in [0.25, 0.5, 0.75]:
        fig.add_trace(
            go.Scatter(x=[None], y=[None], mode="markers",
                       marker=dict(size=size_min + v*(size_max - size_min), color="lightgray"),
                       showlegend=True, name=f"frac {v:g}"),
            row=1, col=1
        )

    strip_scale = _build_discrete_colorscale(palette_colors if len(palette_colors) >= 2 else ["#9aa0a6", "#9aa0a6"])
    z = np.array([strip_vals], dtype=float)
    fig.add_trace(
        go.Heatmap(
            z=z, x=genes, y=["type"], zmin=0, zmax=max(1, len(observed_types)-1),
            colorscale=strip_scale, showscale=False
        ),
        row=2, col=1
    )

    fig.update_xaxes(row=1, col=1, side="top", tickangle=70, automargin=True,
                     ticklabelposition="outside top", tickfont=dict(size=11), showgrid=True, showticklabels=True)
    fig.update_yaxes(row=1, col=1, autorange="reversed", categoryorder="array", categoryarray=groups,
                     tickmode="array", tickvals=groups, ticktext=groups, automargin=True, tickfont=dict(size=12),
                     showgrid=True)
    fig.update_xaxes(row=2, col=1, showticklabels=False)
    fig.update_yaxes(row=2, col=1, visible=False)

    runs = _run_lengths(gene_types)
    for s, e, t in runs:
        if e < s: 
            continue
        cx = (s + e) / 2.0
        fig.add_annotation(
            xref="x2", yref="y2",
            x=genes[int(round(cx))], y="type",
            text=f"<b>{t}</b>",
            showarrow=False, font=dict(size=12),
            yshift=0
        )

    height_px = max(560, int(26 * len(groups) + 240))
    fig.update_layout(template="plotly_white", hovermode="closest",
                      margin=dict(l=10, r=10, t=180, b=8),
                      legend=dict(orientation="h", yanchor="bottom", y=-0.18, font=dict(size=12)),
                      height=height_px, font=dict(size=13))
    title_outside(fig, title, top_margin=180)
    return p1.to_html(fig, full_html=False)

# ---- Main ----

def main():
    ap = argparse.ArgumentParser(description="Build fully-embedded FinalReport.html from pipeline outputs (v4.8).")
    ap.add_argument("-b", "--base", required=True, help="BASE_DIR path")
    ap.add_argument("-o", "--out", default=None, help="Output directory")
    ap.add_argument("-r", "--marker", default=None, help="Marker CSV (columns: type + gene_name/gene_id/gene)")
    ap.add_argument("-n", "--name", default="FinalReport.html", help="Output HTML file name (default: FinalReport.html)")
    args = ap.parse_args()

    base_dir = Path(args.base).resolve()
    paths = resolve_paths(base_dir)

    out_dir = Path(args.out).resolve() if args.out else paths["default_report_dir"]
    ensure_dir(out_dir)
    out_path = out_dir / args.name
    project = base_dir.name

    filter_txt = read_text(paths["filter_summary"])
    tss_cut, umi_pair, df_A, df_B = parse_filter_summary(filter_txt or "")
    filter_cards_html = f"""
<div class="cards">
  <div class="card"><div class="k">TSS cutoff</div><div class="v">{f"{tss_cut:.1f}" if tss_cut is not None else "NA"}</div></div>
  <div class="card"><div class="k">UMI min/max</div><div class="v">{(fmt_int(umi_pair[0])+' / '+fmt_int(umi_pair[1])) if umi_pair else "NA"}</div></div>
  <div class="card"><div class="k">—</div><div class="v"> </div></div>
</div>
"""
    filter_tables_html = ""
    filter_tables_html += df_to_html_table(df_A, "A) Filtering (TSSE/UMI)") if df_A is not None else '<div class="muted">No A) Filtering block.</div>'
    filter_tables_html += df_to_html_table(df_B, "B) Doublet removal (Scrublet)") if df_B is not None else '<div class="muted">No B) Doublet removal block.</div>'

    # Clustering!!
    umap_df, umap_errs = load_umap_csv(paths["umap_csv"])
    cluster_cards_html = ""; s_counts = pd.DataFrame(); l_counts = pd.DataFrame()
    cluster_figs_html = ""
    if umap_errs:
        cluster_figs_html += "<div class='muted'>UMAP CSV errors:</div><ul>" + "".join([f"<li>{e}</li>" for e in umap_errs]) + "</ul>"
    if not umap_df.empty:
        total = len(umap_df); n_samples = umap_df["sample"].nunique(); n_leiden = umap_df["leiden"].nunique()
        cluster_cards_html = f"""
<div class="cards">
  <div class="card"><div class="k">Total cells</div><div class="v">{fmt_int(total)}</div></div>
  <div class="card"><div class="k">Samples</div><div class="v">{fmt_int(n_samples)}</div></div>
  <div class="card"><div class="k">Leiden clusters</div><div class="v">{fmt_int(n_leiden)}</div></div>
</div>
"""
        s_counts = umap_df["sample"].value_counts().reset_index()
        s_counts.columns = ["sample", "cells"]
        l_counts = umap_df["leiden"].value_counts().reset_index()
        l_counts.columns = ["leiden", "cells"]

        p1 = PlotlyOnce()
        fig_depth  = make_umap_by_depth(umap_df, "UMAP by depth (log10 fragments)", p1, fig_size=560)
        fig_sample = make_umap_grouped(umap_df, "sample", "UMAP by sample", p1, fig_size=560)
        fig_leiden = make_umap_grouped(umap_df, "leiden", "UMAP by leiden", p1, fig_size=560)
        cluster_figs_html += f'<div class="grid-1x3"><div class="figure">{fig_depth}</div><div class="figure">{fig_sample}</div><div class="figure">{fig_leiden}</div></div>'
    else:
        cluster_figs_html += "<div class='muted'>No UMAP data found in umap.csv.</div>"

    left_col  = filter_cards_html + filter_tables_html
    right_col = cluster_cards_html
    if not s_counts.empty:
        right_col += df_to_html_table(s_counts, "Cells per sample")
    if not l_counts.empty:
        right_col += df_to_html_table(l_counts, "Cells per leiden", extra_class="narrow tall")
    pair_html = f'<div class="grid-1x2"><div class="vcol">{left_col}</div><div class="vcol">{right_col}</div></div>'
    sec_filter = wrap_section("filter", "Filter & Clustering", "Cards are width-aligned with their tables.", pair_html + cluster_figs_html)

    # palette
    lm = None
    try:
        lm = pd.read_csv(paths["label_map"])
        lm.columns = [c.lower() for c in lm.columns]
        lm = lm[["leiden", "type"]].astype(str)
    except Exception:
        lm = None

    type_levels_all: List[str] = []
    if lm is not None:
        type_levels_all += lm["type"].tolist()

    marker_types_all: List[str] = []
    gene_type_map: Optional[Dict[str, str]] = None
    marker_csv_path = Path(args.marker) if args.marker else None
    df_frac_tmp, df_mean_tmp = None, None
    try:
        df_frac_tmp = pd.read_csv(paths["dot_frac"], index_col=0)
        df_mean_tmp = pd.read_csv(paths["dot_mean"], index_col=0)
    except Exception:
        pass
    genes_for_map = list(df_frac_tmp.columns.astype(str)) if isinstance(df_frac_tmp, pd.DataFrame) else []
    if marker_csv_path and genes_for_map:
        gene_type_map = parse_marker_types(marker_csv_path, genes_for_map)
        if gene_type_map:
            marker_types_all += list(dict.fromkeys(gene_type_map.values()))

    type_levels = list(dict.fromkeys(type_levels_all + marker_types_all + (["NA"] if True else [])))
    type_palette_map = make_type_palette(type_levels)

    # Annotation
    annot_inner = ""
    p2 = PlotlyOnce()

    df_frac, df_mean = None, None
    try:
        df_frac, df_mean = pd.read_csv(paths["dot_frac"], index_col=0), pd.read_csv(paths["dot_mean"], index_col=0)
        rows = [r for r in df_frac.index if r in set(df_mean.index)]
        cols = [c for c in df_frac.columns if c in set(df_mean.columns)]
        df_frac, df_mean = df_frac.loc[rows, cols], df_mean.loc[rows, cols]
    except Exception:
        df_frac, df_mean = None, None

    if df_frac is not None and df_mean is not None:
        dot_html = dotplot_with_typestrip_bottom(
            df_frac, df_mean,
            "Dotplot (size=frac, color=mean) + type strip",
            p2,
            gene_type_map=gene_type_map,
            type_palette=type_palette_map
        )
    else:
        dot_html = "<div class='muted'>Missing or mismatched dot_frac/mean CSVs</div>"

    mid_col = ""
    if lm is not None:
        mid_col += df_to_html_table(lm, "label_map (leiden → type)")
    else:
        mid_col += "<div class='muted'>No label_map.csv found.</div>"

    right_right = ""
    if lm is not None and not umap_df.empty:
        umap_with_type = umap_df.merge(lm, how="left", on="leiden")
        umap_with_type["type"] = umap_with_type["type"].fillna("NA").astype(str)
        umap_type_html = make_umap_grouped(umap_with_type, "type", "UMAP by type (label_map)", p2,
                                           palette_map=type_palette_map, fig_size=560)
        right_right += f'<div class="figure">{umap_type_html}</div>'
    else:
        right_right += "<div class='muted'>UMAP by type unavailable.</div>"

    annot_inner += f'<div class="grid-1x3"><div class="figure">{dot_html}</div><div class="vcol">{mid_col}</div><div class="figure">{right_right}</div></div>'
    sec_annotation = wrap_section("annotation", "Annotation", None, annot_inner)

    # update_h5ad
    update_block = ""
    if paths["update_summary"].exists():
        try:
            upd = pd.read_csv(paths["update_summary"])
            show = upd.copy()
            if show.shape[0] > 30: show = show.head(30)
            if show.shape[1] > 18: show = show.iloc[:, :18]
            update_block = '<div class="block">' + df_to_html_table(show, "update_h5ad summary (preview)", extra_class="wide") + '</div>'
        except Exception:
            update_block = "<div class='muted'>Failed to read update_summary.csv</div>"
    else:
        update_block = "<div class='muted'>No update_summary.csv found.</div>"
    sec_update = wrap_section("update", "update_h5ad", "Full-width preview", update_block)

    html = assemble_html(project=project, base_dir=base_dir, sections=[sec_filter, sec_annotation, sec_update])
    ensure_dir(out_path.parent)
    out_path.write_text(html, encoding="utf-8")
    print(f"[DONE] Saved: {out_path}")

if __name__ == "__main__":
    main()
