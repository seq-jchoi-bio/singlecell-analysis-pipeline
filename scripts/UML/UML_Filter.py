#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch

class GridFlow:
    def __init__(self, figsize=(18, 9.5), cell_w=3.3, cell_h=1.15, gap_x=1.15, gap_y=1.25, font=10):
        self.fig, self.ax = plt.subplots(figsize=figsize)
        self.cell_w, self.cell_h = cell_w, cell_h
        self.gx, self.gy = cell_w + gap_x, cell_h + gap_y
        self.font, self.nodes = font, {}

    def _rc_to_xy(self, row, col, w=None, h=None):
        w = self.cell_w if w is None else w
        h = self.cell_h if h is None else h
        return col * self.gx, -row * self.gy, w, h

    def _center(self, x, y, w, h): return x + w/2, y + h/2

    def add_node(self, key, label, row, col, shape="rect", w=None, h=None):
        x, y, w, h = self._rc_to_xy(row, col, w, h)
        if shape == "rect":
            artist = Rectangle((x, y), w, h, fill=False, linewidth=1.6)
        elif shape == "para":  # I/O
            skew = 0.42; pts = [(x+skew, y), (x+w, y), (x+w-skew, y+h), (x, y+h)]
            artist = Polygon(pts, fill=False, linewidth=1.6)
        elif shape == "diamond":  # decision
            cx, cy = self._center(x, y, w, h)
            pts = [(cx, cy+h/2), (cx+w/2, cy), (cx, cy-h/2), (cx-w/2, cy)]
            artist = Polygon(pts, fill=False, linewidth=1.8)
        else:
            raise ValueError("shape must be 'rect' | 'para' | 'diamond'")
        self.ax.add_patch(artist)
        cx, cy = self._center(x, y, w, h)
        self.ax.text(cx, cy, label, ha="center", va="center", fontsize=self.font)
        self.nodes[key] = {"x": x, "y": y, "w": w, "h": h}
        return key

    def _east(self, nd): return (nd["x"]+nd["w"], nd["y"]+nd["h"]/2)
    def _west(self, nd): return (nd["x"], nd["y"]+nd["h"]/2)
    def _north(self, nd): return (nd["x"]+nd["w"]/2, nd["y"]+nd["h"])
    def _south(self, nd): return (nd["x"]+nd["w"]/2, nd["y"])

    def _pad_point(self, p1, p2, pad=0.16):
        x1, y1 = p1; x2, y2 = p2
        dx, dy = x2-x1, y2-y1
        L = max((dx*dx+dy*dy)**0.5, 1e-6)
        return (x1 + pad*dx/L, y1 + pad*dy/L)

    def connect(self, src, dst, style="solid", prefer="auto", bend=0.06, pad=0.18, dx=0.0, dy=0.38, label=None):
        s, d = self.nodes[src], self.nodes[dst]
        same_row = abs((s["y"]+s["h"]/2)-(d["y"]+d["h"]/2)) < self.gy/2
        same_col = abs((s["x"]+s["w"]/2)-(d["x"]+d["w"]/2)) < self.gx/2
        if prefer == "h" or (prefer == "auto" and same_row and not same_col):
            p1, p2 = self._east(s), self._west(d)
        elif prefer == "v" or (prefer == "auto" and same_col):
            p1, p2 = (self._south(s), self._north(d)) if s["y"] > d["y"] else (self._north(s), self._south(d))
        else:
            p1, p2 = self._east(s), self._west(d)
        p1p, p2p = self._pad_point(p1, p2, pad), self._pad_point(p2, p1, pad)
        ls = "-" if style == "solid" else "--"
        conn = f"arc3,rad={bend:.3f}" if bend else "arc3"
        self.ax.add_patch(FancyArrowPatch(p1p, p2p, arrowstyle="->", mutation_scale=12,
                                          linewidth=1.4, linestyle=ls, connectionstyle=conn))
        if label:
            mx, my = (p1p[0]+p2p[0])/2 + dx, (p1p[1]+p2p[1])/2 + dy
            self.ax.text(mx, my, label, fontsize=self.font-1, ha="center", va="bottom",
                         bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=0.2))

    def show(self, xpad=1.4, ypad=1.4, savepath=None, dpi=300):
        if self.nodes:
            xs=[]; ys=[]; xe=[]; ye=[]
            for nd in self.nodes.values():
                xs.append(nd["x"]); ys.append(nd["y"]); xe.append(nd["x"]+nd["w"]); ye.append(nd["y"]+nd["h"])
            self.ax.set_xlim(min(xs)-xpad, max(xe)+xpad)
            self.ax.set_ylim(min(ys)-ypad, max(ye)+ypad)
        self.ax.axis("off")
        self.fig.subplots_adjust(left=0.035, right=0.985, top=0.92, bottom=0.08)
        if savepath: self.fig.savefig(savepath, dpi=dpi)
        plt.show()

if __name__ == "__main__":
    gf = GridFlow(figsize=(18, 9.5), cell_w=3.3, cell_h=1.15, gap_x=1.15, gap_y=1.25, font=10)

    gf.ax.text(1.7, 0.6,  "Inputs",                 fontsize=12, ha="center")
    gf.ax.text(6.2, 0.6,  "Pre-flight",             fontsize=12, ha="center")
    gf.ax.text(11.1, 0.6, "Per-sample core",        fontsize=12, ha="center")
    gf.ax.text(16.3, 0.6, "Doublet & Outputs",      fontsize=12, ha="center")

    gf.add_node("base", "base_dir\n<sample>/outs/fragments.tsv.gz", row=1, col=0, shape="para")
    gf.add_node("csz",  "chrom.sizes",                               row=2, col=0, shape="para")
    gf.add_node("gtf",  "genes.gtf (plain .gtf)",                    row=3, col=0, shape="para")
    gf.add_node("cfg",  "Params:\nTSS_CUTOFF, UMI_MIN/MAX\nFEATURE_MODE, N_FEATURES\nSTRIP_BARCODE_SUFFIX", row=4, col=0, shape="rect", h=1.6)

    gf.add_node("finds",  "Find samples",            row=1, col=2, shape="rect")
    gf.add_node("loadcs", "Load chrom.sizes → dict", row=2, col=2, shape="rect")
    gf.add_node("chkchr", "Chrom naming sanity",     row=3, col=2, shape="rect")

    gf.add_node("import", "Import fragments (AnnData)",                           row=1, col=4, shape="rect")
    gf.add_node("tsse",   "Compute TSSE (GTF wrapper)",                           row=2, col=4, shape="rect")
    gf.add_node("qc",     "Build QC table\n(is__cell_barcode, passed_filters, TSSE)", row=3, col=4, shape="rect")
    gf.add_node("filt",   "Filter cells\nTSSE ≥ cutoff & UMI in [min,max]",       row=4, col=4, shape="rect")
    gf.add_node("csv",    "Save <sample>_filtered.csv",                            row=5, col=4, shape="para")
    gf.add_node("subset", "Subset kept barcodes\n(+ optional strip suffix)",       row=6, col=4, shape="rect", w=3.8)

    gf.add_node("tile",   "Add tile matrix\nEnsure CSR",                           row=1, col=6, shape="rect")

    gf.add_node("fmode",  "FEATURE_MODE ?",                                        row=2, col=6, shape="diamond")
    gf.add_node("common", "COMMON grid\n(no var subset)",                          row=3, col=6, shape="rect")
    gf.add_node("per",    "PER-SAMPLE select_features\n(n = N_FEATURES) → subset vars", row=4, col=6, shape="rect", w=4.6)

    gf.add_node("scrub",  "Scrublet → filter_doublets",                            row=5, col=6, shape="rect")
    gf.add_node("h5ad",   "Save <sample>_doublets.h5ad",                            row=6, col=6, shape="para")

    gf.add_node("out1", "Filtered CSVs",                               row=2, col=8, shape="para")
    gf.add_node("out2", "Doublet-free h5ad",                           row=3, col=8, shape="para")
    gf.add_node("out3", "summary.txt\n(params + pre/post + skipped)",  row=4, col=8, shape="para", h=1.5)

    gf.connect("base","finds",   prefer="h", dy=0.50)
    gf.connect("csz","loadcs",   prefer="h", dy=0.50)
    gf.connect("gtf","chkchr",   prefer="h", dy=0.50)
    gf.connect("cfg","finds",    style="dashed", prefer="h", dy=-0.25)

    gf.connect("finds","import", prefer="h", dy=0.50)
    gf.connect("loadcs","import",style="dashed", prefer="h", dy=0.20)
    gf.connect("chkchr","import",style="dashed", prefer="h", dy=-0.15)

    gf.connect("import","tsse",  prefer="h", dy=0.50)
    gf.connect("tsse","qc",      prefer="h", dy=0.50)
    gf.connect("qc","filt",      prefer="h", dy=0.50)
    gf.connect("filt","csv",     prefer="h", dy=0.50)
    gf.connect("csv","out1",     prefer="h", dy=0.50)

    gf.connect("filt","subset",  prefer="h", dy=0.50, label="if kept > 0", dx=-0.25)
    gf.connect("subset","tile",  prefer="v", dy=0.32)
    gf.connect("tile","fmode",   prefer="v", dy=0.32)

    gf.connect("fmode","common", prefer="v", dy=0.26, label="COMMON", dx=-0.55)
    gf.connect("fmode","per",    prefer="v", dy=0.26, label="PER-SAMPLE", dx=0.70)
    gf.connect("common","scrub", prefer="v", dy=0.32)
    gf.connect("per","scrub",    prefer="v", dy=0.32)
    gf.connect("scrub","h5ad",   prefer="v", dy=0.32)
    gf.connect("h5ad","out2",    prefer="h", dy=0.50)

    gf.connect("qc","out3",   style="dashed", prefer="h", dy=0.25, bend=0.08)
    gf.connect("scrub","out3",style="dashed", prefer="h", dy=-0.05, bend=-0.06)
    gf.connect("cfg","out3",  style="dashed", prefer="h", dy=-0.35, bend=0.10)

    gf.show(savepath=None)
