#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch

class GridFlow:
    def __init__(self, figsize=(18, 9.5), cell_w=3.2, cell_h=1.25, gap_x=1.3, gap_y=1.4, font=10):
        self.fig, self.ax = plt.subplots(figsize=figsize)
        self.cell_w = cell_w
        self.cell_h = cell_h
        self.gx = cell_w + gap_x
        self.gy = cell_h + gap_y
        self.font = font
        self.nodes = {}

    def _rc_to_xy(self, row, col, w=None, h=None):
        """Bottom-left (x,y) from (row,col). Row 0 at top (invert y)."""
        if w is None: w = self.cell_w
        if h is None: h = self.cell_h
        x = col * self.gx
        y = -row * self.gy
        return x, y, w, h

    def _center(self, x, y, w, h):
        return x + w/2, y + h/2

    def add_node(self, key, label, row, col, shape="rect", w=None, h=None, fontsize=None):
        x, y, w, h = self._rc_to_xy(row, col, w, h)
        artist = None

        if shape == "rect":
            artist = Rectangle((x, y), w, h, fill=False, linewidth=1.6)
            self.ax.add_patch(artist)

        elif shape == "para":
            skew = 0.45
            pts = [(x+skew, y), (x+w, y), (x+w-skew, y+h), (x, y+h)]
            artist = Polygon(pts, fill=False, linewidth=1.6)
            self.ax.add_patch(artist)

        elif shape == "diamond":
            cx, cy = self._center(x, y, w, h)
            pts = [(cx, cy+h/2), (cx+w/2, cy), (cx, cy-h/2), (cx-w/2, cy)]
            artist = Polygon(pts, fill=False, linewidth=1.8)
            self.ax.add_patch(artist)

        else:
            raise ValueError("Unknown shape: choose from 'rect', 'para', 'diamond'.")

        cx, cy = self._center(x, y, w, h)
        self.ax.text(cx, cy, label, ha="center", va="center",
                     fontsize=(fontsize or self.font))
        self.nodes[key] = {"x": x, "y": y, "w": w, "h": h, "shape": shape, "artist": artist}
        return key

    def _east(self, nd):   return (nd["x"]+nd["w"],   nd["y"]+nd["h"]/2)
    def _west(self, nd):   return (nd["x"],           nd["y"]+nd["h"]/2)
    def _north(self, nd):  return (nd["x"]+nd["w"]/2, nd["y"]+nd["h"])
    def _south(self, nd):  return (nd["x"]+nd["w"]/2, nd["y"])

    def _pad_point(self, p1, p2, pad=0.12):
        x1, y1 = p1; x2, y2 = p2
        dx, dy = x2 - x1, y2 - y1
        mag = max((dx*dx + dy*dy) ** 0.5, 1e-6)
        return (x1 + pad*dx/mag, y1 + pad*dy/mag)

    def connect(self, src_key, dst_key, label=None,
                prefer="auto", style="solid", bend=0.0, pad=0.18,
                dx=0.0, dy=0.45, label_box=True):
        s = self.nodes[src_key]; d = self.nodes[dst_key]

        same_row = abs((s["y"]+s["h"]/2) - (d["y"]+d["h"]/2)) < self.gy/2
        same_col = abs((s["x"]+s["w"]/2) - (d["x"]+d["w"]/2)) < self.gx/2

        if prefer == "h" or (prefer == "auto" and same_row and not same_col):
            p1 = self._east(s); p2 = self._west(d)
        elif prefer == "v" or (prefer == "auto" and same_col):
            if s["y"] > d["y"]:
                p1 = self._south(s); p2 = self._north(d)
            else:
                p1 = self._north(s); p2 = self._south(d)
        else:
            p1 = self._east(s); p2 = self._west(d)

        p1p = self._pad_point(p1, p2, pad=pad)
        p2p = self._pad_point(p2, p1, pad=pad)

        ls = "-" if style == "solid" else "--"
        conn = f"arc3,rad={bend:.3f}" if bend else "arc3"

        arr = FancyArrowPatch(p1p, p2p, arrowstyle="->", mutation_scale=12,
                              linewidth=1.4, linestyle=ls,
                              connectionstyle=conn)
        self.ax.add_patch(arr)

        if label:
            mx, my = (p1p[0]+p2p[0])/2 + dx, (p1p[1]+p2p[1])/2 + dy
            if label_box:
                bbox = dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.2)
            else:
                bbox = None
            self.ax.text(mx, my, label, ha="center", va="bottom",
                         fontsize=self.font-1, bbox=bbox)

    def legend(self, x=0.2, y=1.12):
        self.ax.text(x, y, "Legend:", fontsize=self.font+1)
        self.ax.text(x, y-0.35,  "— solid: data flow", fontsize=self.font)
        self.ax.text(x, y-0.70,  "-- dashed: dependency", fontsize=self.font)

    def show(self, xpad=1.5, ypad=1.5):
        if not self.nodes:
            self.ax.axis("off"); plt.show(); return
        xs = []; ys = []; xe = []; ye = []
        for nd in self.nodes.values():
            xs.append(nd["x"]); ys.append(nd["y"])
            xe.append(nd["x"]+nd["w"]); ye.append(nd["y"]+nd["h"])
        self.ax.set_xlim(min(xs)-xpad, max(xe)+xpad)
        self.ax.set_ylim(min(ys)-ypad, max(ye)+ypad)
        self.ax.axis("off")
        self.fig.subplots_adjust(left=0.03, right=0.98, top=0.92, bottom=0.06)
        plt.show()


if __name__ == "__main__":
    gf = GridFlow(figsize=(18, 9.5), cell_w=3.2, cell_h=1.25, gap_x=1.3, gap_y=1.4, font=10)

    gf.ax.text(1.6, 0.6,  "Inputs",     fontsize=12, ha="center")
    gf.ax.text(6.0, 0.6,  "Pre-flight",  fontsize=12, ha="center")
    gf.ax.text(10.8, 0.6, "Processing",  fontsize=12, ha="center")
    gf.ax.text(15.6, 0.6, "Outputs",     fontsize=12, ha="center")

    gf.add_node("base", "base_dir\n<sample>/outs/fragments.tsv.gz", row=1, col=0, shape="para")
    gf.add_node("csz",  "chrom.sizes (rice)",                       row=2, col=0, shape="para")
    gf.add_node("gtf",  "genes.gtf (rice, unzipped)",               row=3, col=0, shape="para")

    gf.add_node("collect", "Collect fragments",         row=1, col=2, shape="rect")
    gf.add_node("loadcs",  "Load chrom.sizes → dict",   row=2, col=2, shape="rect")
    gf.add_node("checkgtf","Check GTF (.gtf, not .gz)", row=3, col=2, shape="rect")

    gf.add_node("fsd_loop", "Per sample:\nimport_fragments + frag_size_distr", row=1, col=4, shape="rect")
    gf.add_node("fsd_grid", "Fragment-size grid (PNG/SVG)",                    row=1, col=5, shape="rect")

    gf.add_node("tsse_loop","Per sample:\nTSSE compute (GTF wrapper)\ncollect valid TSSE", row=2, col=4, shape="rect")
    gf.add_node("tsse_dec", "Any valid TSSE?",                                  row=2, col=5, shape="diamond")
    gf.add_node("violin",   "TSSE violin (PNG/SVG)",                             row=2, col=6, shape="rect")
    gf.add_node("skip_vio", "SKIP violin",                                      row=3, col=6, shape="rect")

    gf.add_node("grid_tsse","Per sample:\nTSSE figure (Plotly→image or fallback)\ncombine to grid", row=3, col=4, shape="rect")
    gf.add_node("tsse_grid","TSSE combined grid (PNG/SVG)",                                         row=3, col=6, shape="rect")

    gf.add_node("o_fsd",   "fragment_size.png/svg", row=0, col=7, shape="para")
    gf.add_node("o_violin","TSSE_violin.png/svg",   row=1, col=7, shape="para")
    gf.add_node("o_grid",  "TSSE_grid.png/svg",     row=2, col=7, shape="para")
    gf.add_node("o_sum",   "qc_summary.txt / .csv", row=3, col=7, shape="para")

    gf.connect("base",   "collect",  prefer="h", style="solid", dy=0.55)
    gf.connect("csz",    "loadcs",   prefer="h", style="solid", dy=0.55)
    gf.connect("gtf",    "checkgtf", prefer="h", style="solid", dy=0.55)

    gf.connect("collect","fsd_loop", prefer="h", style="solid", dy=0.55)
    gf.connect("loadcs", "fsd_loop", prefer="v", style="dashed", bend=0.08, dy=0.35)
    gf.connect("fsd_loop","fsd_grid",prefer="h", style="solid", dy=0.55)
    gf.connect("fsd_grid","o_fsd",   prefer="h", style="solid", dy=0.55)

    gf.connect("collect","tsse_loop", prefer="v", style="solid", bend=0.06, dy=0.30)
    gf.connect("loadcs", "tsse_loop", prefer="h", style="dashed", dy=0.55)
    gf.connect("checkgtf","tsse_loop",prefer="h", style="dashed", dy=0.35)

    gf.connect("tsse_loop","tsse_dec", prefer="h", style="solid", dy=0.55)
    gf.connect("tsse_dec","violin",    prefer="h", style="solid", dy=0.60)
    gf.connect("tsse_dec","skip_vio",  prefer="v", style="solid", dy=0.30)
    gf.connect("violin",  "o_violin",  prefer="h", style="solid", dy=0.55)

    gf.connect("tsse_loop","grid_tsse", prefer="v", style="solid", bend=0.08, dy=0.28)
    gf.connect("grid_tsse","tsse_grid", prefer="h", style="solid", dy=0.60)
    gf.connect("tsse_grid","o_grid",    prefer="h", style="solid", dy=0.55)

    gf.add_node("summary","QC summary\n(build per-sample table,\nwrite TXT/CSV)", row=4, col=6, shape="rect")
    gf.connect("fsd_grid","summary",   prefer="v", style="dashed", bend=0.10, dy=0.30)
    gf.connect("tsse_loop","summary",  prefer="v", style="dashed", bend=0.05, dy=0.30)
    gf.connect("summary","o_sum",      prefer="h", style="solid", dy=0.55)

    gf.legend(x=0.2, y=1.12)
    gf.show()
