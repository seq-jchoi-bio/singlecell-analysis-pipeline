#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch

class GridFlow:
    def __init__(self, figsize=(18, 9), cell_w=3.2, cell_h=1.25, gap_x=1.3, gap_y=1.4, font=10):
        self.fig, self.ax = plt.subplots(figsize=figsize)
        self.cw, self.ch = cell_w, cell_h
        self.gx, self.gy = cell_w + gap_x, cell_h + gap_y
        self.font = font
        self.nodes = {}

    def _xy(self, r, c, w=None, h=None):
        w = self.cw if w is None else w; h = self.ch if h is None else h
        return c*self.gx, -r*self.gy, w, h

    def _center(self, x, y, w, h): return x+w/2, y+h/2

    def add(self, key, label, r, c, shape="rect", w=None, h=None, fs=None):
        x, y, w, h = self._xy(r, c, w, h)
        if shape == "rect":
            artist = Rectangle((x,y), w,h, fill=False, linewidth=1.6)
        elif shape == "para":
            skew=0.45; pts=[(x+skew,y),(x+w,y),(x+w-skew,y+h),(x,y+h)]
            artist = Polygon(pts, fill=False, linewidth=1.6)
        elif shape == "diamond":
            cx,cy=self._center(x,y,w,h); pts=[(cx,cy+h/2),(cx+w/2,cy),(cx,cy-h/2),(cx-w/2,cy)]
            artist = Polygon(pts, fill=False, linewidth=1.8)
        else:
            raise ValueError("shape must be rect|para|diamond")
        self.ax.add_patch(artist)
        cx, cy = self._center(x,y,w,h)
        self.ax.text(cx, cy, label, ha="center", va="center", fontsize=(fs or self.font))
        self.nodes[key] = {"x":x,"y":y,"w":w,"h":h}
        return key

    def _east(self,n):  return (n["x"]+n["w"],   n["y"]+n["h"]/2)
    def _west(self,n):  return (n["x"],         n["y"]+n["h"]/2)
    def _north(self,n): return (n["x"]+n["w"]/2, n["y"]+n["h"])
    def _south(self,n): return (n["x"]+n["w"]/2, n["y"])

    def _pad(self, p1, p2, pad=0.18):
        x1,y1=p1; x2,y2=p2; dx,dy=x2-x1,y2-y1
        L=max((dx*dx+dy*dy)**0.5,1e-6)
        return (x1+pad*dx/L, y1+pad*dy/L)

    def connect(self, src, dst, style="solid", prefer="auto", bend=0.06, pad=0.18, dx=0.0, dy=0.45, label=None):
        s,d=self.nodes[src], self.nodes[dst]
        same_row=abs((s["y"]+s["h"]/2)-(d["y"]+d["h"]/2))<self.gy/2
        same_col=abs((s["x"]+s["w"]/2)-(d["x"]+d["w"]/2))<self.gx/2
        if prefer=="h" or (prefer=="auto" and same_row and not same_col):
            p1=self._east(s); p2=self._west(d)
        elif prefer=="v" or (prefer=="auto" and same_col):
            p1,p2=(self._south(s),self._north(d)) if s["y"]>d["y"] else (self._north(s),self._south(d))
        else:
            p1=self._east(s); p2=self._west(d)
        p1p=self._pad(p1,p2,pad); p2p=self._pad(p2,p1,pad)
        ls="-" if style=="solid" else "--"
        conn=f"arc3,rad={bend:.3f}" if bend else "arc3"
        arr=FancyArrowPatch(p1p,p2p,arrowstyle="->",mutation_scale=12,linewidth=1.4,linestyle=ls,connectionstyle=conn)
        self.ax.add_patch(arr)
        if label:
            mx,my=(p1p[0]+p2p[0])/2+dx,(p1p[1]+p2p[1])/2+dy
            self.ax.text(mx,my,label,ha="center",va="bottom",fontsize=self.font-1,
                         bbox=dict(facecolor="white",edgecolor="none",alpha=0.75,pad=0.2))

    def legend(self, x=0.25, y=1.10):
        self.ax.text(x,y,"Legend:",fontsize=self.font+1)
        self.ax.text(x,y-0.35,"— solid: data flow",fontsize=self.font)
        self.ax.text(x,y-0.70,"-- dashed: dependency",fontsize=self.font)

    def show(self, xpad=1.5, ypad=1.6):
        if not self.nodes: self.ax.axis("off"); plt.show(); return
        xs,ys,xe,ye=[],[],[],[]
        for nd in self.nodes.values():
            xs.append(nd["x"]); ys.append(nd["y"]); xe.append(nd["x"]+nd["w"]); ye.append(nd["y"]+nd["h"])
        self.ax.set_xlim(min(xs)-xpad, max(xe)+xpad)
        self.ax.set_ylim(min(ys)-ypad, max(ye)+ypad)
        self.ax.axis("off")
        self.fig.subplots_adjust(left=0.03,right=0.98,top=0.92,bottom=0.07)
        plt.show()

# Diagram
if __name__ == "__main__":
    gf = GridFlow(figsize=(18,9), cell_w=3.2, cell_h=1.25, gap_x=1.3, gap_y=1.4, font=10)

    gf.ax.text(1.6, 0.6,  "Inputs",                 fontsize=12, ha="center")
    gf.ax.text(7.0, 0.6,  "Load & Preprocess",      fontsize=12, ha="center")
    gf.ax.text(12.4, 0.6, "Gene Activity & Labels", fontsize=12, ha="center")
    gf.ax.text(16.6, 0.6, "Outputs",                fontsize=12, ha="center")

    gf.add("in_h5",  "merged_doublets.h5ads\n(read_dataset)",  1, 0, "para")
    gf.add("in_csz", "chrom.sizes\n(load_chrom_sizes)",        2, 0, "para")
    gf.add("in_gtf", "genes.gtf\n(build_gene_only_gtf_from_any)", 3, 0, "para")
    gf.add("in_csv", "marker CSV\n(load_reference_* )",        4, 0, "para")

    gf.add("load", "Load pre-merged dataset\n(read_dataset)",  1, 2, "rect")
    gf.add("prep", "Select→Spectral→Harmony→KNN→Leiden\n(select_features / spectral / harmony / knn / leiden)", 2, 3, "rect", w=5.8, h=1.25)
    gf.add("umap", "UMAP + plots (optional)\n(umap / plot_umap_matplotlib / plot_umap_with_text)", 3, 3, "rect", w=5.8, h=1.25)

    gf.add("ga",    "Gene activity (auto+fallback)\n(make_gene_matrix_geneid_only / build_gene_activity_fallback)", 2, 6, "rect", w=5.8)
    gf.add("sig",   "Load markers / weights\n(load_reference_markers / load_reference_weights)",                     3, 6, "rect", w=5.8)
    gf.add("score", "Score\n(score_auc_rank / score_weighted_sum)",                                                 4, 6, "rect", w=5.8)
    gf.add("pick",  "Pick labels + confidence\n(pick_labels)",                                                       5, 6, "rect", w=5.8)
    gf.add("post",  "Map clusterName & Leiden-majority\n(build_clustername_mapping / add_leiden_majority_label)",    6, 6, "rect", w=5.8)

    gf.add("o_h5",  "annot_merged_cells.h5ad\n(.write)",             2, 9, "para")
    gf.add("o_csv", "annotation_results.csv\n(to_csv)",              3, 9, "para")
    gf.add("o_bed", "BASE_DIR/bed/*.bed\n(write_bed_per_leiden_and_sample)", 4, 9, "para")
    gf.add("o_sum", "summary.txt\n(write_summary)",                  5, 9, "para")

    gf.connect("in_h5","load", prefer="h", dy=0.55)
    gf.connect("load","prep",  prefer="h", dy=0.55)
    gf.connect("prep","umap",  prefer="v", dy=0.30)

    gf.connect("in_csz","ga",  style="dashed", prefer="h", dy=0.40)
    gf.connect("in_gtf","ga",  style="dashed", prefer="h", dy=0.15)
    gf.connect("load","ga",    prefer="h", dy=0.55)

    gf.connect("in_csv","sig", prefer="h", dy=0.55)
    gf.connect("ga","score",   prefer="v", dy=0.30)
    gf.connect("sig","score",  prefer="v", dy=0.30, bend=-0.06)
    gf.connect("score","pick", prefer="v", dy=0.30)
    gf.connect("pick","post",  prefer="v", dy=0.30)

    gf.connect("post","o_h5",  prefer="h", dy=0.55)
    gf.connect("post","o_csv", prefer="h", dy=0.20, bend=-0.05)
    gf.connect("post","o_bed", prefer="h", dy=-0.15, bend=-0.06)
    gf.connect("post","o_sum", prefer="h", dy=-0.45, bend=-0.08)

    gf.legend(x=0.25, y=1.10)
    gf.show()
