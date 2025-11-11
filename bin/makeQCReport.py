#!/usr/bin/env python3
from __future__ import annotations

import json, os
from pathlib import Path

PLOTLY_LOADER = """<script>
(function(){
  function ready(fn){ if (document.readyState !== 'loading'){ fn(); } else { document.addEventListener('DOMContentLoaded', fn); } }
  function loadPlotly(cb){
    if (window.Plotly) { cb(); return; }
    var s = document.createElement('script');
    s.src = 'https://cdn.plot.ly/plotly-2.27.0.min.js';
    s.onload = cb;
    document.head.appendChild(s);
  }
  window.__qcReady = function(cb){ ready(function(){ loadPlotly(cb); }); };
})();
</script>"""

def _read_json(p: Path):
    try:
        return json.loads(p.read_text(encoding='utf-8'))
    except Exception:
        return None

def _btn_row(btns):
    return '<div style="display:flex;gap:8px;flex-wrap:wrap">' + ''.join(btns) + '</div>'

def build_html_report(out_dir: str, html_title: str = 'QC Summary', embed_assets: bool = False) -> str:
    out = Path(out_dir)
    assets = out / 'QC_reports'
    (assets / 'assets' / 'css').mkdir(parents=True, exist_ok=True)
    (assets / 'assets' / 'js').mkdir(parents=True, exist_ok=True)
    (assets / 'assets' / 'img').mkdir(parents=True, exist_ok=True)
    (assets / 'data').mkdir(parents=True, exist_ok=True)

    txt = out / 'qc_summary.txt'
    csv = out / 'qc_summary_table.csv'
    if txt.exists():
        (assets / 'data' / txt.name).write_bytes(txt.read_bytes())
    if csv.exists():
        (assets / 'data' / csv.name).write_bytes(csv.read_bytes())

    css = (
        'body{font-family:system-ui,Segoe UI,Roboto,Apple SD Gothic Neo,AppleGothic,Malgun Gothic,Arial,sans-serif;margin:0;background:#f6f7fb;color:#111}'
        '.container{max-width:1280px;margin:24px auto;padding:0 16px}'
        '.card{background:#fff;border:1px solid #e5e7eb;border-radius:14px;padding:16px;margin:14px 0;box-shadow:0 1px 2px rgba(16,24,40,0.06)}'
        '.h1{font-size:22px;font-weight:700;margin:0 0 6px}'
        '.h2{font-size:18px;font-weight:700;margin:0 0 6px}'
        '.muted{color:#667085}'
        '.btn{display:inline-block;padding:8px 12px;border:1px solid #d1d5db;border-radius:8px;text-decoration:none;color:#111;background:#fff}'
        '.grid{display:grid;grid-template-columns:1fr;gap:12px}@media(min-width:980px){.grid{grid-template-columns:1fr 1fr}}'
        '.plot{width:100%;height:520px}'
    )
    (assets / 'assets' / 'css' / 'report.css').write_text(css, encoding='utf-8')

    def copy_img(stem: str):
        for ext in ('.png', '.svg'):
            p = out / (stem + ext)
            if p.exists():
                (assets / 'assets' / 'img' / p.name).write_bytes(p.read_bytes())
                return 'QC_reports/assets/img/' + p.name
        return None

    img_fragment = copy_img('fragment_size')
    img_tsseviolin = copy_img('TSSE_violin')
    img_tssegrid = copy_img('TSSE_grid')
    img_fviolin = copy_img('fragment_violin')
    img_knee = copy_img('knee_plot')

    js_fragment = _read_json(assets / 'data' / 'fragment_size.json')
    js_tv = _read_json(assets / 'data' / 'tsse_violin.json')
    js_grid = _read_json(assets / 'data' / 'tsse_grid.json')

    btns = []
    if csv.exists(): btns.append('<a class="btn" href="QC_reports/data/' + csv.name + '" download>Download CSV</a>')
    if txt.exists(): btns.append('<a class="btn" href="QC_reports/data/' + txt.name + '" download>Download TXT</a>')
    header = '<div class="card"><div class="h1">' + html_title + '</div><div class="muted">Interactive QC report</div>' + _btn_row(btns) + '</div>'

    summary_html = ''
    if txt.exists():
        t = txt.read_text(encoding='utf-8', errors='ignore')
        t = t.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
        summary_html = '<div class="card"><div class="h2">Summary (qc_summary.txt)</div><pre style="white-space:pre-wrap;background:#0b1021;color:#e3e7ff;padding:12px;border-radius:8px">' + t + '</pre></div>'

    plots_html = []
    scripts = []

    # 1) Fragment-size interactive
    if js_fragment and js_fragment.get('samples'):
        plots_html.append('<div class="card"><div class="h2">Fragment-size Distribution (interactive)</div><div id="plot_frag" class="plot"></div></div>')
        data_obj = json.dumps(js_fragment)
        scripts.append("""(function(){
  const payload = DATA;
  const traces = [];
  const xs = [...Array(payload.max_recorded_size).keys()].map(i => i+1);
  for (const s of payload.samples) {
    traces.push({x: xs, y: s.hist, type: 'scatter', mode: 'lines', name: s.name});
  }
  const layout = {xaxis:{title:'Fragment length (bp)'}, yaxis:{title:'Count'}, hovermode:'closest'};
  const config = {modeBarButtonsToAdd:[{name:'Download PNG (square)',title:'Download PNG (1024x1024)',icon:Plotly.Icons.camera,click:function(gd){var fn=(gd&&gd._fullLayout&&gd._fullLayout._meta&&gd._fullLayout._meta.filename)||(gd&&gd.layout&&gd.layout.title&&gd.layout.title.text)||'figure';Plotly.downloadImage(gd,{format:'png',width:1024,height:1024,filename:fn});}}], responsive:true, displaylogo:false, toImageButtonOptions:{width:1024,height:1024,format:'svg',filename:'fragment_size'}};
  Plotly.newPlot('plot_frag', traces, layout, config);
})();""".replace('DATA', data_obj))
    elif img_fragment:
        plots_html.append('<div class="card"><div class="h2">Fragment-size Distribution</div><img src="' + img_fragment + '" style="max-width:100%"></div>')

    # 2) TSSE violin interactive
    if js_tv and js_tv.get('labels'):
        plots_html.append('<div class="card"><div class="h2">TSSE (per-sample) Violin (interactive)</div><div id="plot_tsse_violin" class="plot"></div></div>')
        data_obj = json.dumps(js_tv)
        scripts.append("""(function(){
  const p = DATA;
  const traces = [];
  for (let i=0;i<p.labels.length;i++) {
    traces.push({type:'violin', y:p.tsse[i]||[], name:p.labels[i], box:{visible:false}, meanline:{visible:true}, points:false});
  }
  const shapes = [];
  if (p.global_median_tsse) {
    shapes.push({type:'line', xref:'paper', yref:'y', x0:0, x1:1, y0:p.global_median_tsse, y1:p.global_median_tsse, line:{dash:'dash', color:'red'}});
  }
  const layout = {xaxis:{automargin:true}, yaxis:{title:'TSSE'}, hovermode:'closest', shapes:shapes};
  const config = {modeBarButtonsToAdd:[{name:'Download PNG (square)',title:'Download PNG (1024x1024)',icon:Plotly.Icons.camera,click:function(gd){var fn=(gd&&gd._fullLayout&&gd._fullLayout._meta&&gd._fullLayout._meta.filename)||(gd&&gd.layout&&gd.layout.title&&gd.layout.title.text)||'figure';Plotly.downloadImage(gd,{format:'png',width:1024,height:1024,filename:fn});}}], responsive:true, displaylogo:false, toImageButtonOptions:{width:1024,height:1024,format:'svg',filename:'TSSE_violin'}};
  Plotly.newPlot('plot_tsse_violin', traces, layout, config);


// --- Global mean-of-medians toggle for plot_tsse_violin ---
(function(){
  const gd = document.getElementById('plot_tsse_violin');
  if (!gd) return;
  // UI
  const lbl = document.createElement('label');
  lbl.style.fontSize = '12px'; lbl.style.display = 'inline-block'; lbl.style.marginTop = '6px';
  const cb = document.createElement('input'); cb.type = 'checkbox'; cb.id='tv_global_toggle'; cb.checked = true;
  lbl.appendChild(cb); lbl.appendChild(document.createTextNode(' Show global median'));
  // place right under plot container
  if (gd.parentElement) gd.parentElement.appendChild(lbl);

  function _median(arr){
    const v = (arr||[]).filter(a => Number.isFinite(a)).slice().sort((a,b)=>a-b);
    const n = v.length; if (n===0) return null;
    return (n%2 ? v[(n-1)>>1] : 0.5*(v[n/2-1] + v[n/2]));
  }
  // Use current traces' y arrays
  let meds = [];
  if (gd.data && gd.data.length) {
    for (let i=0;i<gd.data.length;i++) {
      const y = gd.data[i].y || [];
      const m = _median(y);
      if (m !== null) meds.push(m);
    }
  }
  const gmm = meds.length ? (meds.reduce((a,b)=>a+b,0)/meds.length) : null;

  function apply(show){
    let shapes = (gd.layout && gd.layout.shapes) ? gd.layout.shapes.slice() : [];
    shapes = shapes.filter(s => s && s._tag !== 'global_median_tv');
    if (show && Number.isFinite(gmm)) {
      shapes.push({type:'line', xref:'paper', yref:'y', x0:0, x1:1, y0:gmm, y1:gmm,
                   line:{color:'red', width:2, dash:'dash'}, _tag:'global_median_tv'});
    }
    Plotly.relayout(gd, {'shapes': shapes});
  }
  apply(cb.checked);
  cb.addEventListener('change', function(){ apply(cb.checked); });
})();

})();""".replace('DATA', data_obj))
    elif img_tsseviolin:
        plots_html.append('<div class="card"><div class="h2">TSSE (per-sample) Violin</div><img src="' + img_tsseviolin + '" style="max-width:100%"></div>')

    # 3) TSSE grid (density) from prebinned json
    if js_grid and js_grid.get('samples'):
        grid_divs = []; grid_scripts = []
        # Optional mapping from label -> vectors for fragment sums
        label_to_idx = {}
        if js_tv and js_tv.get('labels'):
            for j, nm in enumerate(js_tv['labels']):
                label_to_idx[str(nm)] = j
        for idx, s in enumerate(js_grid['samples']):
            div_id = 'tg_' + str(idx)
            grid_divs.append('<div class="plot" id="' + div_id + '"></div>')
            sdat = json.dumps(s)
            rawdat = 'null'
            nm = str(s.get('name',''))
            if label_to_idx and nm in label_to_idx:
                j = label_to_idx[nm]
                rawdat = json.dumps({'nfrag': js_tv['nfrag'][j], 'tsse': js_tv['tsse'][j]})
            grid_scripts.append((r"""(function(){
  const s = S_DATA;
  const raw = RAW_DATA;
  const divId = 'DIV';
  const xe = s.x_edges, ye = s.y_edges;
  const x = xe.slice(0,-1), y = ye.slice(0,-1);
  const Z = s.z;
  // S: sum of fragments per bin (computed if raw vectors available)
  const S = Array.from({length:y.length}, ()=> Array(x.length).fill(0));
  if (raw && raw.nfrag && raw.tsse){
    const nf = raw.nfrag.filter(v=>isFinite(v)&&v>0);
    const ts = raw.tsse.filter(v=>isFinite(v));
    const n = Math.min(nf.length, ts.length);
    const findBin = (val, edges) => { let lo=0, hi=edges.length-1; while (lo<hi-1){ const mid=(lo+hi)>>1; if (val<edges[mid]) hi=mid; else lo=mid; } return Math.max(0, Math.min(edges.length-2, lo)); };
    for (let i=0;i<n;i++){ const xv=nf[i], yv=ts[i]; if (!isFinite(xv)||!isFinite(yv)||xv<=0) continue;
      const cx=findBin(xv, xe), cy=findBin(yv, ye); S[cy][cx]+=xv; }
  }
  
  // --- Area-corrected density to emulate PNG (per-bin) ---
  const xe_log = xe.map(v => Math.log10(v));
  const wx = xe_log.slice(1).map((v,i)=> Math.max(1e-12, v - xe_log[i]));
  const wy = ye.slice(1).map((v,i)=> Math.max(1e-12, v - ye[i]));
  let total = 0;
  for (let r=0;r<Z.length;r++){ for (let c=0;c<Z[0].length;c++){ const v=Z[r][c]; if (isFinite(v) && v>0) total += v; } }
  const ZD = [];
  for (let r=0;r<Z.length;r++){
    const row = new Array(Z[0].length);
    for (let c=0;c<Z[0].length;c++){
      const v = Z[r][c];
      row[c] = (isFinite(v)&&v>0 && total>0) ? (v / (total * wx[c] * wy[r])) : 0;
    }
    ZD.push(row);
  }

  // --- Gaussian smoothing (separable) ---
  function gaussianKernel1D(sigma){
    const rad = Math.max(1, Math.round(3*sigma));
    const k = new Array(2*rad+1);
    const s2 = 2*sigma*sigma;
    let sum = 0;
    for (let i=-rad;i<=rad;i++){ const w = Math.exp(-(i*i)/s2); k[i+rad] = w; sum += w; }
    for (let i=0;i<k.length;i++) k[i] /= sum;
    return {k:k, rad:rad};
  }
  function clamp(n, lo, hi){ return n<lo?lo:(n>hi?hi:n); }
  function convolveRows(M, ker){
    const H = M.length, W = M[0].length, out = new Array(H);
    for (let r=0;r<H;r++){
      const row = new Array(W);
      for (let c=0;c<W;c++){
        let acc = 0;
        for (let i=-ker.rad;i<=ker.rad;i++){
          const cc = clamp(c+i, 0, W-1);
          acc += ker.k[i+ker.rad] * M[r][cc];
        }
        row[c] = acc;
      }
      out[r] = row;
    }
    return out;
  }
  function convolveCols(M, ker){
    const H = M.length, W = M[0].length, out = new Array(H);
    for (let r=0;r<H;r++) out[r] = new Array(W).fill(0);
    for (let r=0;r<H;r++){
      for (let c=0;c<W;c++){
        let acc = 0;
        for (let j=-ker.rad;j<=ker.rad;j++){
          const rr = clamp(r+j, 0, H-1);
          acc += ker.k[j+ker.rad] * M[rr][c];
        }
        out[r][c] = acc;
      }
    }
    return out;
  }
  const kx = gaussianKernel1D(1.0);  // σx (bins)
  const ky = gaussianKernel1D(1.2);  // σy (bins)
  const ZR = convolveRows(ZD, kx);
  const ZS = convolveCols(ZR, ky);

  // --- Percentile clipping (P99) for better contrast with heavy tails ---
  function percentileFrom2D(M, q){
    const arr = [];
    for (let r=0;r<M.length;r++){
      for (let c=0;c<M[0].length;c++){
        const v = M[r][c];
        if (isFinite(v)) arr.push(v);
      }
    }
    if (arr.length === 0) return 0;
    arr.sort((a,b)=>a-b);
    const idx = Math.floor((q/100) * (arr.length-1));
    return arr[idx];
  }
  const zmax99 = percentileFrom2D(ZS, 99);

  const trace = {
    type:'heatmap',
    x:x, y:y, z:ZS,
    colorscale:'Aggrnyl',
    zmin:0,
    zmax:(zmax99>0?zmax99:undefined),
    colorbar:{title:'Density'}
  };

  const layout = { title:s.name, xaxis:{type:'log', title:'n_fragment'}, yaxis:{title:'TSSE'}, margin:{t:40}, dragmode:'select', shapes:[] };
  const config = {modeBarButtonsToAdd:[{name:'Download PNG (square)',title:'Download PNG (1024x1024)',icon:Plotly.Icons.camera,click:function(gd){var fn=(gd&&gd._fullLayout&&gd._fullLayout._meta&&gd._fullLayout._meta.filename)||(gd&&gd.layout&&gd.layout.title&&gd.layout.title.text)||'figure';Plotly.downloadImage(gd,{format:'png',width:1024,height:1024,filename:fn});}}],  responsive:true, displaylogo:false, toImageButtonOptions:{width:1024,height:1024,format:'svg',filename:'TSSE_grid_'+s.name} };
  Plotly.newPlot(divId, [trace], layout, config);

  
  // --- ESC clear (shapes + selections + selectedpoints), preserving select drag ---
  (function(){
    const gd = document.getElementById(divId);
    function clearAll(){
      try { Plotly.relayout(gd, {'shapes': [], 'selections': []}); } catch(e) {}
      try {
        const nTr = (gd && gd.data ? gd.data.length : 0);
        if (nTr > 0) Plotly.restyle(gd, {'selectedpoints': Array(nTr).fill(null)});
      } catch(e) {}
      try { if (panel) panel.textContent = ''; } catch(e) {}
      try { Plotly.relayout(gd, {'dragmode':'select'}); } catch(e) {}
    }
    document.addEventListener('keydown', function(e){
      if (e.key === 'Escape'){ clearAll(); }
    });
    if (gd && gd.on){
      gd.on('plotly_doubleclick', function(){ clearAll(); return false; });
    }
    if (typeof config !== 'undefined'){
      config.modeBarButtonsToAdd = (config.modeBarButtonsToAdd||[]).concat([{
        name:'Clear selection', title:'Clear selection / quadrants',
        icon:Plotly.Icons.trash, click:function(){ clearAll(); }
      }]);
    }
  })();
// Small info panel
  const panel = document.createElement('div');
  panel.style.fontSize = '12px'; panel.style.marginTop = '8px';
  document.getElementById(divId).parentElement.appendChild(panel);

  // Drag selection: sum Z (and S if available) within the selection rectangle using ev.range
  function summarizeSelection(range){
    if (!range || !range.x || !range.y) return;
    const x0 = Math.min(range.x[0], range.x[1]);
    const x1 = Math.max(range.x[0], range.x[1]);
    const y0 = Math.min(range.y[0], range.y[1]);
    const y1 = Math.max(range.y[0], range.y[1]);
    const rect = [{type:'rect', x0:x0, x1:x1, y0:y0, y1:y1, line:{width:2, dash:'dot'}, fillcolor:'rgba(0,0,0,0)'}];
    Plotly.relayout(divId, {'shapes': rect});
    let cellCount = 0, fragSum = 0;
    for (let r=0;r<y.length;r++){
      for (let c=0;c<x.length;c++){
        const xc0 = xe[c], xc1 = xe[c+1];
        const yc0 = ye[r], yc1 = ye[r+1];
        if (xc1>=x0 && xc0<=x1 && yc1>=y0 && yc0<=y1){
          cellCount += (Z[r][c]||0);
          fragSum += (S[r][c]||0);
        }
      }
    }
    panel.innerHTML = 'Selected: <b>' + cellCount.toLocaleString() + '</b> cells' +
      (fragSum > 0 ? ' | Σfrags=<b>' + Math.round(fragSum).toLocaleString() + '</b>' : '');
  }

  const gd = document.getElementById(divId);
  gd.on('plotly_selected', ev => {
    if (ev && ev.range) summarizeSelection(ev.range);
  });

  gd.on('plotly_click', ev => {
    if (!ev || !ev.points || !ev.points.length) return;
    const xp = ev.points[0].x, yp = ev.points[0].y;
    const shapes = [
      {type:'line', x0:xp, x1:xp, y0:ye[0], y1:ye[ye.length-1], line:{dash:'dot', width:2}},
      {type:'line', x0:xe[0], x1:xe[xe.length-1], y0:yp, y1:yp, line:{dash:'dot', width:2}}
    ];
    Plotly.relayout(divId, {'shapes': shapes});

    const findBin = (val, edges) => { let lo=0, hi=edges.length-1; while (lo<hi-1){ const mid=(lo+hi)>>1; if (val<edges[mid]) hi=mid; else lo=mid; } return Math.max(0, Math.min(edges.length-2, lo)); };
    const cx = findBin(xp, xe), cy = findBin(yp, ye);
    const sums = [0,0,0,0], sumsF = [0,0,0,0];
    for (let r=0;r<y.length;r++){
      for (let c=0;c<x.length;c++){
        const z = (Z[r] && Z[r][c]) || 0;
        const ssum = (S[r] && S[r][c]) || 0;
        if (c>=cx && r>=cy){ sums[0]+=z; sumsF[0]+=ssum; }      // Q1
        else if (c<cx && r>=cy){ sums[1]+=z; sumsF[1]+=ssum; }  // Q2
        else if (c<cx && r<cy){ sums[2]+=z; sumsF[2]+=ssum; }   // Q3
        else { sums[3]+=z; sumsF[3]+=ssum; }                    // Q4
      }
    }
    panel.innerHTML = 'Q1: ' + sums[0].toLocaleString() + (sumsF[0]>0 ? ' | Σfrags=' + Math.round(sumsF[0]).toLocaleString() : '') +
                      ' &nbsp; Q2: ' + sums[1].toLocaleString() + (sumsF[1]>0 ? ' | Σfrags=' + Math.round(sumsF[1]).toLocaleString() : '') +
                      '<br>Q3: ' + sums[2].toLocaleString() + (sumsF[2]>0 ? ' | Σfrags=' + Math.round(sumsF[2]).toLocaleString() : '') +
                      ' &nbsp; Q4: ' + sums[3].toLocaleString() + (sumsF[3]>0 ? ' | Σfrags=' + Math.round(sumsF[3]).toLocaleString() : '');
  });
})();""").replace('S_DATA', sdat).replace('RAW_DATA', rawdat).replace('DIV', div_id))
        plots_html.append('<div class="card"><div class="h2">TSSE Grid (pseudocolor, interactive)</div><div class="grid">' + ''.join(grid_divs) + '</div></div>')
        scripts.extend(grid_scripts)

    # 3b) TSSE grid built from raw vectors (fallback)
    elif js_tv and js_tv.get('labels') and js_tv.get('tsse') and js_tv.get('nfrag'):
        grid_divs = []; grid_scripts = []
        for idx, name in enumerate(js_tv['labels']):
            div_id = 'tg_' + str(idx)
            grid_divs.append('<div class="plot" id="' + div_id + '"></div>')
            vec_tsse = json.dumps(js_tv['tsse'][idx] or [])
            vec_nfrag = json.dumps(js_tv['nfrag'][idx] or [])
            nm = json.dumps(name)
            grid_scripts.append(r"""(function(){
  const ts = TS.filter(v=>isFinite(v)&&v>0);
  const nf = NF.filter(v=>isFinite(v)&&v>0);
  const name = NM;
  const bx = 60, by = 60;
  if (nf.length===0 || ts.length===0) { 
    const el = document.getElementById('DIV');
    el.innerHTML = '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:#667085">(no data)</div>'; 
    return; 
  }
  const minX = Math.max(1, Math.min(...nf));
  const maxX = Math.max(...nf);
  const minY = Math.min(...ts);
  const maxY = Math.max(...ts);
  function geomspace(a,b,n){ const r=[]; const la=Math.log(a), lb=Math.log(b), step=(lb-la)/n; for(let i=0;i<=n;i++) r.push(Math.exp(la+i*step)); return r; }
  let xedges = (minX>0 && maxX>minX) ? geomspace(minX, maxX, bx) : Array.from({length:bx+1}, (_,i)=>minX + (maxX-minX)*i/bx);
  let yedges = Array.from({length:by+1}, (_,i)=>minY + (maxY-minY)*i/by);
  const Z = Array.from({length:by}, ()=>Array(bx).fill(0));
  const S = Array.from({length:by}, ()=>Array(bx).fill(0));
  for (let i=0;i<nf.length;i++) {
    const x = nf[i], y = ts[i];
    if (!isFinite(x)||!isFinite(y)) continue;
    let xi = xedges.findIndex((v,idx)=> x>=v && x < xedges[idx+1]); if (xi<0) xi = bx-1;
    let yi = yedges.findIndex((v,idx)=> y>=v && y < yedges[idx+1]); if (yi<0) yi = by-1;
    if (xi>=0 && xi<bx && yi>=0 && yi<by){ Z[yi][xi] += 1; S[yi][xi] += x; }
  }
  const x = xedges.slice(0,-1);
  const y = yedges.slice(0,-1);
  
  // --- Area-corrected density to emulate PNG ---
  const xe_log = xe.map(v => Math.log10(v));
  const wx = xe_log.slice(1).map((v,i)=> Math.max(1e-12, v - xe_log[i]));
  const wy = ye.slice(1).map((v,i)=> Math.max(1e-12, v - ye[i]));
  let total = 0;
  for (let r=0;r<Z.length;r++){ for (let c=0;c<Z[0].length;c++){ const v=Z[r][c]; if (isFinite(v) && v>0) total += v; } }
  const ZD = [];
  for (let r=0;r<Z.length;r++){
    const row = new Array(Z[0].length);
    for (let c=0;c<Z[0].length;c++){
      const v = Z[r][c];
      row[c] = (isFinite(v)&&v>0 && total>0) ? (v / (total * wx[c] * wy[r])) : 0;
    }
    ZD.push(row);
  }

  // --- Gaussian smoothing (σx=1.0, σy=1.2 bins) ---
  function gaussianKernel1D(sigma){
    const rad = Math.max(1, Math.round(3*sigma));
    const k = new Array(2*rad+1);
    const s2 = 2*sigma*sigma;
    let sum = 0;
    for (let i=-rad;i<=rad;i++){ const w = Math.exp(-(i*i)/s2); k[i+rad] = w; sum += w; }
    for (let i=0;i<k.length;i++) k[i] /= sum;
    return {k:k, rad:rad};
  }
  function clamp(n, lo, hi){ return n<lo?lo:(n>hi?hi:n); }
  function convolveRows(M, ker){
    const H = M.length, W = M[0].length, out = new Array(H);
    for (let r=0;r<H;r++){
      const row = new Array(W);
      for (let c=0;c<W;c++){
        let acc = 0;
        for (let i=-ker.rad;i<=ker.rad;i++){
          const cc = clamp(c+i, 0, W-1);
          acc += ker.k[i+ker.rad] * M[r][cc];
        }
        row[c] = acc;
      }
      out[r] = row;
    }
    return out;
  }
  function convolveCols(M, ker){
    const H = M.length, W = M[0].length, out = new Array(H);
    for (let r=0;r<H;r++) out[r] = new Array(W).fill(0);
    for (let r=0;r<H;r++){
      for (let c=0;c<W;c++){
        let acc = 0;
        for (let j=-ker.rad;j<=ker.rad;j++){
          const rr = clamp(r+j, 0, H-1);
          acc += ker.k[j+ker.rad] * M[rr][c];
        }
        out[r][c] = acc;
      }
    }
    return out;
  }
  const kx = gaussianKernel1D(1.0);
  const ky = gaussianKernel1D(1.2);
  const ZR = convolveRows(ZD, kx);
  const ZS = convolveCols(ZR, ky);

  // --- Optional upsampling for higher visual resolution ---
  const UPSAMPLE = 3; // set to 1 to disable
  function upsampleRow(row, f){
    if (f<=1) return row.slice();
    const W = row.length;
    const out = new Array(W*f);
    for (let c=0;c<W-1;c++){
      const a = row[c], b = row[c+1];
      for (let t=0;t<f;t++){
        const alpha = t/f;
        out[c*f + t] = a + alpha*(b-a);
      }
    }
    out[out.length-1] = row[W-1];
    return out;
  }
  function upsample2D(M, fx, fy){
    if (fx<=1 && fy<=1) return M.map(r=>r.slice());
    const H = M.length, W = M[0].length;
    // cols
    const cols = new Array(H);
    for (let r=0;r<H;r++) cols[r] = upsampleRow(M[r], fx);
    // rows
    const H2 = H*fy, W2 = cols[0].length;
    const out = new Array(H2);
    for (let r=0;r<H-1;r++){
      const ra = cols[r], rb = cols[r+1];
      for (let t=0;t<fy;t++){
        const alpha = t/fy;
        const row = new Array(W2);
        for (let c=0;c<W2;c++) row[c] = ra[c] + alpha*(rb[c]-ra[c]);
        out[r*fy + t] = row;
      }
    }
    out[out.length-1] = cols[H-1].slice();
    return out;
  }
  function centersFromEdgesLogX(xe, ye, fx, fy){
    // x centers on log10 scale; y centers linear
    const xc = [], yc = [];
    for (let i=0;i<xe.length-1;i++){
      const lx0 = Math.log10(xe[i]), lx1 = Math.log10(xe[i+1]);
      for (let t=0;t<fx;t++){
        const lx = lx0 + (t+0.5)/fx*(lx1-lx0);
        xc.push(Math.pow(10, lx));
      }
    }
    for (let j=0;j<ye.length-1;j++){
      const y0 = ye[j], y1 = ye[j+1];
      for (let t=0;t<fy;t++){
        yc.push(y0 + (t+0.5)/fy*(y1-y0));
      }
    }
    return {xc:xc, yc:yc};
  }

  const Zplot = (UPSAMPLE>1) ? upsample2D(ZS, UPSAMPLE, UPSAMPLE) : ZS;
  const xy = (UPSAMPLE>1) ? centersFromEdgesLogX(xe, ye, UPSAMPLE, UPSAMPLE) : {xc:x, yc:y};

  // --- Percentile clipping (P99) ---
  function percentileFrom2D(M, q){
    const arr = [];
    for (let r=0;r<M.length;r++) for (let c=0;c<M[0].length;c++){ const v=M[r][c]; if (isFinite(v)) arr.push(v); }
    if (!arr.length) return 0;
    arr.sort((a,b)=>a-b);
    const idx = Math.floor((q/100) * (arr.length-1));
    return arr[idx];
  }
  const zmax99 = percentileFrom2D(Zplot, 99);

  const trace = {
    type:'heatmap',
    x:xy.xc, y:xy.yc, z:Zplot,
    colorscale:'Aggrnyl',
    zmin:0,
    zmax:(zmax99>0?zmax99:undefined)
  };

  const layout = { title:name, xaxis:{type:'log', title:'n_fragment'}, yaxis:{title:'TSSE'}, margin:{t:40}, dragmode:'select', shapes:[] };
  const config = {modeBarButtonsToAdd:[{name:'Download PNG (square)',title:'Download PNG (1024x1024)',icon:Plotly.Icons.camera,click:function(gd){var fn=(gd&&gd._fullLayout&&gd._fullLayout._meta&&gd._fullLayout._meta.filename)||(gd&&gd.layout&&gd.layout.title&&gd.layout.title.text)||'figure';Plotly.downloadImage(gd,{format:'png',width:1024,height:1024,filename:fn});}}],  responsive:true, displaylogo:false, toImageButtonOptions:{width:1024,height:1024,format:'svg',filename:'TSSE_grid_'+name} };
  Plotly.newPlot('DIV', [trace], layout, config);

  // Small info panel + selection updater
  const panel = document.createElement('div');
  panel.style.fontSize = '12px'; panel.style.marginTop = '8px';
  document.getElementById('DIV').parentElement.appendChild(panel);

  function summarizeSelection(range){
    if (!range || !range.x || !range.y) return;
    const x0 = Math.min(range.x[0], range.x[1]);
    const x1 = Math.max(range.x[0], range.x[1]);
    const y0 = Math.min(range.y[0], range.y[1]);
    const y1 = Math.max(range.y[0], range.y[1]);
    const rect = [{type:'rect', x0:x0, x1:x1, y0:y0, y1:y1, line:{width:2, dash:'dot'}, fillcolor:'rgba(0,0,0,0)'}];
    Plotly.relayout('DIV', {'shapes': rect});
    let cellCount = 0, fragSum = 0;
    for (let r=0;r<y.length;r++){
      for (let c=0;c<x.length;c++){
        const xc0 = xedges[c], xc1 = xedges[c+1];
        const yc0 = yedges[r], yc1 = yedges[r+1];
        if (xc1>=x0 && xc0<=x1 && yc1>=y0 && yc0<=y1){
          cellCount += (Z[r][c]||0);
          fragSum += (S[r][c]||0);
        }
      }
    }
    panel.innerHTML = 'Selected: <b>' + cellCount.toLocaleString() + '</b> cells' +
      (fragSum > 0 ? ' | Σfrags=<b>' + Math.round(fragSum).toLocaleString() + '</b>' : '');
  }

  const gd = document.getElementById('DIV');
  gd.on('plotly_selected', ev => {
    if (ev && ev.range) summarizeSelection(ev.range);
  });

  // Quadrant summary on click
  const findBin = (val, edges) => { let lo=0, hi=edges.length-1; while (lo<hi-1){ const mid=(lo+hi)>>1; if (val<edges[mid]) hi=mid; else lo=mid; } return Math.max(0, Math.min(edges.length-2, lo)); };

  gd.on('plotly_click', ev => {
    if (!ev || !ev.points || !ev.points.length) return;
    const xp = ev.points[0].x, yp = ev.points[0].y;
    const shapes = [
      {type:'line', x0:xp, x1:xp, y0:yedges[0], y1:yedges[yedges.length-1], line:{dash:'dot', width:2}},
      {type:'line', x0:xedges[0], x1:xedges[xedges.length-1], y0:yp, y1:yp, line:{dash:'dot', width:2}}
    ];
    Plotly.relayout('DIV', {'shapes': shapes});

    const cx = findBin(xp, xedges), cy = findBin(yp, yedges);
    const sums = [0,0,0,0], sumsF = [0,0,0,0];
    for (let r=0;r<y.length;r++){
      for (let c=0;c<x.length;c++){
        const z = (Z[r] && Z[r][c]) || 0;
        const ssum = (S[r] && S[r][c]) || 0;
        if (c>=cx && r>=cy){ sums[0]+=z; sumsF[0]+=ssum; }
        else if (c<cx && r>=cy){ sums[1]+=z; sumsF[1]+=ssum; }
        else if (c<cx && r<cy){ sums[2]+=z; sumsF[2]+=ssum; }
        else { sums[3]+=z; sumsF[3]+=ssum; }
      }
    }
    panel.innerHTML = 'Q1: ' + sums[0].toLocaleString() +
                      ' &nbsp; Q2: ' + sums[1].toLocaleString() +
                      '<br>Q3: ' + sums[2].toLocaleString() +
                      ' &nbsp; Q4: ' + sums[3].toLocaleString();
  });
})();""".replace('TS', vec_tsse).replace('NF', vec_nfrag).replace('NM', nm).replace('DIV', div_id))
        plots_html.append('<div class="card"><div class="h2">TSSE Grid (density)</div><div class="grid">' + ''.join(grid_divs) + '</div></div>')
        scripts.extend(grid_scripts)

    # 4) Fragment-count violin (n_fragment)
    if js_tv and js_tv.get('labels') and js_tv.get('nfrag'):
        plots_html.append('<div class="card"><div class="h2">Fragment Counts (per-sample) Violin (interactive)</div><div id="plot_frag_violin" class="plot"></div></div>')
        data_obj = json.dumps({'labels': js_tv['labels'], 'nfrag': js_tv['nfrag']})
        scripts.append(r'''(function(){
  const p = DATA;
  const traces = [];
  for (let i = 0; i < p.labels.length; i++) {
    const raw = p.nfrag[i] || [];
    const y = raw.filter(v => Number.isFinite(v) && v >= 0).map(v => Math.log10(v + 1));
    if (y.length) {
      traces.push({
        type: 'violin',
        y: y,
        name: p.labels[i],
        box: { visible: false },
        meanline: { visible: true },
        points: false
      });
    }
  }
  const layout = { yaxis: { title: 'log10(n_fragment + 1)' }, hovermode: 'closest' };
  const config = {
    modeBarButtonsToAdd: [{
      name: 'Download PNG (square)',
      title: 'Download PNG (1024x1024)',
      icon: Plotly.Icons.camera,
      click: function(gd){
        var fn = (gd && gd._fullLayout && gd._fullLayout._meta && gd._fullLayout._meta.filename)
                 || (gd && gd.layout && gd.layout.title && gd.layout.title.text)
                 || 'figure';
        Plotly.downloadImage(gd, { format: 'png', width: 1024, height: 1024, filename: fn });
      }
    }],
    responsive: true,
    displaylogo: false,
    toImageButtonOptions: { width: 1024, height: 1024, format: 'svg', filename: 'fragment_violin' }
  };
  Plotly.newPlot('plot_frag_violin', traces, layout, config);

  // --- Global mean-of-medians toggle (kept inside the same string) ---
  (function(){
    const gd = document.getElementById('plot_frag_violin');
    if (!gd) return;
    const lbl = document.createElement('label');
    lbl.style.fontSize = '12px'; lbl.style.display = 'inline-block'; lbl.style.marginTop = '6px';
    const cb = document.createElement('input'); cb.type = 'checkbox'; cb.id='fv_global_toggle'; cb.checked = true;
    lbl.appendChild(cb); lbl.appendChild(document.createTextNode(' Show global median'));
    if (gd.parentElement) gd.parentElement.appendChild(lbl);

    function _median(arr){
      const v = (arr||[]).filter(a => Number.isFinite(a)).slice().sort((a,b)=>a-b);
      const n = v.length; if (n===0) return null;
      return (n%2 ? v[(n-1)>>1] : 0.5*(v[n/2-1] + v[n/2]));
    }
    let meds = [];
    if (gd.data && gd.data.length) {
      for (let i=0;i<gd.data.length;i++) {
        const y = gd.data[i].y || [];
        const m = _median(y);
        if (m !== null) meds.push(m);
      }
    }
    const gmm = meds.length ? (meds.reduce((a,b)=>a+b,0)/meds.length) : null;

    function apply(show){
      let shapes = (gd.layout && gd.layout.shapes) ? gd.layout.shapes.slice() : [];
      shapes = shapes.filter(s => s && s._tag !== 'global_median_fv');
      if (show && Number.isFinite(gmm)) {
        shapes.push({type:'line', xref:'paper', yref:'y', x0:0, x1:1, y0:gmm, y1:gmm,
                     line:{color:'red', width:2, dash:'dash'}, _tag:'global_median_fv'});
      }
      Plotly.relayout(gd, {'shapes': shapes});
    }
    apply(cb.checked);
    cb.addEventListener('change', function(){ apply(cb.checked); });
  })();
})();'''.replace('DATA', data_obj))
    elif img_fviolin:
        plots_html.append('<div class="card"><div class="h2">Fragment Counts (per-sample) Violin</div><img src="' + img_fviolin + '" style="max-width:100%"></div>')


    # 5) Knee plot
    if js_tv and js_tv.get('nfrag'):
        plots_html.append('<div class="card"><div class="h2">Knee Plot (rank–count, log–log)</div><div id="plot_knee" class="plot"></div></div>')
        data_obj = json.dumps({'labels': js_tv['labels'], 'nfrag': js_tv['nfrag']})
        scripts.append("""(function(){
  const p = DATA;
  const traces = [];
  for (let i=0;i<p.labels.length;i++) {
    const vec = (p.nfrag[i]||[]).filter(v=>v>0).sort((a,b)=>b-a);
    const ranks = Array.from({length: vec.length}, (_,k)=>k+1);
    traces.push({x:ranks, y:vec, type:'scatter', mode:'lines', name:p.labels[i]});
  }
  const layout = {xaxis:{type:'log', title:'Rank (cells)'}, yaxis:{type:'log', title:'n_fragment'}, hovermode:'closest'};
  const config = {modeBarButtonsToAdd:[{name:'Download PNG (square)',title:'Download PNG (1024x1024)',icon:Plotly.Icons.camera,click:function(gd){var fn=(gd&&gd._fullLayout&&gd._fullLayout._meta&&gd._fullLayout._meta.filename)||(gd&&gd.layout&&gd.layout.title&&gd.layout.title.text)||'figure';Plotly.downloadImage(gd,{format:'png',width:1024,height:1024,filename:fn});}}], responsive:true, displaylogo:false, toImageButtonOptions:{width:1024,height:1024,format:'svg',filename:'knee_plot'}};
  Plotly.newPlot('plot_knee', traces, layout, config);
})();""".replace('DATA', data_obj))
    elif img_knee:
        plots_html.append('<div class="card"><div class="h2">Knee Plot</div><img src="' + img_knee + '" style="max-width:100%"></div>')

    html = []
    html.append('<!DOCTYPE html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">')
    html.append('<title>' + html_title + '</title>')
    html.append('<link rel="stylesheet" href="QC_reports/assets/css/report.css">')
    html.append(PLOTLY_LOADER)
    html.append('</head><body><div class="container">')
    html.append(header)
    html.append(summary_html)
    html.append(''.join(plots_html))
    html.append('</div><script>__qcReady(function(){')
    html.append(''.join(scripts))
    html.append('});</script></body></html>')

    out_html = out / 'qc_summary.html'
    out_html.write_text(''.join(html), encoding='utf-8')
    return str(out_html)

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Build interactive QC HTML report.")
    ap.add_argument("-o", "--out_dir", required=True, help="Output directory containing QC_reports/")
    ap.add_argument("--title", default="QC Summary", help="HTML title")
    ap.add_argument("--embed_assets", action="store_true", help="(reserved)")
    args = ap.parse_args()
    try:
        path = build_html_report(args.out_dir, html_title=args.title, embed_assets=args.embed_assets)
        print("[OK] HTML report generated at:", path)
    except Exception as e:
        print(f"[WARN] HTML report generation failed: {e}")
