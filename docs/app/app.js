/* app.js — the interactive layer.
 *
 * No framework and no build step: this is served straight from GitHub Pages, so everything
 * is vanilla ES modules and hand-drawn SVG. Charts are small enough that a plotting library
 * would cost more than it saves, and avoiding one keeps the page working with no network
 * beyond this origin.
 */

import { solveOperatingPoint, DEFAULT_PARAMS } from "./fvcb.js";

/* Chart colours come from the stylesheet, not from constants here, so the SVGs follow the
 * viewer's light/dark theme along with everything else. Hard-coding them painted ink-black
 * axis lines onto a dark background. Re-read whenever the scheme changes. */
const C = {};
function readTheme() {
  const cs = getComputedStyle(document.documentElement);
  const v = (name, fallback) => (cs.getPropertyValue(name).trim() || fallback);
  Object.assign(C, {
    ink: v("--ink", "#1a1d21"),
    soft: v("--ink-soft", "#4a5058"),
    faint: v("--ink-faint", "#8a8f98"),
    line: v("--line", "#dfe2e6"),
    accent: v("--accent", "#0072b2"),
    vermilion: v("--vermilion", "#d55e00"),
    orange: v("--orange", "#e69f00"),
    green: v("--green", "#009e73"),
    sky: v("--sky", "#56b4e9"),
    purple: v("--purple", "#cc79a7"),
    grey: v("--ink-faint", "#8a8f98"),
    // Marker outlines need the page background, not white, or they ring dark-mode dots.
    bg: v("--bg", "#ffffff"),
  });
}
readTheme();
const SVG = "http://www.w3.org/2000/svg";

const el = (tag, attrs = {}, text) => {
  const n = document.createElementNS(SVG, tag);
  for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, v);
  if (text !== undefined) n.textContent = text;
  return n;
};
const clear = (node) => { while (node.firstChild) node.removeChild(node.firstChild); };
const json = (name) => fetch(`data/${name}`).then((r) => r.json());

/* ───────────────────────────────────────────────── a minimal chart helper */

function axes(svg, { w, h, pad, xlim, ylim, xlabel, ylabel, yticks = 5, xticks = 5 }) {
  const [x0, x1] = xlim, [y0, y1] = ylim;
  const X = (v) => pad.l + ((v - x0) / (x1 - x0)) * (w - pad.l - pad.r);
  const Y = (v) => h - pad.b - ((v - y0) / (y1 - y0)) * (h - pad.t - pad.b);

  for (let i = 0; i <= yticks; i++) {
    const v = y0 + ((y1 - y0) * i) / yticks;
    svg.append(el("line", { x1: pad.l, x2: w - pad.r, y1: Y(v), y2: Y(v),
                            stroke: C.line, "stroke-dasharray": "2 3" }));
    svg.append(el("text", { x: pad.l - 7, y: Y(v) + 3.5, "text-anchor": "end",
                            "font-size": 10, fill: C.faint }, fmt(v)));
  }
  for (let i = 0; i <= xticks; i++) {
    const v = x0 + ((x1 - x0) * i) / xticks;
    svg.append(el("text", { x: X(v), y: h - pad.b + 15, "text-anchor": "middle",
                            "font-size": 10, fill: C.faint }, fmt(v)));
  }
  svg.append(el("line", { x1: pad.l, x2: w - pad.r, y1: h - pad.b, y2: h - pad.b,
                          stroke: C.faint }));
  svg.append(el("line", { x1: pad.l, x2: pad.l, y1: pad.t, y2: h - pad.b, stroke: C.faint }));
  if (xlabel) svg.append(el("text", { x: (pad.l + w - pad.r) / 2, y: h - 4,
                                      "text-anchor": "middle", "font-size": 11,
                                      fill: C.soft }, xlabel));
  if (ylabel) {
    const t = el("text", { x: 12, y: (pad.t + h - pad.b) / 2, "text-anchor": "middle",
                           "font-size": 11, fill: C.soft,
                           transform: `rotate(-90 12 ${(pad.t + h - pad.b) / 2})` }, ylabel);
    svg.append(t);
  }
  return { X, Y };
}

const fmt = (v) => {
  const a = Math.abs(v);
  if (a >= 1000) return v.toFixed(0);
  if (a >= 10) return v.toFixed(0);
  if (a >= 1) return v.toFixed(1);
  return v.toFixed(2);
};

const path = (pts, { X, Y }) =>
  pts.map(([x, y], i) => `${i ? "L" : "M"}${X(x).toFixed(1)},${Y(y).toFixed(1)}`).join(" ");

/* ─────────────────────────────────────────────────────────── the model tab */

const SLIDERS = [
  { k: "Ca", label: "canister CO₂", unit: "µmol mol⁻¹", min: 60, max: 420, step: 5, val: 100 },
  { k: "g_bl", label: "boundary layer g_bl", unit: "mol m⁻² s⁻¹",
    min: 0.05, max: 1.2, step: 0.001, val: 0.291 },
  { k: "Q", label: "light (PPFD)", unit: "µmol m⁻² s⁻¹", min: 0, max: 2000, step: 25, val: 1000 },
  { k: "TleafC", label: "leaf temperature", unit: "°C", min: 10, max: 40, step: 0.5, val: 25 },
  { k: "Vcmax25", label: "Vcmax at 25 °C", unit: "µmol m⁻² s⁻¹", min: 20, max: 180, step: 1,
    val: 90 },
  { k: "Jmax25", label: "Jmax at 25 °C", unit: "µmol m⁻² s⁻¹", min: 40, max: 300, step: 1,
    val: 153 },
  { k: "g_s", label: "stomatal g_s", unit: "mol m⁻² s⁻¹", min: 0.02, max: 0.6, step: 0.005,
    val: 0.20 },
  { k: "g_m", label: "mesophyll g_m", unit: "mol m⁻² s⁻¹", min: 0.05, max: 0.8, step: 0.005,
    val: 0.30 },
];

const PRESETS = {
  vented:          { Ca: 400, g_bl: 0.540, Q: 1000 },
  "sealed-ground": { Ca: 100, g_bl: 0.540, Q: 1000 },
  "sealed-flight": { Ca: 100, g_bl: 0.291, Q: 1000 },
  "canopy-ug":     { Ca: 100, g_bl: 0.109, Q: 1000 },
};

const state = Object.fromEntries(SLIDERS.map((s) => [s.k, s.val]));

function params() {
  return { ...DEFAULT_PARAMS, Vcmax25: state.Vcmax25, Jmax25: state.Jmax25,
           g_s: state.g_s, g_m: state.g_m };
}
const solveAt = (gBl, Ca) =>
  solveOperatingPoint(gBl, params(), { Ca, TleafC: state.TleafC, Q: state.Q });

function buildControls() {
  const host = document.getElementById("controls");
  for (const s of SLIDERS) {
    const wrap = document.createElement("label");
    wrap.innerHTML =
      `<span class="row"><span class="name">${s.label}</span>` +
      `<span class="val" id="v_${s.k}"></span></span>` +
      `<input type="range" id="s_${s.k}" min="${s.min}" max="${s.max}" step="${s.step}">`;
    host.append(wrap);
    const input = wrap.querySelector("input");
    input.value = state[s.k];
    input.addEventListener("input", () => {
      state[s.k] = parseFloat(input.value);
      renderModel();
    });
  }
}

function renderModel() {
  for (const s of SLIDERS) {
    document.getElementById(`v_${s.k}`).textContent =
      `${state[s.k] >= 100 ? state[s.k].toFixed(0) : state[s.k].toFixed(3).replace(/0+$/, "").replace(/\.$/, "")} ${s.unit}`;
    document.getElementById(`s_${s.k}`).value = state[s.k];
  }

  const r = solveAt(state.g_bl, state.Ca);
  const cards = [
    ["Cc", r.Cc.toFixed(1), "µmol mol⁻¹"],
    ["oxygenation φ", (r.phi * 100).toFixed(1), "%"],
    ["Vo / Vc", r.Vo_over_Vc.toFixed(3), ""],
    ["assimilation A", r.A.toFixed(2), "µmol m⁻² s⁻¹"],
    ["photoresp. Rp", r.Rp.toFixed(2), "µmol m⁻² s⁻¹"],
    ["limited by", r.limiting, ""],
  ];
  const host = document.getElementById("readout");
  host.innerHTML = cards.map(([k, v, u]) =>
    `<div class="card"><div class="k">${k}</div><div class="v">${v}</div>` +
    `<div class="u">${u}</div></div>`).join("");

  drawCurve(r);
  drawGap();
}

const CA_GRID = Array.from({ length: 80 }, (_, i) => 60 + (i * (420 - 60)) / 79);

function drawCurve(op) {
  const svg = document.getElementById("curve");
  clear(svg);
  const w = 560, h = 320, pad = { l: 52, r: 52, t: 14, b: 40 };
  const phi = CA_GRID.map((ca) => [ca, solveAt(state.g_bl, ca).phi * 100]);
  const A = CA_GRID.map((ca) => [ca, solveAt(state.g_bl, ca).A]);
  const amax = Math.max(4, ...A.map((d) => d[1]));

  const ax = axes(svg, { w, h, pad, xlim: [60, 420], ylim: [0, 70],
                         xlabel: "canister CO₂ (µmol mol⁻¹)", ylabel: "oxygenation φ (%)" });
  const AY = (v) => h - pad.b - (v / amax) * (h - pad.t - pad.b);

  svg.append(el("path", { d: path(A.map(([x, y]) => [x, y]), { X: ax.X, Y: AY }),
                          fill: "none", stroke: C.green, "stroke-width": 1.6,
                          "stroke-dasharray": "5 3" }));
  svg.append(el("path", { d: path(phi, ax), fill: "none", stroke: C.vermilion,
                          "stroke-width": 2.2 }));

  svg.append(el("circle", { cx: ax.X(state.Ca), cy: ax.Y(op.phi * 100), r: 5.5,
                            fill: C.vermilion, stroke: C.bg, "stroke-width": 1.6 }));
  svg.append(el("circle", { cx: ax.X(state.Ca), cy: AY(op.A), r: 4.5,
                            fill: C.green, stroke: C.bg, "stroke-width": 1.4 }));

  for (let i = 0; i <= 4; i++) {
    const v = (amax * i) / 4;
    svg.append(el("text", { x: w - pad.r + 7, y: AY(v) + 3.5, "font-size": 10,
                            fill: C.green }, fmt(v)));
  }
  svg.append(el("text", { x: w - 8, y: (pad.t + h - pad.b) / 2, "text-anchor": "middle",
                          "font-size": 11, fill: C.green,
                          transform: `rotate(90 ${w - 8} ${(pad.t + h - pad.b) / 2})` },
                         "assimilation A (µmol m⁻² s⁻¹)"));
  svg.append(el("text", { x: pad.l + 8, y: pad.t + 12, "font-size": 10.5,
                          fill: C.vermilion, "font-weight": 600 }, "φ oxygenation"));
  svg.append(el("text", { x: pad.l + 8, y: pad.t + 26, "font-size": 10.5,
                          fill: C.green, "font-weight": 600 }, "A assimilation (dashed)"));
}

function drawGap() {
  const svg = document.getElementById("gap");
  clear(svg);
  const w = 560, h = 210, pad = { l: 52, r: 16, t: 14, b: 40 };
  const gaps = CA_GRID.map((ca) =>
    [ca, (solveAt(0.291, ca).phi - solveAt(0.540, ca).phi) * 100]);
  const ymax = Math.max(0.6, ...gaps.map((d) => d[1])) * 1.25;
  const ax = axes(svg, { w, h, pad, xlim: [60, 420], ylim: [0, ymax], yticks: 4,
                         xlabel: "canister CO₂ (µmol mol⁻¹)",
                         ylabel: "φ(µg) − φ(1 g)  (pp)" });
  svg.append(el("path", { d: path(gaps, ax), fill: "none", stroke: C.accent,
                          "stroke-width": 2.2 }));
  const here = (solveAt(0.291, state.Ca).phi - solveAt(0.540, state.Ca).phi) * 100;
  svg.append(el("circle", { cx: ax.X(state.Ca), cy: ax.Y(here), r: 5,
                            fill: C.accent, stroke: C.bg, "stroke-width": 1.5 }));
  svg.append(el("text", { x: ax.X(state.Ca) + 9, y: ax.Y(here) - 7, "font-size": 10.5,
                          fill: C.accent, "font-weight": 600 },
                         `${here.toFixed(2)} pp at your CO₂`));
}

/* ──────────────────────────────────────────────────────── enclosures tab */

async function buildEnclosure() {
  const ts = await json("enclosure_timeseries.json");
  const svg = document.getElementById("encl");
  const w = 620, h = 320, pad = { l: 62, r: 130, t: 14, b: 42 };
  const all = Object.values(ts).flatMap((d) => d.co2).filter((v) => v !== null);
  const ax = axes(svg, { w, h, pad, xlim: [0, 40], ylim: [Math.min(...all) * 1.1,
                         Math.max(...all) * 1.15],
                         xlabel: "model step (×10³)",
                         ylabel: "enclosure CO₂ excess (model units)" });
  const colour = { "BRIC light": C.vermilion, "BRIC dark": "#8c2e0f",
                   "CARA tape": C.orange, "VEGGIE vented": C.accent, "open (ref)": C.grey };
  let i = 0;
  for (const [name, d] of Object.entries(ts)) {
    // Vented cases have no enclosure mean to report and arrive as null; they are listed in
    // the legend as flat-at-zero rather than drawn as a broken line.
    const pts = d.step.map((s, k) => [s / 1000, d.co2[k]]).filter((p) => p[1] !== null);
    if (!pts.length) {
      svg.append(el("text", { x: w - pad.r + 10, y: pad.t + 16 + i * 16, "font-size": 10.5,
                              fill: colour[name] || C.grey },
                             `${name} (vents, no drift)`));
      i++;
      continue;
    }
    svg.append(el("path", { d: path(pts, ax), fill: "none",
                            stroke: colour[name] || C.grey, "stroke-width": 1.8,
                            "stroke-dasharray": name.includes("open") ? "5 3" : "" }));
    svg.append(el("text", { x: w - pad.r + 10, y: pad.t + 16 + i * 16, "font-size": 10.5,
                            fill: colour[name] || C.grey }, name));
    i++;
  }

  const ret = await json("carbon_retention.json");
  document.getElementById("retention").innerHTML = ret.map((r) =>
    `<div class="card"><div class="k">${r.enclosure}</div>` +
    `<div class="v">${r["12h carbon (% Earth)"]}%</div>` +
    `<div class="u">of Earth's 12 h carbon gain</div></div>`).join("");
}

/* ───────────────────────────────────────────────────────────── omics tab */

let volcanoes = {}, geneSets = {}, layer = "transcriptome", highlight = "", query = "";

async function buildOmics() {
  volcanoes = {
    transcriptome: await json("volcano_transcriptome.json"),
    proteome: await json("volcano_proteome.json"),
  };
  geneSets = await json("gene_sets.json");

  const pick = document.getElementById("setpick");
  pick.innerHTML = `<option value="">none</option>` +
    Object.keys(geneSets).map((k) => `<option value="${k}">${k}</option>`).join("");
  pick.addEventListener("change", () => { highlight = pick.value; drawVolcano(); });

  document.getElementById("genesearch").addEventListener("input", (e) => {
    query = e.target.value.trim().toUpperCase();
    drawVolcano();
  });
  document.querySelectorAll("[data-layer]").forEach((b) =>
    b.addEventListener("click", () => {
      document.querySelectorAll("[data-layer]").forEach((x) => x.classList.remove("on"));
      b.classList.add("on");
      layer = b.dataset.layer;
      drawVolcano();
    }));
  drawVolcano();
}

function drawVolcano() {
  const v = volcanoes[layer];
  const svg = document.getElementById("volcano");
  clear(svg);
  const w = 620, h = 400, pad = { l: 56, r: 18, t: 14, b: 42 };
  const xs = v.points.map((p) => p[1]), ys = v.points.map((p) => p[2]);
  const xm = Math.max(...xs.map(Math.abs)) * 1.05;
  const ax = axes(svg, { w, h, pad, xlim: [-xm, xm], ylim: [0, Math.max(...ys) * 1.06],
                         xlabel: "log₂ fold change (flight / ground)",
                         ylabel: layer === "proteome" ? "−log₁₀ adjusted p" : "−log₁₀ FDR" });

  // Highlighted set is drawn last so it is never buried; the gene-set index is keyed on
  // AGI, so it only applies to the transcriptome layer.
  const set = highlight && layer === "transcriptome"
    ? new Set(geneSets[highlight]) : new Set();
  const hits = [];
  const frag = document.createDocumentFragment();
  for (const [id, x, y, sig] of v.points) {
    const inSet = set.has(id);
    const isHit = query && id.toUpperCase().includes(query);
    if (inSet || isHit) { hits.push([id, x, y, inSet, isHit]); continue; }
    frag.append(el("circle", { cx: ax.X(x), cy: ax.Y(y), r: sig ? 1.9 : 1.3,
                               fill: sig ? (x > 0 ? C.vermilion : C.accent) : C.line,
                               "fill-opacity": sig ? 0.75 : 0.5 }));
  }
  svg.append(frag);
  for (const [id, x, y, inSet, isHit] of hits) {
    svg.append(el("circle", { cx: ax.X(x), cy: ax.Y(y), r: isHit ? 5 : 3,
                              fill: isHit ? C.green : C.purple,
                              stroke: C.bg, "stroke-width": isHit ? 1.4 : 0.7 }));
    if (isHit) {
      svg.append(el("text", { x: ax.X(x) + 8, y: ax.Y(y) - 6, "font-size": 10.5,
                              fill: C.green, "font-weight": 600 }, id));
    }
  }

  const thresh = -Math.log10(v.alpha);
  svg.append(el("line", { x1: pad.l, x2: w - pad.r, y1: ax.Y(thresh), y2: ax.Y(thresh),
                          stroke: C.ink, "stroke-dasharray": "4 3", "stroke-width": 0.8 }));

  const matched = query ? hits.filter((hh) => hh[4]).length : 0;
  document.getElementById("volcap").textContent =
    `${v.n_sig.toLocaleString()} of ${v.n_total.toLocaleString()} features significant at ` +
    `${v.alpha}. ${v.n_plotted.toLocaleString()} plotted.` +
    (highlight && layer === "transcriptome"
      ? `  Highlighted: ${highlight} (${set.size} measured).` : "") +
    (query ? `  Search "${query}": ${matched} match${matched === 1 ? "" : "es"}.` : "");
}

/* ──────────────────────────────────────────────────────── tables (shared) */

function table(node, rows, cols, opts = {}) {
  clear(node);
  const thead = document.createElement("thead");
  const tr = document.createElement("tr");
  for (const c of cols) {
    const th = document.createElement("th");
    th.textContent = c.label;
    th.addEventListener("click", () => {
      const dir = node._sort === c.key ? -(node._dir || 1) : 1;
      node._sort = c.key; node._dir = dir;
      rows.sort((a, b) => {
        const x = a[c.key], y = b[c.key];
        const nx = parseFloat(x), ny = parseFloat(y);
        if (!isNaN(nx) && !isNaN(ny)) return (nx - ny) * dir;
        return String(x).localeCompare(String(y)) * dir;
      });
      table(node, rows, cols, opts);
    });
    tr.append(th);
  }
  thead.append(tr); node.append(thead);

  const tbody = document.createElement("tbody");
  for (const r of rows) {
    const row = document.createElement("tr");
    if (opts.rowClass) row.className = opts.rowClass(r);
    for (const c of cols) {
      const td = document.createElement("td");
      if (c.html) td.innerHTML = c.html(r); else td.textContent = c.get ? c.get(r) : r[c.key];
      if (c.num) td.className = "num";
      row.append(td);
    }
    tbody.append(row);
  }
  node.append(tbody);
}

const sci = (v) => {
  const n = parseFloat(v);
  if (isNaN(n)) return String(v);
  return n < 0.001 ? n.toExponential(2) : n.toFixed(4);
};

/* ──────────────────────────────────────────────────────── pathways tab */

async function buildPathways() {
  const sig = await json("enrichment.json");
  const carbon = await json("carbon_pathways.json");
  const seen = new Set(sig.map((r) => r.pathway));
  const rows = sig.concat(carbon.filter((r) => !seen.has(r.pathway)));

  const cols = [
    { key: "db", label: "DB" },
    { key: "pathway", label: "Pathway" },
    { key: "p_gene", label: "genes", num: true, get: (r) => sci(r.p_gene) },
    { key: "p_protein", label: "proteins", num: true, get: (r) => sci(r.p_protein) },
    { key: "p_metabolite", label: "metabolites", num: true, get: (r) => sci(r.p_metabolite) },
    { key: "p_combined_fisher", label: "combined", num: true,
      html: (r) => {
        const v = parseFloat(r.p_combined_fisher);
        return `<span class="tag ${v < 0.05 ? "sig" : "no"}">${sci(v)}</span>`;
      } },
  ];
  const node = document.getElementById("pathtable");
  const render = (q) => table(node, rows.filter((r) =>
    !q || r.pathway.toLowerCase().includes(q)), cols);
  render("");
  document.getElementById("pathfilter").addEventListener("input", (e) =>
    render(e.target.value.trim().toLowerCase()));

  const hubs = await json("metabolite_hubs.json");
  table(document.getElementById("hubtable"), hubs, [
    { key: "compound", label: "Compound (PREDICTED value)" },
    { key: "kegg_id", label: "KEGG" },
    { key: "block", label: "Block" },
    { key: "DE_neighbours", label: "measured DE genes nearby", num: true },
    { key: "FDR", label: "FDR", num: true,
      html: (r) => `<span class="tag ${parseFloat(r.FDR) < 0.05 ? "sig" : "no"}">` +
                   `${r.FDR}</span>` },
  ]);
}

/* ────────────────────────────────────────────────────────── ladder tab */

async function buildLadder() {
  const con = (await json("ladder_contrasts.json"))
    .filter((r) => r.contrast === "starvation vs photosynthesis" && r.separation !== null)
    .sort((a, b) => a.separation - b.separation);

  const svg = document.getElementById("ladderfig");
  const w = 620, h = 300, pad = { l: 108, r: 130, t: 14, b: 44 };
  const lim = Math.max(...con.map((r) => Math.abs(r.separation))) * 1.15;
  const ax = axes(svg, { w, h, pad, xlim: [-lim, lim], ylim: [-0.6, con.length - 0.4],
                         yticks: 1, xlabel: "starvation − photosynthesis (log₂FC difference)" });
  svg.append(el("line", { x1: ax.X(0), x2: ax.X(0), y1: pad.t, y2: h - pad.b,
                          stroke: C.ink, "stroke-width": 1 }));
  con.forEach((r, i) => {
    const y = ax.Y(i), bw = ax.X(r.separation) - ax.X(0);
    const col = r.light === "light" ? C.orange : C.accent;
    svg.append(el("rect", { x: Math.min(ax.X(0), ax.X(r.separation)), y: y - 9,
                            width: Math.abs(bw), height: 18, fill: col }));
    svg.append(el("text", { x: pad.l - 9, y: y + 4, "text-anchor": "end", "font-size": 10.5,
                            fill: C.ink }, r.accession.replace("OSD-", "")));
    svg.append(el("text", { x: w - pad.r + 10, y: y + 4, "font-size": 10,
                            fill: C.faint }, `${r.enclosure_class} · ${r.light}`));
  });
  svg.append(el("text", { x: pad.l, y: pad.t + 4, "font-size": 10.5, fill: C.soft },
                         "every lit study sits above every dark study"));

  const studies = await json("studies.json");
  table(document.getElementById("studytable"), studies, [
    { key: "accession", label: "Study" },
    { key: "hardware", label: "Hardware" },
    { key: "enclosure_class", label: "Class" },
    { key: "tissue", label: "Tissue" },
    { key: "light", label: "Light" },
    { key: "include", label: "Used?",
      html: (r) => `<span class="tag ${r.include === "yes" ? "sig" : "no"}">` +
                   `${r.include === "yes" ? "yes" : "excluded"}</span>` },
    { key: "reason", label: "Why", get: (r) => r.reason },
  ], { rowClass: (r) => (r.include === "yes" ? "" : "excluded") });
}

/* ────────────────────────────────────────────────────────────── wiring */

document.querySelectorAll("#tabs button").forEach((b) =>
  b.addEventListener("click", () => {
    document.querySelectorAll("#tabs button").forEach((x) => x.classList.remove("on"));
    document.querySelectorAll(".panel").forEach((p) => p.classList.remove("on"));
    b.classList.add("on");
    document.getElementById(b.dataset.tab).classList.add("on");
  }));

document.querySelectorAll("[data-preset]").forEach((b) =>
  b.addEventListener("click", () => {
    Object.assign(state, PRESETS[b.dataset.preset]);
    document.querySelectorAll("[data-preset]").forEach((x) => x.classList.remove("on"));
    b.classList.add("on");
    renderModel();
  }));

/* Follow the OS scheme live: re-read the tokens and repaint every SVG. */
if (window.matchMedia) {
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => {
    readTheme();
    renderModel();
    buildEnclosure().catch(() => {});
    drawVolcano();
    buildLadder().catch(() => {});
  });
}

buildControls();
renderModel();
buildEnclosure().catch((e) => console.error("enclosure", e));
buildOmics().catch((e) => console.error("omics", e));
buildPathways().catch((e) => console.error("pathways", e));
buildLadder().catch((e) => console.error("ladder", e));
