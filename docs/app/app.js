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

/* Human labels for the gene sets. Mirrors PRETTY in scripts/figures/style.py so the site
 * and the manuscript figures call the same set by the same name. */
const PRETTY = {
  ath00630: "Glyoxylate & dicarboxylate",
  ath00710: "Carbon fixation",
  ath00500: "Starch & sucrose",
  ath00010: "Glycolysis",
  photorespiration_core: "Photorespiration C2 enzymes",
  carbon_starvation_DIN: "Carbon starvation (DIN)",
  fermentation: "Fermentation",
  hypoxia_responsive: "Hypoxia-responsive",
  photosynthesis_apparatus: "Photosystem & light harvesting",
  rubisco: "Rubisco",
};

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

let enclosureData = null, enclosureWhere = "surface";

async function buildEnclosure() {
  enclosureData = await json("enclosure_timeseries.json");
  document.querySelectorAll("[data-encl]").forEach((chip) => {
    chip.classList.toggle("on", chip.dataset.encl === enclosureWhere);
    chip.addEventListener("click", () => {
      enclosureWhere = chip.dataset.encl;
      document.querySelectorAll("[data-encl]").forEach(
        (c) => c.classList.toggle("on", c.dataset.encl === enclosureWhere));
      drawEnclosure();
    });
  });
  drawEnclosure();

  const ret = await json("carbon_retention.json");
  document.getElementById("retention").innerHTML = ret.map((r) =>
    `<div class="card"><div class="k">${r.enclosure}</div>` +
    `<div class="v">${r["12h carbon (% Earth)"]}%</div>` +
    `<div class="u">of Earth's 12 h carbon gain</div></div>`).join("");
}

const ENCL_COLOUR = { "BRIC light": () => C.vermilion, "BRIC dark": () => "#8c2e0f",
                      "CARA light": () => C.orange, "CARA dark": () => "#b05a00", "VEGGIE vented": () => C.accent,
                      "open (ref)": () => C.grey };

function drawEnclosure() {
  const svg = document.getElementById("encl");
  if (!svg || !enclosureData) return;
  clear(svg);
  const ts = enclosureData, key = enclosureWhere;
  const w = 620, h = 320, pad = { l: 62, r: 130, t: 14, b: 42 };
  const all = Object.values(ts).flatMap((d) => d[key]).filter((v) => v !== null);
  const ax = axes(svg, { w, h, pad, xlim: [0, 40],
                         ylim: [Math.min(...all) * 1.1, Math.max(...all) * 1.15],
                         xlabel: "model step (×10³)",
                         ylabel: key === "surface" ? "leaf-surface CO₂ excess (model units)"
                                                   : "enclosure-mean CO₂ excess (model units)" });
  let i = 0;
  for (const [name, d] of Object.entries(ts)) {
    const pts = d.step.map((s, k) => [s / 1000, d[key][k]]).filter((p) => p[1] !== null);
    const colour = (ENCL_COLOUR[name] || (() => C.grey))();
    // A vented case has no closed volume, so it has no enclosure mean — not "no drift",
    // which is what the legend used to say. Say which series is missing, and why.
    const label = pts.length ? name : `${name} — vents, no enclosure volume`;
    if (pts.length) {
      svg.append(el("path", { d: path(pts, ax), fill: "none", stroke: colour,
                              "stroke-width": 1.8,
                              "stroke-dasharray": name.includes("open") ? "5 3" : "" }));
    }
    svg.append(el("text", { x: w - pad.r + 10, y: pad.t + 16 + i * 16, "font-size": 10.5,
                            fill: pts.length ? colour : C.faint }, label));
    i++;
  }

  const cap = document.getElementById("enclcap");
  if (cap) {
    const final = (n) => {
      const v = ts[n][key].filter((x) => x !== null);
      return v.length ? v[v.length - 1] : null;
    };
    cap.textContent = key === "surface"
      ? `CO₂ excess at the leaf surface, the quantity the photosynthesis model consumes. `
        + `All five cases are defined here. By the end of the run the lit sealed canister `
        + `is lowest (${final("BRIC light").toFixed(3)}), the open reference and taped `
        + `canister are close behind (${final("open (ref)").toFixed(3)} and `
        + `${final("CARA light").toFixed(3)}), VEGGIE is the shallowest drawdown `
        + `(${final("VEGGIE vented").toFixed(3)}), and the dark canister rises `
        + `(${final("BRIC dark").toFixed(3)} and ${final("CARA dark").toFixed(3)}) because respiration has no uptake to offset it.`
      : `CO₂ excess averaged over the enclosure volume. Only the three closed cases have `
        + `one: a vented case has no closed volume to average, so VEGGIE and the open `
        + `reference are absent here rather than flat. Switch to the leaf surface to see `
        + `all five.`;
  }
}

/* ───────────────────────────────────────────────────────────── omics tab */

let volcanoes = {}, geneSets = {}, layer = "transcriptome", highlight = "", query = "";

async function buildOmics() {
  volcanoes = {
    transcriptome: await json("volcano_transcriptome.json"),
    proteome: await json("volcano_proteome.json"),
  };
  // Two membership indexes share one namespace and one control: the curated KEGG-derived
  // sets, and the pathways PaintOmics called significant. The Sankey below highlights
  // through the same control, which is why they are merged here rather than kept apart.
  const curated = await json("gene_sets.json");
  const sigPaths = await json("pathway_membership.json");
  geneSets = { ...curated, ...sigPaths };

  const group = (label, keys) => keys.length
    ? `<optgroup label="${label}">` +
      keys.map((k) => `<option value="${k}">${k}</option>`).join("") + `</optgroup>`
    : "";
  const pick = document.getElementById("setpick");
  pick.innerHTML = `<option value="">none</option>` +
    group("Curated gene sets", Object.keys(curated)) +
    group("Significant pathways", Object.keys(sigPaths));
  pick.addEventListener("change", () => { highlight = pick.value; drawVolcano(); });
  // The heatmap below also drives this, so keep the two in step.
  window.__setHighlight = (name) => {
    highlight = name;
    pick.value = name;
    drawVolcano();
    document.getElementById("volcano").scrollIntoView({ block: "center",
                                                       behavior: "smooth" });
  };

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

  // Highlighted members are drawn last so they are never buried. The index carries a
  // separate id list per layer — AGI for transcripts, UniProt for proteins — so the
  // control works on both. It used to be gated to the transcriptome and did nothing at
  // all on the proteome tab.
  const entry = highlight ? geneSets[highlight] : null;
  const set = entry ? new Set(entry[layer]) : new Set();
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
  let note = "";
  if (entry) {
    const measured = entry.measured[layer];
    const gap = measured - set.size;
    // Say what is NOT shown as well as what is. A member can be measured but carry no
    // adjusted p after DESeq2 independent filtering, so it has no y-value to plot.
    note = measured === 0
      ? `  ${highlight}: none of its ${entry.total} members was measured on this layer.`
      : `  Highlighted ${highlight}: ${set.size} of ${measured} measured on this layer` +
        (gap ? ` (${gap} measured but with no adjusted p, so not plottable)` : "") + ".";
  }
  document.getElementById("volcap").textContent =
    `${v.n_sig.toLocaleString()} of ${v.n_total.toLocaleString()} features significant at ` +
    `${v.alpha}. ${v.n_plotted.toLocaleString()} plotted.` + note +
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

let ladderContrasts = null, ladderPick = "starvation vs photosynthesis";

async function buildLadder() {
  // Three set-vs-set contrasts ship, and the chart used to draw one of them and discard
  // the other twelve rows — including both contrasts that test the photorespiration
  // claim, which is the paper's central one. All three are selectable now.
  ladderContrasts = await json("ladder_contrasts.json");
  const names = [...new Set(ladderContrasts.map((r) => r.contrast))];
  const row = document.getElementById("ladderpick");
  if (row) {
    row.innerHTML = names.map((n) =>
      `<button class="chip${n === ladderPick ? " on" : ""}" data-contrast="${n}">` +
      `${n}</button>`).join("");
    row.querySelectorAll("[data-contrast]").forEach((chip) => {
      chip.addEventListener("click", () => {
        ladderPick = chip.dataset.contrast;
        row.querySelectorAll("[data-contrast]").forEach(
          (c) => c.classList.toggle("on", c.dataset.contrast === ladderPick));
        drawLadderFig();
      });
    });
  }
  drawLadderFig();
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

function drawLadderFig() {
  const svg = document.getElementById("ladderfig");
  if (!svg || !ladderContrasts) return;
  clear(svg);
  const con = ladderContrasts
    .filter((r) => r.contrast === ladderPick && r.separation !== null)
    .sort((a, b) => a.separation - b.separation);
  if (!con.length) return;

  const w = 620, h = 300, pad = { l: 108, r: 130, t: 14, b: 44 };
  const lim = Math.max(...con.map((r) => Math.abs(r.separation))) * 1.15;
  const ax = axes(svg, { w, h, pad, xlim: [-lim, lim], ylim: [-0.6, con.length - 0.4],
                         yticks: 1, xlabel: `${ladderPick} (log₂FC difference)` });
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

  // Whether illumination separates the studies is TRUE for two of the three contrasts and
  // false for the third, so it is read off the data rather than asserted.
  const lit = con.filter((r) => r.light === "light").map((r) => r.separation);
  const dark = con.filter((r) => r.light === "dark").map((r) => r.separation);
  const clean = lit.length && dark.length && Math.min(...lit) > Math.max(...dark);
  svg.append(el("text", { x: pad.l, y: pad.t + 4, "font-size": 10.5,
                          fill: clean ? C.soft : C.vermilion },
                         clean ? "every lit study sits above every dark study"
                               : "lit and dark overlap on this contrast"));

  const cap = document.getElementById("laddercap");
  if (cap) {
    const f = (v) => v.toFixed(2);
    // .sort() on the formatted strings orders "-0.05" before "-1.13"; sort the numbers.
    const list = (v) => [...v].sort((a, b) => a - b).map(f).join(", ");
    cap.textContent = clean
      ? `Every illuminated study sits above every dark one on this contrast `
        + `(lit ${list(lit)}; dark ${list(dark)}). `
        + `OSD-678 light against OSD-678 dark is the same hardware, genotype and flight — `
        + `an internal control with no confound.`
      : `Illumination does NOT separate the studies on this contrast: lit spans `
        + `${f(Math.min(...lit))} to ${f(Math.max(...lit))} and dark `
        + `${f(Math.min(...dark))} to ${f(Math.max(...dark))}, so they overlap. `
        + `The separation the paper reports is on the other two contrasts.`;
  }
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
    drawHeatmap();
    drawSankey();
    buildLadder().catch(() => {});
  });
}

buildControls();
renderModel();
buildEnclosure().catch((e) => console.error("enclosure", e));
buildOmics().catch((e) => console.error("omics", e));
buildPathways().catch((e) => console.error("pathways", e));
buildLadder().catch((e) => console.error("ladder", e));
buildHeatmap().catch((e) => console.error("heatmap", e));
buildSankey().catch((e) => console.error("sankey", e));

/* ────────────────────────────────────────────────── heatmap: sets x studies */

let ladderRows = [];

/* Diverging blue -> white -> red, symmetric about zero, in the page's own tokens so it
 * follows light and dark. Mid-point is white in light mode and the panel colour in dark,
 * otherwise "no change" cells glow against a dark background. */
function diverging(v, max) {
  const t = Math.max(-1, Math.min(1, v / (max || 1)));
  const mid = C.bg;
  const end = t >= 0 ? C.vermilion : C.accent;
  return mix(mid, end, Math.abs(t));
}
function mix(a, b, t) {
  const p = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16));
  const [ar, ag, ab] = p(a.length === 4 ? `#${a[1]}${a[1]}${a[2]}${a[2]}${a[3]}${a[3]}` : a);
  const [br, bg, bb] = p(b.length === 4 ? `#${b[1]}${b[1]}${b[2]}${b[2]}${b[3]}${b[3]}` : b);
  const c = (x, y) => Math.round(x + (y - x) * t);
  return `rgb(${c(ar, br)},${c(ag, bg)},${c(ab, bb)})`;
}

const tip = () => {
  let el = document.getElementById("tooltip");
  if (!el) {
    el = document.createElement("div");
    el.id = "tooltip";
    el.className = "tooltip";
    document.body.append(el);
  }
  return el;
};

async function buildHeatmap() {
  ladderRows = await json("ladder.json");
  drawHeatmap();
}

function drawHeatmap() {
  const svg = document.getElementById("heatmap");
  if (!svg || !ladderRows.length) return;
  clear(svg);

  const setOrder = ["carbon_starvation_DIN", "photorespiration_core", "ath00630",
                    "fermentation", "hypoxia_responsive", "ath00500", "ath00010",
                    "ath00710", "photosynthesis_apparatus", "rubisco"];
  const studies = [...new Map(ladderRows.map((r) => [r.accession, r])).values()]
    .sort((a, b) => (a.light === b.light ? 0 : a.light === "light" ? -1 : 1));
  const byKey = new Map(ladderRows.map((r) => [`${r.set}|${r.accession}`, r]));
  const sets = setOrder.filter((s) => ladderRows.some((r) => r.set === s));

  const w = 620, pad = { l: 168, r: 16, t: 34, b: 52 };
  const cw = (w - pad.l - pad.r) / studies.length;
  const ch = 21;
  const h = pad.t + sets.length * ch + pad.b;
  svg.setAttribute("viewBox", `0 0 ${w} ${h}`);

  const max = Math.max(...ladderRows.map((r) => Math.abs(r.shift || 0)));

  studies.forEach((st, i) => {
    const label = st.accession.replace("OSD-", "").replace("-light", " lit")
      .replace("-dark", " dark");
    svg.append(el("text", {
      x: pad.l + i * cw + cw / 2, y: pad.t - 10, "text-anchor": "middle",
      "font-size": 9.5, "font-weight": 600,
      fill: st.light === "light" ? C.orange : C.accent,
    }, label));
  });

  sets.forEach((name, r) => {
    const y = pad.t + r * ch;
    const g = el("g", { class: "hm-row", style: "cursor:pointer" });
    g.append(el("text", {
      x: pad.l - 8, y: y + ch / 2 + 3.5, "text-anchor": "end", "font-size": 9, fill: C.ink,
    }, PRETTY[name] || name));

    studies.forEach((st, i) => {
      const rec = byKey.get(`${name}|${st.accession}`);
      const cell = el("rect", {
        x: pad.l + i * cw + 0.6, y: y + 0.6, width: cw - 1.2, height: ch - 1.2, rx: 2,
        fill: rec && rec.shift !== null ? diverging(rec.shift, max) : "none",
        stroke: C.line, "stroke-width": 0.6,
      });
      if (rec) {
        cell.addEventListener("mousemove", (e) => {
          const t = tip();
          t.innerHTML =
            `<b>${PRETTY[name] || name}</b><br>${st.accession} · ${st.hardware} · ` +
            `${st.light}<br>shift ${rec.shift >= 0 ? "+" : ""}${(+rec.shift).toFixed(4)}` +
            ` · n=${rec.n} · p=${rec.p === null ? "n/a" : (+rec.p).toExponential(2)}`;
          t.style.display = "block";
          t.style.left = `${e.pageX + 12}px`;
          t.style.top = `${e.pageY - 10}px`;
        });
        cell.addEventListener("mouseleave", () => { tip().style.display = "none"; });
      }
      g.append(cell);
    });

    // The row doubles as a control for the volcano above — the discoverable route to
    // the thing that was broken.
    g.addEventListener("click", () => {
      if (window.__setHighlight) window.__setHighlight(name);
    });
    svg.append(g);
  });

  // legend
  const lw = 120, lx = pad.l, ly = pad.t + sets.length * ch + 26;
  for (let i = 0; i <= 40; i++) {
    const v = (i / 40) * 2 - 1;
    svg.append(el("rect", { x: lx + (i / 41) * lw, y: ly - 8, width: lw / 41 + 0.6,
                            height: 9, fill: diverging(v * max, max) }));
  }
  svg.append(el("text", { x: lx - 6, y: ly, "text-anchor": "end", "font-size": 8.5,
                          fill: C.faint }, `−${max.toFixed(1)}`));
  svg.append(el("text", { x: lx + lw + 6, y: ly, "font-size": 8.5, fill: C.faint },
                         `+${max.toFixed(1)}`));
  svg.append(el("text", { x: lx + lw + 52, y: ly, "font-size": 8.5, fill: C.soft },
                         "shift vs all other genes (log₂FC) · click a row to highlight it above"));
}

/* ────────────────────────────────── sankey: transcript -> protein -> pathway */

let sankeyData = null, sankeyN = 10, sankeyBroad = false;

async function buildSankey() {
  sankeyData = await json("sankey.json");
  const sel = document.getElementById("sankeyN");
  if (sel) {
    sel.innerHTML = [6, 8, 10, 14, 20, 31]
      .filter((n) => n <= sankeyData.pathways.length)
      .map((n) => `<option value="${n}"${n === sankeyN ? " selected" : ""}>` +
                  `top ${n} pathways</option>`).join("");
    sel.addEventListener("change", () => { sankeyN = +sel.value; drawSankey(); });
  }
  const broad = document.getElementById("sankeyBroad");
  if (broad) {
    broad.checked = sankeyBroad;
    broad.addEventListener("change", () => {
      sankeyBroad = broad.checked; drawSankey();
    });
  }
  drawSankey();
}

const BIN_COLOUR = { up: () => C.vermilion, down: () => C.accent,
                     "not significant": () => C.grey,
                     "not measured": () => C.line };
// Pathways come from two databases with two membership sources, so the column says which.
// Deliberately NOT --green: the shared CoSE kit retints that token onto --accent, which is
// also the "down" bin, so a green pathway node renders the same blue as the column beside
// it. Orange and purple are the two Okabe-Ito tokens the kit leaves alone.
const DB_COLOUR = { KEGG: () => C.orange, MapMan: () => C.purple };

function drawSankey() {
  const svg = document.getElementById("sankey");
  if (!svg || !sankeyData) return;
  clear(svg);
  const d = sankeyData;
  const v = d.views[sankeyBroad ? "all" : "core"];
  const eligible = d.pathways.filter((p) => sankeyBroad || !p.broad);
  const paths = eligible.slice(0, sankeyN);
  const keep = new Set(paths.map((p) => p.name));
  const byName = new Map(d.pathways.map((p) => [p.name, p]));

  const w = 620, pad = { l: 96, r: 210, t: 26, b: 20 };
  const colX = [pad.l, (pad.l + (w - pad.r)) / 2, w - pad.r];
  const gap = 7;

  // Column 3 totals are memberships, not genes, so each column is scaled to its own sum.
  const p2p = d.pr_to_path.filter((l) => keep.has(l.to));
  const pathTotals = new Map(paths.map((p) => [p.name,
    p2p.filter((l) => l.to === p.name).reduce((a, l) => a + l.n, 0)]));

  const cols = [
    d.bins.map((b) => ({ key: b, label: b, n: v.tx_bins[b] || 0 })),
    d.bins.map((b) => ({ key: b, label: b, n: v.pr_bins[b] || 0 })),
    paths.map((p) => ({ key: p.name, label: p.name, n: pathTotals.get(p.name) || 0 })),
  ].map((nodes) => nodes.filter((n) => n.n > 0));

  const h = Math.max(300, 34 + Math.max(...cols.map((c) => c.length)) * 26);
  svg.setAttribute("viewBox", `0 0 ${w} ${h}`);
  const usable = h - pad.t - pad.b;

  // Lay out each column: node height proportional to its share of that column's total,
  // on top of a 3px floor so a 17-locus pathway beside a 4,812-locus one is still a
  // visible target. The floors are subtracted from the free space before it is shared
  // out — take them afterwards and the column overruns its own box, which it did.
  const FLOOR = 3;
  const layout = cols.map((nodes) => {
    const total = nodes.reduce((a, n) => a + n.n, 0) || 1;
    const free = Math.max(0, usable - gap * (nodes.length - 1) - FLOOR * nodes.length);
    let y = pad.t;
    return nodes.map((n) => {
      const nh = FLOOR + (n.n / total) * free;
      const box = { ...n, y, h: nh, inUsed: 0, outUsed: 0 };
      y += nh + gap;
      return box;
    });
  });
  const find = (ci, key) => layout[ci].find((n) => n.key === key);

  const ribbon = (x0, y0, h0, x1, y1, h1, colour, title) => {
    const cx = (x0 + x1) / 2;
    const dd = `M${x0},${y0} C${cx},${y0} ${cx},${y1} ${x1},${y1} ` +
               `L${x1},${y1 + h1} C${cx},${y1 + h1} ${cx},${y0 + h0} ${x0},${y0 + h0} Z`;
    const p = el("path", { d: dd, fill: colour, "fill-opacity": 0.34 });
    p.addEventListener("mousemove", (e) => {
      const t = tip(); t.innerHTML = title; t.style.display = "block";
      t.style.left = `${e.pageX + 12}px`; t.style.top = `${e.pageY - 10}px`;
      p.setAttribute("fill-opacity", 0.62);
    });
    p.addEventListener("mouseleave", () => {
      tip().style.display = "none"; p.setAttribute("fill-opacity", 0.34);
    });
    svg.append(p);
  };

  const nodeW = 11;
  for (const l of v.tx_to_pr) {
    const a = find(0, l.from), b = find(1, l.to);
    if (!a || !b) continue;
    const hh = (l.n / a.n) * a.h, hb = (l.n / b.n) * b.h;
    ribbon(colX[0] + nodeW, a.y + a.outUsed, hh, colX[1], b.y + b.inUsed, hb,
           BIN_COLOUR[l.from](), `${l.from} → ${l.to}<br><b>${l.n}</b> loci`);
    a.outUsed += hh; b.inUsed += hb;
  }
  for (const l of p2p) {
    const a = find(1, l.from), b = find(2, l.to);
    if (!a || !b) continue;
    const hh = (l.n / a.n) * a.h, hb = (l.n / b.n) * b.h;
    ribbon(colX[1] + nodeW, a.y + a.outUsed, hh, colX[2], b.y + b.inUsed, hb,
           BIN_COLOUR[l.from](), `${l.from} → ${l.to}<br><b>${l.n}</b> loci`);
    a.outUsed += hh; b.inUsed += hb;
  }

  layout.forEach((nodes, ci) => {
    for (const n of nodes) {
      const meta = ci === 2 ? byName.get(n.key) : null;
      const rect = el("rect", {
        x: colX[ci], y: n.y, width: nodeW, height: n.h, rx: 2,
        fill: meta ? DB_COLOUR[meta.db]() : BIN_COLOUR[n.key](),
      });
      svg.append(rect);
      if (meta) {
        // PaintOmics does not publish MapMan membership, so it is reconstructed from the
        // vendored diagrams. Where that disagrees with the server's own feature count the
        // node says so rather than presenting one number as settled.
        const disagrees = Math.abs(meta.measured - meta.features) > 3;
        const t = `<b>${meta.name}</b><br>${meta.db} · p = ${meta.p.toExponential(1)}` +
                  `<br>${n.n.toLocaleString()} loci in this view` +
                  `<br>${meta.measured.toLocaleString()} of ${meta.size.toLocaleString()} ` +
                  `members measured here` +
                  (disagrees ? `<br><i>PaintOmics placed ${meta.features.toLocaleString()}` +
                               ` features on this map</i>` : "") +
                  (meta.broad ? `<br><i>whole-ontology map: ` +
                                `${Math.round(100 * meta.measured / d.n_measured)} % of all ` +
                                `measured loci</i>` : "");
        rect.addEventListener("mousemove", (e) => {
          const tt = tip(); tt.innerHTML = t; tt.style.display = "block";
          tt.style.left = `${e.pageX + 12}px`; tt.style.top = `${e.pageY - 10}px`;
        });
        rect.addEventListener("mouseleave", () => { tip().style.display = "none"; });
        // A pathway node is a control, the same one the heatmap rows and the dropdown
        // drive. Broad maps are absent from the highlight index on purpose (an overlay
        // marking 44 % of the points marks nothing), so they stay inert.
        if (!meta.broad) {
          rect.style.cursor = "pointer";
          rect.addEventListener("click", () => window.__setHighlight(meta.name));
        }
      }
      const right = ci === 2;
      svg.append(el("text", {
        x: right ? colX[ci] + nodeW + 6 : colX[ci] - 6,
        y: n.y + n.h / 2 + 3, "text-anchor": right ? "start" : "end",
        "font-size": 8.6, fill: C.ink,
      }, `${n.label.length > 30 ? `${n.label.slice(0, 28)}…` : n.label} (${n.n})`));
    }
  });

  ["Transcriptome", "Proteome", "Pathway"].forEach((t, i) => {
    svg.append(el("text", { x: colX[i] + (i === 2 ? nodeW : 0), y: pad.t - 11,
                            "text-anchor": i === 2 ? "start" : "middle",
                            "font-size": 9.5, "font-weight": 600, fill: C.soft }, t));
  });

  // Database key, under the pathway column heading.
  let kx = colX[2] + nodeW;
  for (const db of ["KEGG", "MapMan"]) {
    const n = paths.filter((p) => byName.get(p.name).db === db).length;
    if (!n) continue;
    svg.append(el("rect", { x: kx, y: pad.t - 6, width: 7, height: 7, rx: 1.5,
                            fill: DB_COLOUR[db]() }));
    svg.append(el("text", { x: kx + 10, y: pad.t, "font-size": 8, fill: C.faint },
                           `${db} (${n})`));
    kx += 62;
  }

  const cap = document.getElementById("sankeycap");
  if (cap) {
    // Columns 1 and 2 cover every pathway in the current view, not just the top N drawn
    // in column 3, so the caption has to name both numbers or the totals read as wrong.
    const nk = eligible.filter((p) => p.db === "KEGG").length;
    const nm = eligible.length - nk;
    const hidden = d.pathways.filter((p) => p.broad).map((p) => p.name);
    cap.textContent =
      `${v.n_genes.toLocaleString()} loci sit in at least one of the ${eligible.length} ` +
      `pathways in this view (${nk} KEGG, ${nm} MapMan), out of the ` +
      `${d.n_pathways_total} PaintOmics called significant; the pathway column draws the ` +
      `${paths.length} most significant of them. Transcript and protein columns conserve ` +
      `loci exactly; the pathway column counts memberships, so the ` +
      `${(v.n_memberships - v.n_genes).toLocaleString()} loci sitting in more than one ` +
      `pathway are counted more than once. KEGG membership comes from the KEGG REST API; ` +
      `MapMan publishes none, so it is reconstructed from the diagram layouts — matching ` +
      `PaintOmics' own feature count for 9 of the 13 MapMan maps and running larger on the ` +
      `other 4, with both numbers in every node's tooltip. ` +
      (sankeyBroad
        ? `Whole-ontology maps are included: ${hidden.join(", ")} each cover more than ` +
          `${Math.round(100 * d.broad_map_share)} % of the measured loci, so their ribbons ` +
          `are close to the marginal distribution rather than a statement about them.`
        : `${hidden.length} whole-ontology MapMan maps — ${hidden.join(", ")} — are off by ` +
          `default because each covers more than ${Math.round(100 * d.broad_map_share)} % ` +
          `of the measured loci; tick the box to include them.`);
  }
}
