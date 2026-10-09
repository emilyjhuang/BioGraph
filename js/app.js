/* BioGraph UI: D3 force-directed rendering and interaction. */
(() => {
  const D = window.BIOGRAPH_DATA;
  const A = window.BioGraphAnalysis;

  // ---------- model ----------
  const g = A.buildGraph(D.nodes, D.edges);
  const metrics = A.centralities(g);
  const loops = A.feedbackLoops(g);
  const byId = new Map(D.nodes.map((n) => [n.id, n]));
  const pathwayIdx = new Map(D.PATHWAYS.map((p, i) => [p.id, i]));
  const pathwayById = new Map(D.PATHWAYS.map((p) => [p.id, p]));
  const drugsByTarget = new Map();
  D.DRUGS.forEach((d) => d.targets.forEach((t) => {
    if (!drugsByTarget.has(t)) drugsByTarget.set(t, []);
    drugsByTarget.get(t).push(d);
  }));
  const rank = (key) => {
    const ids = D.nodes.map((n) => n.id).sort((a, b) => metrics[b][key] - metrics[a][key]);
    return new Map(ids.map((id, i) => [id, i + 1]));
  };
  const ranks = { betweenness: rank("betweenness"), pagerank: rank("pagerank"), degree: rank("degree") };
  const sortedNodes = D.nodes.slice().sort((a, b) => a.label.localeCompare(b.label));

  const state = {
    colorBy: "pathway",
    sizeBy: "betweenness",
    hidden: new Set(),
    selected: null,
    hover: null,
    highlight: null,     // { nodes:Set<id>, edges:Set<idx> }
    sim: null,           // { values:{id:x}, frames, clamped:Set }
    simFrame: null,      // Float64Array currently displayed
    perturb: new Map(),  // id -> ±1
    scenario: null,
  };

  // ---------- helpers ----------
  const $ = (s) => document.querySelector(s);
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  const glyph = (sign) => (sign > 0 ? "→" : "⊣");
  const fmt = (v, d = 2) => (v >= 0 ? "+" : "−") + Math.abs(v).toFixed(d);
  const pathwayColor = (pid) => css(`--s${pathwayIdx.get(pid) + 1}`);
  const lbl = (id) => byId.get(id).label;

  let palette = {};
  function readPalette() {
    palette = {
      edge: css("--edge"), text: css("--text-primary"), muted: css("--text-muted"),
      neg: css("--div-neg"), mid: css("--div-mid"), pos: css("--div-pos"),
      onc: css("--onc"), tsg: css("--tsg"), surface2: css("--surface-2"),
    };
    palette.diverge = (v) => (v < 0 ? d3.interpolateRgb(palette.mid, palette.neg)(Math.min(1, -v / 0.7))
                                     : d3.interpolateRgb(palette.mid, palette.pos)(Math.min(1, v / 0.7)));
  }
  readPalette();

  // ---------- header stats ----------
  const negLoops = loops.filter((l) => l.sign < 0).length;
  $("#stats").innerHTML = [
    [D.nodes.length, "proteins"], [D.edges.length, "interactions"], [D.PATHWAYS.length, "pathways"],
    [loops.length, "feedback loops"], [D.DRUGS.length, "drugs"],
  ].map(([v, k]) => `<div class="stat"><b>${v}</b><span>${k}</span></div>`).join("");
  document.querySelectorAll("[data-stat=nodes]").forEach((el) => (el.textContent = D.nodes.length));
  document.querySelectorAll("[data-stat=edges]").forEach((el) => (el.textContent = D.edges.length));
  $("#keggList").textContent = D.PATHWAYS.map((p) => `${p.kegg} ${p.name}`).join("; ");

  // ---------- SVG scaffolding ----------
  const svg = d3.select("#graph");
  const defs = svg.append("defs");
  const root = svg.append("g");
  const hullLayer = root.append("g");
  const linkLayer = root.append("g");
  const nodeLayer = root.append("g");
  const hullLabelLayer = root.append("g");

  const markerVariants = ["base", "hi", "pos", "neg"];
  function buildMarkers() {
    defs.selectAll("*").remove();
    const colors = { base: palette.edge, hi: palette.text, pos: palette.pos, neg: palette.neg };
    for (const v of markerVariants) {
      defs.append("marker").attr("id", `arr-${v}`).attr("viewBox", "0 -5 10 10").attr("refX", 9).attr("refY", 0)
        .attr("markerWidth", 7).attr("markerHeight", 7).attr("orient", "auto").attr("markerUnits", "userSpaceOnUse")
        .append("path").attr("d", "M0,-4L10,0L0,4Z").attr("fill", colors[v]);
      defs.append("marker").attr("id", `inh-${v}`).attr("viewBox", "0 -6 4 12").attr("refX", 3).attr("refY", 0)
        .attr("markerWidth", 4).attr("markerHeight", 10).attr("orient", "auto").attr("markerUnits", "userSpaceOnUse")
        .append("rect").attr("x", 0).attr("y", -5).attr("width", 2.5).attr("height", 10).attr("fill", colors[v]);
    }
  }
  buildMarkers();

  // ---------- simulation data ----------
  const simNodes = D.nodes.map((n) => ({ ...n }));
  const simLinks = D.edges.map((e, i) => ({ ...e, i }));
  const nodeById = new Map(simNodes.map((n) => [n.id, n]));

  const metricExtent = {
    betweenness: d3.max(D.nodes, (n) => metrics[n.id].betweenness),
    pagerank: d3.max(D.nodes, (n) => metrics[n.id].pagerank),
    degree: d3.max(D.nodes, (n) => metrics[n.id].degree),
  };
  const radius = (n) => {
    if (state.sizeBy === "none") return 8;
    return 5 + 13 * Math.sqrt(metrics[n.id][state.sizeBy] / metricExtent[state.sizeBy]);
  };

  // Pathway anchor points on an ellipse
  const anchors = new Map();
  function layoutAnchors(w, h) {
    const order = ["rtk", "pi3k", "wnt", "cycle", "p53", "apop", "nfkb"];
    order.forEach((pid, i) => {
      const a = -Math.PI / 2 + (i / order.length) * 2 * Math.PI;
      anchors.set(pid, { x: Math.cos(a) * w * 0.4, y: Math.sin(a) * h * 0.4 });
    });
  }
  layoutAnchors(1100, 800);

  const force = d3.forceSimulation(simNodes)
    .force("link", d3.forceLink(simLinks).id((d) => d.id).distance((l) => (l.source.pathway === l.target.pathway ? 42 : 90)).strength(0.35))
    .force("charge", d3.forceManyBody().strength(-210).distanceMax(400))
    .force("collide", d3.forceCollide().radius((d) => radius(d) + 9))
    .force("x", d3.forceX((d) => anchors.get(d.pathway).x).strength(0.14))
    .force("y", d3.forceY((d) => anchors.get(d.pathway).y).strength(0.14))
    .stop();
  for (let i = 0; i < 260; i++) force.tick();

  // ---------- links ----------
  const linkSel = linkLayer.selectAll("path").data(simLinks).join("path")
    .attr("class", (d) => "link" + (d.indirect ? " indirect" : ""));

  function linkPath(d) {
    const s = d.source, t = d.target;
    const dx = t.x - s.x, dy = t.y - s.y, dist = Math.hypot(dx, dy) || 1;
    const k = 0.14;
    const cx = (s.x + t.x) / 2 - dy * k, cy = (s.y + t.y) / 2 + dx * k;
    const towards = (px, py, qx, qy, r) => {
      const ddx = qx - px, ddy = qy - py, dd = Math.hypot(ddx, ddy) || 1;
      return [px + (ddx / dd) * r, py + (ddy / dd) * r];
    };
    const [sx, sy] = towards(s.x, s.y, cx, cy, radius(s));
    const [tx, ty] = towards(t.x, t.y, cx, cy, radius(t) + 2);
    return `M${sx},${sy}Q${cx},${cy} ${tx},${ty}`;
  }

  // ---------- nodes ----------
  const nodeSel = nodeLayer.selectAll("g").data(simNodes, (d) => d.id).join("g").attr("class", "node");
  nodeSel.append("circle").attr("class", "halo");
  nodeSel.append("circle").attr("class", "body");
  nodeSel.append("circle").attr("class", "ring").attr("opacity", 0);
  nodeSel.append("text").attr("text-anchor", "middle");

  function applySizes() {
    nodeSel.select("circle.body").attr("r", radius);
    nodeSel.select("circle.halo").attr("r", (d) => radius(d) + 4).attr("opacity", (d) => (drugsByTarget.has(d.id) ? 1 : 0));
    nodeSel.select("circle.ring").attr("r", (d) => radius(d) + 7);
    nodeSel.select("text").attr("dy", (d) => radius(d) + 12).text((d) => d.label);
    force.force("collide").radius((d) => radius(d) + 9);
  }
  applySizes();

  // ---------- hulls ----------
  const hullLine = d3.line().curve(d3.curveCatmullRomClosed.alpha(0.8));
  const hullSel = hullLayer.selectAll("path").data(D.PATHWAYS).join("path").attr("class", "hull");
  const hullLabelSel = hullLabelLayer.selectAll("text").data(D.PATHWAYS).join("text").attr("class", "hull-label").attr("text-anchor", "middle").text((p) => p.name);

  function updateHulls() {
    hullSel.attr("d", (p) => {
      if (state.hidden.has(p.id)) return null;
      const pts = [];
      simNodes.forEach((n) => {
        if (n.pathway !== p.id) return;
        const r = radius(n) + 18;
        for (let a = 0; a < 8; a++) pts.push([n.x + r * Math.cos((a * Math.PI) / 4), n.y + r * Math.sin((a * Math.PI) / 4)]);
      });
      const hull = d3.polygonHull(pts);
      return hull ? hullLine(hull) : null;
    });
    hullLabelSel.each(function (p) {
      const ns = simNodes.filter((n) => n.pathway === p.id);
      const x = d3.mean(ns, (n) => n.x), y = d3.min(ns, (n) => n.y - radius(n)) - 26;
      d3.select(this).attr("x", x).attr("y", y).attr("display", state.hidden.has(p.id) ? "none" : null);
    });
  }

  function ticked() {
    linkSel.attr("d", linkPath);
    nodeSel.attr("transform", (d) => `translate(${d.x},${d.y})`);
    updateHulls();
  }
  force.on("tick", ticked);
  ticked();

  // ---------- painting ----------
  function nodeFill(n, i) {
    if (state.simFrame) return palette.diverge(state.simFrame[i]);
    if (state.colorBy === "role") return n.role === "onc" ? palette.onc : n.role === "tsg" ? palette.tsg : palette.surface2;
    if (state.colorBy === "drug") return drugsByTarget.has(n.id) ? pathwayColor(n.pathway) : palette.surface2;
    return pathwayColor(n.pathway);
  }

  function focusSet() {
    // nodes/edges to emphasise; null = everything
    if (state.highlight) return state.highlight;
    const id = state.hover || state.selected;
    if (!id) return null;
    const i = g.index.get(id);
    const nodes = new Set([id]), edges = new Set();
    g.out[i].forEach((e) => { nodes.add(g.nodes[e.to].id); edges.add(e.edge); });
    g.inc[i].forEach((e) => { nodes.add(g.nodes[e.from].id); edges.add(e.edge); });
    return { nodes, edges };
  }

  function paint() {
    const focus = focusSet();
    const visible = (n) => !state.hidden.has(n.pathway);
    const frame = state.simFrame;

    nodeSel
      .attr("display", (d) => (visible(d) ? null : "none"))
      .classed("clamped", (d) => !!(state.sim && state.sim.clamped.has(d.id)))
      .attr("opacity", (d) => (!focus || focus.nodes.has(d.id) ? 1 : 0.13));
    nodeSel.select("circle.body").attr("fill", (d) => nodeFill(d, g.index.get(d.id)));
    nodeSel.select("circle.ring").attr("opacity", (d) => (d.id === state.selected ? 1 : 0));

    linkSel
      .attr("display", (d) => (visible(d.source) && visible(d.target) ? null : "none"))
      .each(function (d) {
        let color = palette.edge, variant = "base", opacity = 0.75, width = 1.2;
        if (frame) {
          const flow = d.sign * frame[g.index.get(d.source.id)];
          const mag = Math.min(1, Math.abs(flow) / 0.6);
          if (mag > 0.04) {
            variant = flow > 0 ? "pos" : "neg";
            color = flow > 0 ? palette.pos : palette.neg;
            opacity = 0.25 + 0.75 * mag; width = 1 + 2.2 * mag;
          } else opacity = 0.18;
        }
        if (focus) {
          if (focus.edges.has(d.i)) {
            if (!frame) { color = palette.text; variant = "hi"; }
            opacity = 1; width = Math.max(width, 2);
          } else opacity = 0.05;
        }
        d3.select(this).attr("stroke", color).attr("opacity", opacity).attr("stroke-width", width)
          .attr("marker-end", `url(#${d.sign > 0 ? "arr" : "inh"}-${variant})`);
      });

    hullSel.attr("fill", (p) => pathwayColor(p.id)).attr("fill-opacity", frame || focus ? 0.03 : 0.07)
      .attr("stroke", (p) => pathwayColor(p.id)).attr("stroke-opacity", frame || focus ? 0.1 : 0.25);
    hullLabelSel.attr("fill", (p) => pathwayColor(p.id)).attr("opacity", focus ? 0.35 : 0.9);
    nodeSel.select("circle.halo").attr("stroke", palette.muted);
    updateNodeLegend();
  }

  function updateNodeLegend() {
    const el = $("#nodeLegendExtra");
    if (state.simFrame) el.innerHTML = `<span class="sw" style="background:${palette.neg}"></span>down <span class="sw" style="background:${palette.pos}"></span>up`;
    else if (state.colorBy === "role") el.innerHTML = `<span class="sw" style="background:${palette.onc}"></span>oncogene <span class="sw" style="background:${palette.tsg}"></span>tumor suppressor`;
    else el.innerHTML = `<span class="halo"></span>drug target`;
  }

  // ---------- zoom ----------
  const zoom = d3.zoom().scaleExtent([0.2, 5]).on("zoom", (ev) => root.attr("transform", ev.transform));
  svg.call(zoom).on("dblclick.zoom", null);
  function fit(duration = 600, ids = null) {
    const ns = simNodes.filter((n) => (ids ? ids.has(n.id) : !state.hidden.has(n.pathway)));
    if (!ns.length) return;
    const { width, height } = svg.node().getBoundingClientRect();
    const detailOpen = !$("#detail").hidden && width > 900;
    const usableW = width - (detailOpen ? 340 : 0);
    const x0 = d3.min(ns, (n) => n.x) - 60, x1 = d3.max(ns, (n) => n.x) + 60;
    const y0 = d3.min(ns, (n) => n.y) - 110, y1 = d3.max(ns, (n) => n.y) + 50;
    const k = Math.min(ids ? 2 : 1.6, 0.95 / Math.max((x1 - x0) / usableW, (y1 - y0) / height));
    const t = d3.zoomIdentity.translate(usableW / 2 - (k * (x0 + x1)) / 2, height / 2 - (k * (y0 + y1)) / 2).scale(k);
    svg.transition().duration(duration).call(zoom.transform, t);
  }
  $("#zoomIn").onclick = () => svg.transition().call(zoom.scaleBy, 1.4);
  $("#zoomOut").onclick = () => svg.transition().call(zoom.scaleBy, 1 / 1.4);
  $("#zoomFit").onclick = () => fit();

  // ---------- drag ----------
  nodeSel.call(d3.drag()
    .on("start", (ev, d) => { if (!ev.active) force.alphaTarget(0.2).restart(); d.fx = d.x; d.fy = d.y; })
    .on("drag", (ev, d) => { d.fx = ev.x; d.fy = ev.y; })
    .on("end", (ev, d) => { if (!ev.active) force.alphaTarget(0); d.fx = null; d.fy = null; }));

  // ---------- tooltip & hover ----------
  const tip = $("#tooltip");
  nodeSel
    .on("mouseenter", (ev, d) => {
      state.hover = d.id; paint();
      const m = metrics[d.id];
      const simLine = state.sim ? `<div>Predicted change: <b>${fmt(state.sim.values[d.id])}</b>${state.sim.clamped.has(d.id) ? " (clamped)" : ""}</div>` : "";
      tip.innerHTML = `<b>${esc(d.label)}</b> <span class="sub mono">${esc(d.id)}</span>
        <div class="sub">${esc(pathwayById.get(d.pathway).name)} · ${esc(d.type)}</div>
        ${simLine}<div class="sub">in ${m.inDegree} · out ${m.outDegree} · betweenness #${ranks.betweenness.get(d.id)}</div>`;
      tip.hidden = false;
    })
    .on("mousemove", (ev) => {
      const r = $(".stage").getBoundingClientRect();
      let x = ev.clientX - r.left + 14, y = ev.clientY - r.top + 14;
      if (x + 270 > r.width) x -= 290;
      tip.style.left = x + "px"; tip.style.top = y + "px";
    })
    .on("mouseleave", () => { state.hover = null; tip.hidden = true; paint(); })
    .on("click", (ev, d) => {
      ev.stopPropagation();
      if (ev.shiftKey || currentTab === "paths") return pathPick(d.id);
      if (currentTab === "simulate" && $("#customBox").open) { $("#pertNode").value = d.id; }
      select(d.id);
    });
  svg.on("click", () => { state.highlight = null; clearCardHighlights(); select(null); });

  linkSel.on("mouseenter", (ev, d) => {
    tip.innerHTML = `<b>${esc(d.source.label)} ${glyph(d.sign)} ${esc(d.target.label)}</b>
      <div class="sub">${esc(d.mechanism)}${d.indirect ? " (indirect)" : ""}</div>${d.note ? `<div>${esc(d.note)}</div>` : ""}`;
    tip.hidden = false;
  }).on("mousemove", (ev) => {
    const r = $(".stage").getBoundingClientRect();
    tip.style.left = ev.clientX - r.left + 14 + "px"; tip.style.top = ev.clientY - r.top + 14 + "px";
  }).on("mouseleave", () => (tip.hidden = true)).style("cursor", "help");

  // ---------- detail drawer ----------
  function select(id) {
    state.selected = id;
    const el = $("#detail");
    if (!id) { el.hidden = true; paint(); return; }
    const n = byId.get(id), m = metrics[id], i = g.index.get(id);
    const nb = (list, dir) => list.map((e) => {
      const ed = D.edges[e.edge], other = g.nodes[dir === "in" ? e.from : e.to];
      return `<li data-id="${other.id}"><span class="g ${ed.sign < 0 ? "inh" : ""}">${glyph(ed.sign)}</span><b>${esc(other.label)}</b><span class="mech">${esc(ed.mechanism)}${ed.indirect ? "*" : ""}</span></li>`;
    }).join("") || `<li class="sub" style="color:var(--text-muted)">none</li>`;
    const drugs = drugsByTarget.get(id) || [];
    const roleTag = n.role === "onc" ? `<span class="tag onc">oncogene</span>` : n.role === "tsg" ? `<span class="tag tsg">tumor suppressor</span>` : "";
    const nLoops = loops.filter((l) => l.nodes.includes(id)).length;
    el.innerHTML = `
      <button class="btn ghost small close" aria-label="Close">✕</button>
      <h2>${esc(n.label)}</h2>
      <div class="gene">${esc(n.id)}</div>
      <div class="tags">
        <span class="tag" style="color:#fff;background:${pathwayColor(n.pathway)}">${esc(pathwayById.get(n.pathway).name)}</span>
        <span class="tag">${esc(n.type)}</span>${roleTag}
      </div>
      <p class="desc">${esc(n.desc)}</p>
      ${state.sim ? `<div class="metric" style="margin-bottom:10px"><b style="color:${palette.diverge(state.sim.values[id])}">${fmt(state.sim.values[id])}</b><span>predicted change in current simulation</span></div>` : ""}
      <div class="metrics">
        <div class="metric"><b>${m.inDegree}/${m.outDegree}</b><span>in / out degree</span></div>
        <div class="metric"><b>#${ranks.betweenness.get(id)}</b><span>betweenness</span></div>
        <div class="metric"><b>${nLoops}</b><span>feedback loops</span></div>
      </div>
      ${drugs.length ? `<h4>Targeted by</h4>${drugs.map((d) => `<div class="drug"><b>${esc(d.name)}</b><span>${esc(d.cls)} · ${d.status}</span></div>`).join("")}` : ""}
      <h4>Regulated by (${g.inc[i].length})</h4><ul class="nbrs">${nb(g.inc[i], "in")}</ul>
      <h4>Regulates (${g.out[i].length})</h4><ul class="nbrs">${nb(g.out[i], "out")}</ul>
      <div class="actions">
        <button class="btn small" data-act="inhibit">⊣ Simulate inhibition</button>
        <button class="btn small" data-act="activate">→ Simulate activation</button>
        <button class="btn small" data-act="downstream">Show downstream</button>
      </div>`;
    el.hidden = false;
    el.querySelector(".close").onclick = () => select(null);
    el.querySelectorAll(".nbrs li[data-id]").forEach((li) => (li.onclick = () => { select(li.dataset.id); centerOn(li.dataset.id); }));
    el.querySelector("[data-act=inhibit]").onclick = () => runPerturbation({ [id]: -1 }, `${n.label} inhibited`);
    el.querySelector("[data-act=activate]").onclick = () => runPerturbation({ [id]: 1 }, `${n.label} activated`);
    el.querySelector("[data-act=downstream]").onclick = () => {
      const ds = new Set(A.downstream(g, id)); ds.add(id);
      const edges = new Set(D.edges.map((e, k) => (ds.has(e.source) && ds.has(e.target) ? k : -1)).filter((k) => k >= 0));
      state.highlight = { nodes: ds, edges }; paint();
    };
    paint();
  }

  function centerOn(id) {
    const n = nodeById.get(id);
    const { width, height } = svg.node().getBoundingClientRect();
    const k = Math.max(d3.zoomTransform(svg.node()).k, 1.2);
    const offset = width > 900 && !$("#detail").hidden ? 170 : 0;
    svg.transition().duration(600).call(zoom.transform, d3.zoomIdentity.translate(width / 2 - offset - k * n.x, height / 2 - k * n.y).scale(k));
  }

  // ---------- tabs ----------
  let currentTab = "explore";
  document.querySelectorAll(".tabs button").forEach((b) => (b.onclick = () => {
    currentTab = b.dataset.tab;
    document.querySelectorAll(".tabs button").forEach((x) => x.setAttribute("aria-selected", x === b));
    document.querySelectorAll(".tab-body").forEach((s) => (s.hidden = s.dataset.body !== currentTab));
    if (currentTab !== "paths" && currentTab !== "loops") { state.highlight = null; clearCardHighlights(); paint(); }
  }));
  const clearCardHighlights = () => document.querySelectorAll(".result-card.on, .looplist li.on").forEach((e) => e.classList.remove("on"));

  // ---------- explore tab ----------
  $("#nodeList").innerHTML = sortedNodes.map((n) => `<option value="${esc(n.label)}">${esc(n.id)}</option>`).join("");
  $("#search").addEventListener("change", (ev) => {
    const q = ev.target.value.trim().toLowerCase();
    const hit = D.nodes.find((n) => n.label.toLowerCase() === q || n.id.toLowerCase() === q)
      || D.nodes.find((n) => n.label.toLowerCase().startsWith(q) || n.id.toLowerCase().startsWith(q));
    if (hit) { state.hidden.delete(hit.pathway); renderLegend(); select(hit.id); centerOn(hit.id); }
  });

  function segmented(elId, onChange) {
    const el = $(elId);
    el.querySelectorAll("button").forEach((b) => (b.onclick = () => {
      el.querySelectorAll("button").forEach((x) => x.classList.toggle("on", x === b));
      onChange(b.dataset.v);
    }));
  }
  segmented("#colorBy", (v) => { state.colorBy = v; paint(); });
  $("#sizeBy").onchange = (ev) => { state.sizeBy = ev.target.value; applySizes(); force.alpha(0.3).restart(); renderHubs(); };

  function renderLegend() {
    $("#pathwayLegend").innerHTML = D.PATHWAYS.map((p) => {
      const count = D.nodes.filter((n) => n.pathway === p.id).length;
      return `<li data-id="${p.id}" class="${state.hidden.has(p.id) ? "off" : ""}">
        <span class="sw" style="background:${pathwayColor(p.id)}"></span>${esc(p.name)}
        <span class="kegg">${p.kegg}</span><span class="count">${count}</span></li>`;
    }).join("");
    document.querySelectorAll("#pathwayLegend li").forEach((li) => (li.onclick = () => {
      const id = li.dataset.id;
      state.hidden.has(id) ? state.hidden.delete(id) : state.hidden.add(id);
      renderLegend(); paint(); updateHulls();
    }));
  }
  renderLegend();

  function renderHubs() {
    const key = state.sizeBy === "none" ? "betweenness" : state.sizeBy;
    const top = D.nodes.slice().sort((a, b) => metrics[b.id][key] - metrics[a.id][key]).slice(0, 10);
    const max = metrics[top[0].id][key];
    document.querySelector("#hubList").previousElementSibling.innerHTML = `Top nodes <small>(${key})</small>`;
    $("#hubList").innerHTML = top.map((n) => {
      const v = metrics[n.id][key];
      return `<li data-id="${n.id}"><b>${esc(n.label)}</b><span class="bar" style="width:${(100 * v) / max}%;background:${pathwayColor(n.pathway)}"></span><span class="val">${key === "degree" ? v : v.toFixed(3)}</span></li>`;
    }).join("");
    document.querySelectorAll("#hubList li").forEach((li) => (li.onclick = () => { select(li.dataset.id); centerOn(li.dataset.id); }));
  }
  renderHubs();

  // ---------- simulate tab ----------
  let animTimer = null;
  function runPerturbation(perturb, title) {
    const res = A.propagate(g, perturb);
    state.sim = { ...res, clamped: new Set(Object.keys(perturb)), title };
    state.highlight = null;
    $("#simBanner").hidden = false;
    animate();
    renderSimResults();
    if (currentTab !== "simulate") document.querySelector('.tabs button[data-tab="simulate"]').click();
    if (state.selected) select(state.selected);
  }

  function animate() {
    if (animTimer) animTimer.stop();
    const frames = state.sim.frames;
    const duration = 2600;
    const ease = d3.easeCubicOut;
    animTimer = d3.timer((elapsed) => {
      const t = Math.min(1, elapsed / duration);
      const idx = Math.min(frames.length - 1, Math.round(ease(t) * (frames.length - 1)));
      state.simFrame = frames[idx];
      $("#simBannerText").innerHTML = `<b>${esc(state.sim.title)}</b> <span class="step">step ${idx}/${frames.length - 1}${state.sim.converged && idx === frames.length - 1 ? " · steady state" : ""}</span>`;
      paint();
      if (t >= 1) animTimer.stop();
    });
  }

  function renderSimResults() {
    const vals = Object.entries(state.sim.values).filter(([, v]) => Math.abs(v) > 0.01);
    vals.sort((a, b) => Math.abs(b[1]) - Math.abs(a[1]));
    const top = vals.filter(([id]) => !state.sim.clamped.has(id)).slice(0, 16);
    const affected = vals.length - state.sim.clamped.size;
    $("#simMeta").textContent = `· ${affected} nodes affected · ${state.sim.steps} iterations`;
    $("#simBars").innerHTML = top.map(([id, v]) => {
      const w = Math.min(50, Math.abs(v) * 50);
      const style = v >= 0 ? `left:50%;width:${w}%` : `right:50%;width:${w}%`;
      return `<li data-id="${id}" title="${esc(byId.get(id).id)}: ${fmt(v, 3)}"><span class="name">${esc(lbl(id))}</span>
        <span class="track"><span class="fill" style="${style};background:${palette.diverge(v * 1.5)}"></span></span>
        <span class="val">${fmt(v)}</span></li>`;
    }).join("");
    document.querySelectorAll("#simBars li").forEach((li) => (li.onclick = () => { select(li.dataset.id); centerOn(li.dataset.id); }));
    $("#simResults").hidden = false;
  }

  function exitSim() {
    if (animTimer) animTimer.stop();
    state.sim = null; state.simFrame = null; state.scenario = null;
    $("#simBanner").hidden = true; $("#simResults").hidden = true; $("#scenarioText").hidden = true;
    document.querySelectorAll(".scenario").forEach((b) => b.classList.remove("on"));
    paint();
    if (state.selected) select(state.selected);
  }
  $("#simExit").onclick = exitSim;
  $("#simReplay").onclick = () => state.sim && animate();

  $("#scenarios").innerHTML = D.SCENARIOS.map((s) => {
    const tag = Object.entries(s.perturb).map(([id, v]) => `${lbl(id)}${v < 0 ? "↓" : "↑"}`).join(" ");
    return `<button class="scenario" data-id="${s.id}"><b>${esc(s.title)}</b><span class="tag">${esc(tag)}</span></button>`;
  }).join("");
  document.querySelectorAll(".scenario").forEach((b) => (b.onclick = () => {
    const s = D.SCENARIOS.find((x) => x.id === b.dataset.id);
    document.querySelectorAll(".scenario").forEach((x) => x.classList.toggle("on", x === b));
    state.scenario = s;
    runPerturbation(s.perturb, s.title);
    const box = $("#scenarioText");
    box.innerHTML = `${esc(s.text)}<div class="watch">${s.watch.map((id) => {
      const v = state.sim.values[id];
      return `<span class="chip" data-id="${id}"><span class="dot" style="background:${palette.diverge(v * 1.5)}"></span>${esc(lbl(id))} <span class="v">${v > 0 ? "↑" : "↓"} ${fmt(v)}</span></span>`;
    }).join("")}</div>`;
    box.hidden = false;
    box.querySelectorAll(".chip").forEach((c) => (c.onclick = () => { select(c.dataset.id); centerOn(c.dataset.id); }));
    select(null);
    fit();
  }));

  $("#pertNode").innerHTML = sortedNodes.map((n) => `<option value="${n.id}">${esc(n.label)}${n.label !== n.id ? ` (${esc(n.id)})` : ""}</option>`).join("");
  $("#drugSelect").innerHTML = `<option value="">Choose a drug…</option>` + D.DRUGS.map((d) =>
    `<option value="${d.id}">${esc(d.name)} — ${esc(d.cls)}${d.status !== "approved" ? " (investigational)" : ""}</option>`).join("");

  function renderPertChips() {
    $("#pertChips").innerHTML = [...state.perturb].map(([id, v]) =>
      `<span class="chip"><span class="dot" style="background:${v < 0 ? palette.neg : palette.pos}"></span>${esc(lbl(id))} ${v < 0 ? "⊣" : "→"}<button data-id="${id}" aria-label="Remove">✕</button></span>`).join("")
      || `<span style="color:var(--text-muted);font-size:12px">No perturbations yet.</span>`;
    document.querySelectorAll("#pertChips button").forEach((b) => (b.onclick = () => { state.perturb.delete(b.dataset.id); renderPertChips(); }));
  }
  renderPertChips();
  $("#addInhibit").onclick = () => { state.perturb.set($("#pertNode").value, -1); renderPertChips(); };
  $("#addActivate").onclick = () => { state.perturb.set($("#pertNode").value, 1); renderPertChips(); };
  $("#drugSelect").onchange = (ev) => {
    const d = D.DRUGS.find((x) => x.id === ev.target.value);
    if (d) d.targets.forEach((t) => state.perturb.set(t, -1));
    ev.target.value = ""; renderPertChips();
  };
  $("#runSim").onclick = () => {
    if (!state.perturb.size) return;
    document.querySelectorAll(".scenario").forEach((x) => x.classList.remove("on"));
    $("#scenarioText").hidden = true;
    const title = [...state.perturb].map(([id, v]) => `${lbl(id)}${v < 0 ? "↓" : "↑"}`).join(" + ");
    runPerturbation(Object.fromEntries(state.perturb), title);
  };
  $("#clearSim").onclick = () => { state.perturb.clear(); renderPertChips(); exitSim(); };

  // ---------- paths tab ----------
  const opts = sortedNodes.map((n) => `<option value="${n.id}">${esc(n.label)}</option>`).join("");
  $("#pathFrom").innerHTML = opts; $("#pathTo").innerHTML = opts;
  $("#pathFrom").value = "EGFR"; $("#pathTo").value = "BAD";
  let pathPickNext = "from";
  function pathPick(id) {
    $(pathPickNext === "from" ? "#pathFrom" : "#pathTo").value = id;
    if (pathPickNext === "to") findPaths();
    pathPickNext = pathPickNext === "from" ? "to" : "from";
    if (currentTab !== "paths") document.querySelector('.tabs button[data-tab="paths"]').click();
  }
  const chainHtml = (ids, edges) => ids.map((id, k) => {
    if (k === 0) return `<span class="n">${esc(lbl(id))}</span>`;
    const s = D.edges[edges[k - 1]].sign;
    return `<span class="g ${s < 0 ? "inh" : ""}">${glyph(s)}</span><span class="n">${esc(lbl(id))}</span>`;
  }).join("");
  const effectBadge = (sign) => sign > 0 ? `<span class="badge pos">net activation</span>` : `<span class="badge neg">net inhibition</span>`;

  function highlightPath(nodes, edges) {
    state.highlight = { nodes: new Set(nodes), edges: new Set(edges) };
    if (state.selected) select(null);
    nodes.forEach((id) => state.hidden.delete(byId.get(id).pathway));
    renderLegend(); paint(); fit(600, state.highlight.nodes);
  }

  function findPaths() {
    const from = $("#pathFrom").value, to = $("#pathTo").value;
    const out = $("#pathResults");
    if (from === to) { out.innerHTML = `<p class="summary">Pick two different proteins.</p>`; return; }
    const res = A.shortestPaths(g, from, to, { undirected: $("#pathUndirected").checked });
    if (!res.paths.length) {
      out.innerHTML = `<p class="summary">No directed path from <b>${esc(lbl(from))}</b> to <b>${esc(lbl(to))}</b>. Information doesn't flow that way in this network. Try ignoring edge direction.</p>`;
      state.highlight = null; paint(); return;
    }
    const signs = new Set(res.paths.map((p) => p.sign));
    const verdict = signs.size > 1
      ? `<span class="badge mixed">mixed</span> Shortest routes disagree in sign, so the net effect depends on context.`
      : `${effectBadge([...signs][0])} Every shortest route agrees.`;
    out.innerHTML = `<p class="summary"><b>${res.paths.length}</b> shortest path${res.paths.length > 1 ? "s" : ""} of length <b>${res.distance}</b>.<br>${verdict}</p>` +
      res.paths.map((p, k) => `<div class="result-card" data-k="${k}"><div class="chain">${chainHtml(p.nodes, p.edges)}</div>
        <div style="margin-top:6px">${$("#pathUndirected").checked ? "" : effectBadge(p.sign)}</div></div>`).join("");
    const allNodes = new Set(res.paths.flatMap((p) => p.nodes)), allEdges = res.paths.flatMap((p) => p.edges);
    highlightPath([...allNodes], allEdges);
    out.querySelectorAll(".result-card").forEach((c) => (c.onclick = () => {
      clearCardHighlights(); c.classList.add("on");
      const p = res.paths[+c.dataset.k]; highlightPath(p.nodes, p.edges);
    }));
  }
  $("#findPath").onclick = findPaths;

  // ---------- loops tab ----------
  let loopFilter = "all";
  function renderLoops() {
    const list = loops.filter((l) => loopFilter === "all" || l.sign === +loopFilter);
    $("#loopList").innerHTML =
      `<p class="summary" style="margin-top:0">${loops.length} loops found: <b>${negLoops}</b> negative, <b>${loops.length - negLoops}</b> positive.</p>` +
      list.map((l) => `<li data-k="${loops.indexOf(l)}"><div class="meta"><span>${l.nodes.length}-node loop</span>
        <span class="badge ${l.sign < 0 ? "neg" : "pos"}">${l.sign < 0 ? "negative" : "positive"} feedback</span></div>
        <div class="chain">${chainHtml([...l.nodes, l.nodes[0]], l.edges)}</div></li>`).join("");
    document.querySelectorAll("#loopList li").forEach((li) => (li.onclick = () => {
      clearCardHighlights(); li.classList.add("on");
      const l = loops[+li.dataset.k]; highlightPath(l.nodes, l.edges);
    }));
  }
  segmented("#loopFilter", (v) => { loopFilter = v; renderLoops(); });
  renderLoops();

  // ---------- data tab ----------
  function renderEdges(q = "") {
    q = q.toLowerCase();
    const rows = D.edges.map((e, k) => ({ e, k })).filter(({ e }) => !q || [e.source, e.target, lbl(e.source), lbl(e.target), e.mechanism,
      pathwayById.get(byId.get(e.source).pathway).name].some((s) => s.toLowerCase().includes(q)));
    $("#edgeRows").innerHTML = rows.map(({ e, k }) => `<tr data-k="${k}"><td><b>${esc(lbl(e.source))}</b></td><td>${glyph(e.sign)}</td>
      <td><b>${esc(lbl(e.target))}</b></td><td>${esc(e.mechanism)}${e.indirect ? "*" : ""}</td></tr>`).join("");
    document.querySelectorAll("#edgeRows tr").forEach((tr) => (tr.onclick = () => {
      const e = D.edges[+tr.dataset.k]; highlightPath([e.source, e.target], [+tr.dataset.k]);
    }));
  }
  renderEdges();
  $("#edgeFilter").oninput = (ev) => renderEdges(ev.target.value);
  const download = (name, text, type) => {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([text], { type })); a.download = name; a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  };
  $("#dlJson").onclick = () => download("biograph-network.json", JSON.stringify({
    pathways: D.PATHWAYS, nodes: D.nodes.map((n) => ({ ...n, ...metrics[n.id] })), edges: D.edges, drugs: D.DRUGS,
  }, null, 2), "application/json");
  $("#dlCsv").onclick = () => download("biograph-edges.csv",
    "source,target,sign,mechanism,indirect,note\n" + D.edges.map((e) =>
      [e.source, e.target, e.sign, e.mechanism, e.indirect, `"${e.note.replace(/"/g, '""')}"`].join(",")).join("\n"), "text/csv");

  // ---------- theme, dialog, keyboard ----------
  function refreshTheme() { readPalette(); buildMarkers(); renderLegend(); renderHubs(); paint(); if (state.sim) renderSimResults(); if (state.selected) select(state.selected); }
  $("#themeBtn").onclick = () => {
    const dark = document.documentElement.dataset.theme
      ? document.documentElement.dataset.theme === "dark"
      : matchMedia("(prefers-color-scheme: dark)").matches;
    document.documentElement.dataset.theme = dark ? "light" : "dark";
    try { localStorage.setItem("biograph-theme", document.documentElement.dataset.theme); } catch (_) {}
    refreshTheme();
  };
  try { const t = localStorage.getItem("biograph-theme"); if (t) { document.documentElement.dataset.theme = t; refreshTheme(); } } catch (_) {}
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", refreshTheme);

  $("#aboutBtn").onclick = () => $("#about").showModal();
  document.addEventListener("keydown", (ev) => {
    if (ev.key === "/" && document.activeElement.tagName !== "INPUT") { ev.preventDefault(); $("#search").focus(); }
    if (ev.key === "Escape") { state.highlight = null; clearCardHighlights(); select(null); }
  });
  window.addEventListener("resize", () => fit(0));

  paint();
  fit(0);
  // Gentle settle animation on load
  force.alpha(0.15).restart();
})();
