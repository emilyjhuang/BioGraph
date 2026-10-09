/*
 * BioGraph analysis engine: graph algorithms over a signed, directed network.
 * Pure functions with no DOM access, so they run in the browser and in Node tests.
 */
(function (root) {
  /** Build adjacency structures from {nodes, edges}. Edge endpoints are node ids. */
  function buildGraph(nodes, edges) {
    const index = new Map(nodes.map((n, i) => [n.id, i]));
    const n = nodes.length;
    const out = Array.from({ length: n }, () => []);
    const inc = Array.from({ length: n }, () => []);
    edges.forEach((e, k) => {
      const s = index.get(typeof e.source === "object" ? e.source.id : e.source);
      const t = index.get(typeof e.target === "object" ? e.target.id : e.target);
      if (s === undefined || t === undefined) throw new Error(`Unknown node in edge ${k}`);
      out[s].push({ to: t, sign: e.sign, edge: k });
      inc[t].push({ from: s, sign: e.sign, edge: k });
    });
    return { nodes, edges, index, n, out, inc };
  }

  /** Brandes' algorithm for betweenness centrality on an unweighted directed graph (normalised to [0,1]). */
  function betweenness(g) {
    const { n, out } = g;
    const cb = new Float64Array(n);
    for (let s = 0; s < n; s++) {
      const stack = [];
      const pred = Array.from({ length: n }, () => []);
      const sigma = new Float64Array(n); sigma[s] = 1;
      const dist = new Int32Array(n).fill(-1); dist[s] = 0;
      const queue = [s];
      for (let qi = 0; qi < queue.length; qi++) {
        const v = queue[qi];
        stack.push(v);
        for (const { to: w } of out[v]) {
          if (dist[w] < 0) { dist[w] = dist[v] + 1; queue.push(w); }
          if (dist[w] === dist[v] + 1) { sigma[w] += sigma[v]; pred[w].push(v); }
        }
      }
      const delta = new Float64Array(n);
      while (stack.length) {
        const w = stack.pop();
        for (const v of pred[w]) delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w]);
        if (w !== s) cb[w] += delta[w];
      }
    }
    const norm = n > 2 ? 1 / ((n - 1) * (n - 2)) : 1;
    return Array.from(cb, (v) => v * norm);
  }

  /** PageRank on the directed graph (dangling mass redistributed uniformly). */
  function pagerank(g, damping = 0.85, iters = 100, tol = 1e-10) {
    const { n, out } = g;
    let pr = new Float64Array(n).fill(1 / n);
    for (let it = 0; it < iters; it++) {
      const next = new Float64Array(n).fill((1 - damping) / n);
      let dangling = 0;
      for (let v = 0; v < n; v++) {
        if (out[v].length === 0) { dangling += pr[v]; continue; }
        const share = (damping * pr[v]) / out[v].length;
        for (const { to } of out[v]) next[to] += share;
      }
      for (let v = 0; v < n; v++) next[v] += (damping * dangling) / n;
      let diff = 0;
      for (let v = 0; v < n; v++) diff += Math.abs(next[v] - pr[v]);
      pr = next;
      if (diff < tol) break;
    }
    return Array.from(pr);
  }

  /** Degree, betweenness and PageRank for every node, keyed by node id. */
  function centralities(g) {
    const bc = betweenness(g);
    const pr = pagerank(g);
    const result = {};
    g.nodes.forEach((node, i) => {
      result[node.id] = {
        inDegree: g.inc[i].length,
        outDegree: g.out[i].length,
        degree: g.inc[i].length + g.out[i].length,
        betweenness: bc[i],
        pagerank: pr[i],
      };
    });
    return result;
  }

  /**
   * Signed network propagation of a perturbation.
   *
   * Perturbed nodes are clamped to a value in [-1, 1] (-1 = inhibited/knocked out,
   * +1 = activated). Every other node relaxes towards
   *     x_v = decay * tanh( gain * sum_{u->v} sign(u,v) * x_u / sqrt(indeg(v)) )
   * with damped synchronous updates. Values are changes relative to the
   * unperturbed baseline: a qualitative estimate of direction and relative
   * magnitude, not a kinetic simulation.
   *
   * Returns { values: {id: x}, frames: [Float64Array], steps, converged }.
   */
  function propagate(g, perturb, opts = {}) {
    const { gain = 1.2, decay = 0.85, damping = 0.5, maxSteps = 200, tol = 1e-5, removed = [] } = opts;
    const { n, inc, index } = g;
    const clamp = new Map();
    for (const [id, v] of Object.entries(perturb)) {
      if (!index.has(id)) throw new Error(`Unknown node: ${id}`);
      clamp.set(index.get(id), Math.max(-1, Math.min(1, v)));
    }
    const gone = new Set(removed.map((id) => index.get(id)));
    let x = new Float64Array(n);
    for (const [i, v] of clamp) x[i] = v;
    const frames = [Float64Array.from(x)];
    let converged = false, step = 0;
    for (step = 1; step <= maxSteps; step++) {
      const next = new Float64Array(n);
      let diff = 0;
      for (let v = 0; v < n; v++) {
        if (clamp.has(v)) { next[v] = clamp.get(v); continue; }
        if (gone.has(v)) continue;
        const ins = inc[v].filter((e) => !gone.has(e.from));
        if (ins.length === 0) continue;
        let s = 0;
        for (const e of ins) s += e.sign * x[e.from];
        const target = decay * Math.tanh((gain * s) / Math.sqrt(ins.length));
        next[v] = (1 - damping) * x[v] + damping * target;
        diff = Math.max(diff, Math.abs(next[v] - x[v]));
      }
      x = next;
      frames.push(Float64Array.from(x));
      if (diff < tol) { converged = true; break; }
    }
    const values = {};
    g.nodes.forEach((node, i) => { values[node.id] = x[i]; });
    return { values, frames, steps: step, converged };
  }

  /**
   * All shortest directed paths from `src` to `dst` (capped), each with its net sign
   * (product of edge signs). Returns { distance, paths: [{nodes:[ids], edges:[idx], sign}] }.
   */
  function shortestPaths(g, src, dst, { cap = 25, undirected = false } = {}) {
    const s = g.index.get(src), t = g.index.get(dst);
    if (s === undefined || t === undefined) throw new Error("Unknown node");
    const nbrs = (v) => {
      const fwd = g.out[v].map((e) => ({ w: e.to, edge: e.edge, sign: e.sign }));
      if (!undirected) return fwd;
      return fwd.concat(g.inc[v].map((e) => ({ w: e.from, edge: e.edge, sign: e.sign })));
    };
    const dist = new Int32Array(g.n).fill(-1); dist[s] = 0;
    const pred = Array.from({ length: g.n }, () => []);
    const queue = [s];
    for (let qi = 0; qi < queue.length; qi++) {
      const v = queue[qi];
      if (v === t) continue;
      for (const { w, edge, sign } of nbrs(v)) {
        if (dist[w] < 0) { dist[w] = dist[v] + 1; queue.push(w); }
        if (dist[w] === dist[v] + 1) pred[w].push({ v, edge, sign });
      }
    }
    if (dist[t] < 0) return { distance: Infinity, paths: [] };
    const paths = [];
    (function walk(v, nodesAcc, edgesAcc, sign) {
      if (paths.length >= cap) return;
      if (v === s) {
        paths.push({
          nodes: [s, ...nodesAcc].map((i) => g.nodes[i].id),
          edges: edgesAcc.slice(),
          sign,
        });
        return;
      }
      for (const p of pred[v]) walk(p.v, [v, ...nodesAcc], [p.edge, ...edgesAcc], sign * p.sign);
    })(t, [], [], 1);
    return { distance: dist[t], paths };
  }

  /**
   * Enumerate elementary directed cycles up to `maxLen` (feedback loops).
   * Each cycle is reported once, rooted at its lowest-index node.
   * sign = product of edge signs: +1 positive feedback, -1 negative feedback.
   */
  function feedbackLoops(g, maxLen = 8, cap = 500) {
    const loops = [];
    const onPath = new Uint8Array(g.n);
    for (let s = 0; s < g.n && loops.length < cap; s++) {
      const nodePath = [s], edgePath = [];
      onPath[s] = 1;
      (function dfs(v, sign) {
        if (loops.length >= cap) return;
        for (const { to, edge, sign: es } of g.out[v]) {
          if (to < s) continue;
          if (to === s) {
            loops.push({
              nodes: nodePath.map((i) => g.nodes[i].id),
              edges: [...edgePath, edge],
              sign: sign * es,
            });
            continue;
          }
          if (onPath[to] || nodePath.length >= maxLen) continue;
          onPath[to] = 1; nodePath.push(to); edgePath.push(edge);
          dfs(to, sign * es);
          onPath[to] = 0; nodePath.pop(); edgePath.pop();
        }
      })(s, 1);
      onPath[s] = 0;
    }
    loops.sort((a, b) => a.nodes.length - b.nodes.length);
    return loops;
  }

  /** Nodes reachable downstream of `id` (directed). */
  function downstream(g, id) {
    const s = g.index.get(id);
    const seen = new Uint8Array(g.n); seen[s] = 1;
    const queue = [s];
    for (let qi = 0; qi < queue.length; qi++)
      for (const { to } of g.out[queue[qi]]) if (!seen[to]) { seen[to] = 1; queue.push(to); }
    return queue.slice(1).map((i) => g.nodes[i].id);
  }

  const api = { buildGraph, betweenness, pagerank, centralities, propagate, shortestPaths, feedbackLoops, downstream };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.BioGraphAnalysis = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
