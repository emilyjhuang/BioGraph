const test = require("node:test");
const assert = require("node:assert/strict");
const A = require("../js/analysis.js");
const D = require("../js/data.js");

const g = A.buildGraph(D.nodes, D.edges);

test("data integrity: unique ids, valid pathways, no duplicate edges", () => {
  const ids = new Set(D.nodes.map((n) => n.id));
  assert.equal(ids.size, D.nodes.length);
  const pw = new Set(D.PATHWAYS.map((p) => p.id));
  for (const n of D.nodes) assert.ok(pw.has(n.pathway), n.id);
  const keys = D.edges.map((e) => `${e.source}>${e.target}`);
  assert.equal(new Set(keys).size, keys.length);
  for (const d of D.DRUGS) for (const t of d.targets) assert.ok(ids.has(t), `${d.id} -> ${t}`);
  for (const s of D.SCENARIOS) for (const id of [...Object.keys(s.perturb), ...s.watch]) assert.ok(ids.has(id), id);
});

test("betweenness and pagerank on a toy chain", () => {
  const toy = A.buildGraph(
    [{ id: "a" }, { id: "b" }, { id: "c" }],
    [{ source: "a", target: "b", sign: 1 }, { source: "b", target: "c", sign: 1 }]
  );
  const bc = A.betweenness(toy);
  assert.deepEqual(bc.map((v) => +v.toFixed(3)), [0, 0.5, 0]);
  const pr = A.pagerank(toy);
  assert.ok(Math.abs(pr.reduce((a, b) => a + b) - 1) < 1e-9);
  assert.ok(pr[2] > pr[1] && pr[1] > pr[0]);
});

test("shortest paths report net sign", () => {
  const r = A.shortestPaths(g, "MDM2", "BAX");
  assert.equal(r.distance, 2);
  assert.ok(r.paths.every((p) => p.sign === -1));
  assert.equal(A.shortestPaths(g, "CASP3", "EGF").paths.length, 0);
});

test("feedback loops include the classic ones with the right sign", () => {
  const loops = A.feedbackLoops(g);
  const find = (a, b) => loops.find((l) => l.nodes.length === 2 && l.nodes.includes(a) && l.nodes.includes(b));
  assert.equal(find("TP53", "MDM2").sign, -1);
  assert.equal(find("RELA", "NFKBIA").sign, -1);
  assert.equal(find("MAPK1", "DUSP6").sign, -1);
  const s6k = loops.find((l) => ["IRS1", "PIK3CA", "PDPK1", "RPS6KB1"].every((x) => l.nodes.includes(x)) && l.nodes.length === 4);
  assert.equal(s6k.sign, -1);
});

test("propagation converges and preserves clamps", () => {
  const r = A.propagate(g, { BRAF: -1 });
  assert.ok(r.converged);
  assert.equal(r.values.BRAF, -1);
  for (const v of Object.values(r.values)) assert.ok(v >= -1 && v <= 1);
});

// Each guided scenario must reproduce the qualitative biology described in its text.
const expectations = {
  braf: { MAPK1: -1, FOS: -1, SOS1: 1, KRAS: 1 },
  rapalog: { RPS6KB1: -1, IRS1: 1, AKT1: 1, MTORC2: 1 },
  mdm2: { TP53: 1, CDKN1A: 1, BBC3: 1, BAX: 1 },
  dsb: { TP53: 1, CDC25A: -1, CDK2: -1, CDKN1A: 1 },
  pten: { AKT1: 1, FOXO3: -1, MTORC1: 1, BCL2L11: -1 },
  apc: { CTNNB1: 1, MYC: 1, CCND1: 1, RB1: -1 },
  combo: { MAPK1: -1, MAPK3: -1, ELK1: -1, BCL2L11: 1 },
};

for (const sc of D.SCENARIOS) {
  test(`scenario "${sc.title}" gives the expected direction of change`, () => {
    const { values, converged } = A.propagate(g, sc.perturb);
    assert.ok(converged);
    for (const [id, dir] of Object.entries(expectations[sc.id])) {
      assert.ok(Math.sign(values[id]) === dir && Math.abs(values[id]) > 0.01, `${id}: ${values[id].toFixed(3)}`);
    }
  });
}
