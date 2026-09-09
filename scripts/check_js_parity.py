#!/usr/bin/env python3
"""Assert that docs/app/fvcb.js reproduces scripts/fvcb.py to within tolerance.

The website lets visitors drive the model themselves. If the JavaScript drifts from the
Python that produced the figures and the manuscript numbers, the site quietly starts telling
a different story from the paper. This runs both over the same parameter grid and fails on
any disagreement.

Needs node on PATH.
"""
from __future__ import annotations
import json, os, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fvcb import LeafParams, solve_operating_point  # noqa: E402
from paths import DOCS, REPO  # noqa: E402

TOL = 1e-6
GRID = [
    dict(g_bl=g, Ca=ca, TleafC=t, Q=q)
    for g in (1.0, 0.54, 0.291, 0.109)
    for ca in (400.0, 250.0, 150.0, 100.0, 70.0)
    for t in (20.0, 25.0, 30.0)
    for q in (400.0, 1000.0, 1500.0)
]
FIELDS = ["Cc", "A", "phi", "Vo_over_Vc", "Vc", "Vo", "Rp", "Gamma_star", "Ci"]


def run_js() -> list[dict]:
    script = f"""
import {{ solveOperatingPoint, DEFAULT_PARAMS }} from '{os.path.join(DOCS, "app", "fvcb.js")}';
const grid = {json.dumps(GRID)};
const out = grid.map(g => solveOperatingPoint(g.g_bl, DEFAULT_PARAMS,
    {{ Ca: g.Ca, TleafC: g.TleafC, Q: g.Q }}));
console.log(JSON.stringify(out));
"""
    with tempfile.NamedTemporaryFile("w", suffix=".mjs", delete=False,
                                     dir=REPO) as fh:
        fh.write(script)
        path = fh.name
    try:
        res = subprocess.run(["node", path], capture_output=True, text=True, timeout=120)
        if res.returncode != 0:
            raise SystemExit(f"node failed:\n{res.stderr[:800]}")
        return json.loads(res.stdout)
    finally:
        os.unlink(path)


def main() -> int:
    print(f"Parity check: docs/app/fvcb.js vs scripts/fvcb.py over {len(GRID)} points")
    js = run_js()
    p = LeafParams()
    worst = 0.0
    worst_where = None
    failures = 0

    for g, j in zip(GRID, js):
        py = solve_operating_point(g["g_bl"], p, Ca=g["Ca"], Tleaf_C=g["TleafC"], Q=g["Q"])
        for f in FIELDS:
            a, b = py[f], j[f]
            if a != a and b != b:      # both NaN
                continue
            denom = max(1.0, abs(a))
            d = abs(a - b) / denom
            if d > worst:
                worst, worst_where = d, (f, g)
            if d > TOL:
                failures += 1
                if failures <= 5:
                    print(f"  MISMATCH {f} at {g}: python {a!r} vs js {b!r} (rel {d:.2e})")

    print(f"  worst relative difference {worst:.3e} on {worst_where[0] if worst_where else '-'}")
    if failures:
        print(f"FAILED — {failures} value(s) exceed {TOL}")
        return 1
    print(f"PASS — every value agrees to better than {TOL}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
