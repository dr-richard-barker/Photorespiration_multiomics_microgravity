#!/usr/bin/env python3
"""Fail if any figure carries a legend entry it never draws.

Written after Figure 1b shipped in the manuscript with a five-entry legend over three
lines. It plotted `dishmean_co2`, which is empty for the two vented hardware cases because
a vented enclosure has no closed volume to average, and matplotlib draws an all-NaN series
as nothing at all — silently, with the legend entry intact. An annotation then explained
the two missing lines as "VEGGIE and the open reference both sit on zero", which the
surface data contradicts.

So: an artist with no finite points is an error, and a legend that promises more series
than the axes deliver is an error. Neither is visible in a rendered PNG without counting.

The figure has to be caught on its way out. `style.save` closes it, so inspecting
`plt.get_fignums()` after `main()` returns finds nothing and passes everything — which is
what the first version of this file did, including on the very figure it was written for.
Wrapping `style.save` is what makes the check real.
"""

from __future__ import annotations

import importlib.util
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
sys.path.insert(0, HERE)
sys.path.insert(0, FIGDIR)


def finite_points(artist) -> int:
    """How many points of this artist would actually appear."""
    try:
        if hasattr(artist, "get_xydata"):
            xy = np.asarray(artist.get_xydata(), dtype=float)
        elif hasattr(artist, "get_offsets"):
            xy = np.asarray(artist.get_offsets(), dtype=float)
        else:
            return -1
    except (TypeError, ValueError):
        return -1
    if xy.size == 0:
        return 0
    return int(np.isfinite(xy).all(axis=1).sum())


def inspect(fig, _name: str) -> list[str]:
    problems = []
    for ax in fig.get_axes():
        drawn, empty = 0, []
        for artist in list(ax.get_lines()) + list(ax.collections):
            label = artist.get_label() or ""
            n = finite_points(artist)
            if n == -1:
                continue
            if n == 0:
                # An unlabelled empty artist is usually a placeholder; a labelled one is a
                # promise in the legend that nothing keeps.
                if not label.startswith("_"):
                    empty.append(label or "(unlabelled)")
            else:
                drawn += 1
        for label in empty:
            problems.append(f"'{label}' is labelled for the legend but has no finite points")
        leg = ax.get_legend()
        if leg is not None:
            promised = [t.get_text() for t in leg.get_texts()]
            if len(promised) > drawn + len(ax.patches):
                problems.append(f"legend promises {len(promised)} series "
                                f"({', '.join(promised)}) but the axes draw {drawn}")
    return problems


def check_module(path: str) -> list[str]:
    """Run one figure script, inspecting each figure as it is saved."""
    import style

    name = os.path.basename(path)[:-3]
    found: list[str] = []
    saved = []
    real_save = style.save

    def capturing_save(fig, stem, outdir):
        saved.append(stem)
        found.extend(f"{stem}: {p}" for p in inspect(fig, stem))
        return real_save(fig, stem, outdir)

    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    plt.close("all")
    style.save = capturing_save
    try:
        spec.loader.exec_module(mod)
        mod.main()
    finally:
        style.save = real_save
        plt.close("all")

    if not saved:
        found.append("style.save was never called, so nothing was inspected")
    return found


def main() -> int:
    scripts = sorted(f for f in os.listdir(FIGDIR)
                     if f.startswith("fig") and f.endswith(".py"))
    print(f"Figure completeness — {len(scripts)} figure scripts\n")
    failures: list[str] = []
    for f in scripts:
        found = check_module(os.path.join(FIGDIR, f))
        print(f"  {'FAIL' if found else 'ok  '} {f}")
        for p in found:
            print(f"       {p}")
        failures += found

    print()
    if failures:
        print(f"FAILED — {len(failures)} problem(s)")
        return 1
    print("every labelled series in every figure has data behind it")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
