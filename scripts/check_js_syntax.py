#!/usr/bin/env python3
"""Parse every browser module the site ships, as a module.

`node --check app.js` looks like a syntax gate and is not one: for a file containing
`import`, Node parses it as a *script*, and a duplicate `const` in module scope passes
silently. It did — a redeclared binding reached the page and blanked a whole tab, with
`node --check` reporting nothing. Copying to `.mjs` first makes Node use the module
grammar, which catches it.

The page has no build step, so this is the only thing standing between a typo and a blank
panel. If node is not installed the check says so and does not pretend to have passed.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import DOCS  # noqa: E402

APP = os.path.join(DOCS, "app")


def main() -> int:
    node = shutil.which("node")
    if not node:
        print("node is not installed — JavaScript syntax was NOT checked")
        return 1

    files = sorted(f for f in os.listdir(APP) if f.endswith(".js"))
    if not files:
        print(f"no .js files under {APP}")
        return 1

    print(f"JavaScript syntax — {len(files)} module(s) under docs/app\n")
    failures = 0
    with tempfile.TemporaryDirectory() as tmp:
        for name in files:
            src = os.path.join(APP, name)
            mjs = os.path.join(tmp, name[:-3] + ".mjs")
            shutil.copyfile(src, mjs)
            r = subprocess.run([node, "--check", mjs], capture_output=True, text=True)
            if r.returncode == 0:
                print(f"  ok   {name}")
            else:
                failures += 1
                detail = (r.stderr.strip().splitlines() or ["(no detail)"])
                print(f"  FAIL {name}")
                for line in detail[:6]:
                    print(f"       {line.replace(tmp, 'docs/app')}")

    print()
    if failures:
        print(f"FAILED — {failures} module(s) do not parse")
        return 1
    print("all modules parse as ES modules")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
