"""Run the browser-side test suite from pytest.

The drawing lives in ``src/plastr/web/js`` and is most of what a user actually
sees, but it used to have no tests at all: a regression there was only ever
caught by someone looking at a picture. These run under Node's built-in test
runner, so there is no package.json, no dependency to install and no build
step -- the rule the web/ directory is written to.

Node is optional. Without it these skip, and the Python suite is unaffected.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

JS_TESTS = Path(__file__).resolve().parent / "js"


def _node() -> str | None:
    return shutil.which("node")


requires_node = pytest.mark.skipif(
    _node() is None,
    reason="node is not installed; the browser-side suite needs it",
)


def _cases() -> list[Path]:
    return sorted(JS_TESTS.glob("*.test.mjs"))


@requires_node
@pytest.mark.parametrize("case", _cases(), ids=lambda p: p.stem)
def test_web_js(case: Path) -> None:
    node = _node()
    assert node is not None
    proc = subprocess.run(
        [node, "--test", str(case)],
        cwd=str(case.parent),
        capture_output=True,
        text=True,
        timeout=300,
    )
    if proc.returncode != 0:
        # Node's TAP output carries the failing assertion; showing it here saves
        # re-running the command by hand to find out what broke.
        pytest.fail(
            f"{case.name} failed\n\n{proc.stdout[-8000:]}\n{proc.stderr[-4000:]}",
            pytrace=False,
        )


def test_js_suite_is_discovered() -> None:
    """Guard the glob: a renamed directory would silently run nothing."""
    assert _cases(), f"no *.test.mjs found under {JS_TESTS}"
