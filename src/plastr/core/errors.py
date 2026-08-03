"""Typed errors used across Plastr.

Parse errors carry enough context (file, line number, offending text) that a user
staring at a truncated GFA can find the problem without reading a traceback.
"""

from __future__ import annotations


class PlastrError(Exception):
    """Base class for every error Plastr raises deliberately."""


class PlastrFormatError(PlastrError):
    """A file could not be parsed."""

    def __init__(
        self,
        message: str,
        path: str | None = None,
        line_no: int | None = None,
        line: str | None = None,
    ) -> None:
        self.message = message
        self.path = path
        self.line_no = line_no
        self.line = line
        super().__init__(str(self))

    def __str__(self) -> str:  # pragma: no cover - formatting only
        where = ""
        if self.path:
            where = f" in {self.path}"
            if self.line_no is not None:
                where += f" at line {self.line_no}"
        text = f"{self.message}{where}"
        if self.line:
            snippet = self.line.rstrip("\n")
            if len(snippet) > 200:
                snippet = snippet[:200] + "..."
            text += f"\n  > {snippet}"
        return text


class MissingDependencyError(PlastrError):
    """An optional third-party tool or module is needed but not installed."""

    def __init__(self, name: str, hint: str = "") -> None:
        self.name = name
        self.hint = hint
        msg = f"'{name}' is required for this feature but is not available."
        if hint:
            msg += f" {hint}"
        super().__init__(msg)


class GraphOperationError(PlastrError):
    """An edit to the graph could not be applied."""
