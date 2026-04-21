"""table_writer.py — Standard tabular output with metadata envelope."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


@dataclass
class TableMetadata:
    """Metadata envelope attached to every tutorial table output.

    Attributes
    ----------
    tutorial:
        Human-readable tutorial name, e.g. ``"NiedererEtAl2012"``.
    units:
        Mapping of column name → unit string, e.g.
        ``{"activationTime": "ms", "DX": "mm"}``.
    generated_at:
        UTC ISO-8601 timestamp string.  Auto-filled on construction if empty.
    """

    tutorial: str
    units: dict[str, str] = field(default_factory=dict)
    generated_at: str = ""

    def __post_init__(self) -> None:
        if not self.generated_at:
            self.generated_at = datetime.now(timezone.utc).isoformat()


class TableWriter:
    """Write tutorial summary tables as CSV (with comment envelope) and HTML."""

    @staticmethod
    def write(
        rows: list[dict[str, Any]],
        output_dir: str | Path,
        filename_stem: str,
        label: str,
        metadata: TableMetadata,
    ) -> list[dict[str, Any]]:
        """Write *rows* as ``<filename_stem>.csv`` and ``<filename_stem>.html``.

        Parameters
        ----------
        rows:
            List of dicts where every dict has the same keys (tutorial-defined
            columns).  An empty list is allowed — only the envelope is written.
        output_dir:
            Directory into which the files are written.
        filename_stem:
            Base name without extension, e.g. ``"NiedererEtAl2012_summary"``.
        label:
            Human-readable description used in the artifact entry.
        metadata:
            :class:`TableMetadata` providing tutorial name, units, and timestamp.

        Returns
        -------
        Two artifact dicts (CSV and HTML), each with keys ``path``, ``label``,
        ``kind``, ``format``.  Paths are relative filenames (not absolute).
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        csv_path = output_dir / f"{filename_stem}.csv"
        html_path = output_dir / f"{filename_stem}.html"

        fieldnames: list[str] = list(rows[0].keys()) if rows else []

        # ------------------------------------------------------------------
        # CSV: comment-line envelope then data
        # ------------------------------------------------------------------
        lines: list[str] = [
            f"# tutorial: {metadata.tutorial}",
            f"# generated_at: {metadata.generated_at}",
            f"# units: {json.dumps(metadata.units)}",
        ]
        if fieldnames:
            lines.append(",".join(fieldnames))
            for row in rows:
                lines.append(",".join(str(row.get(f, "")) for f in fieldnames))
        csv_path.write_text("\n".join(lines) + "\n")

        # ------------------------------------------------------------------
        # HTML: styled table with metadata header block
        # ------------------------------------------------------------------
        parts: list[str] = [
            "<!DOCTYPE html><html><head>",
            "<style>",
            "body{font-family:Arial,sans-serif;margin:20px}",
            "table{border-collapse:collapse;width:100%}",
            "th,td{border:1px solid #ccc;padding:6px 10px;text-align:left}",
            "th{background:#f0f0f0}",
            ".meta{color:#555;margin-bottom:14px;font-size:13px;line-height:1.6}",
            "</style></head><body>",
            "<div class='meta'>",
            f"<strong>tutorial:</strong> {metadata.tutorial}&nbsp;&nbsp;",
            f"<strong>generated_at:</strong> {metadata.generated_at}<br>",
            f"<strong>units:</strong> {json.dumps(metadata.units)}",
            "</div>",
            "<table><thead><tr>",
        ]
        for f in fieldnames:
            parts.append(f"<th>{f}</th>")
        parts.append("</tr></thead><tbody>")
        for row in rows:
            parts.append("<tr>")
            for f in fieldnames:
                parts.append(f"<td>{row.get(f, '')}</td>")
            parts.append("</tr>")
        parts.append("</tbody></table></body></html>")
        html_path.write_text("".join(parts))

        return [
            {
                "path": csv_path.name,
                "label": label,
                "kind": "table",
                "format": "csv",
            },
            {
                "path": html_path.name,
                "label": label,
                "kind": "table",
                "format": "html",
            },
        ]
