from __future__ import annotations

import re
import sqlite3
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path


@dataclass
class Annotation:
    case_id: str
    notes: str = ""
    tags: list[str] = field(default_factory=list)


def _norm_tag(t: str) -> str:
    t = t.strip().lower()
    t = re.sub(r"\s+", "-", t)
    return re.sub(r"[^a-z0-9._-]", "", t)


def _norm_tags(tags: list[str]) -> list[str]:
    seen: list[str] = []
    for t in tags:
        n = _norm_tag(t)
        if n and n not in seen:
            seen.append(n)
    return seen


def _now() -> str:
    return datetime.now(tz=timezone.utc).isoformat()


class AnnotationStore:
    def __init__(self, db_path: Path):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_schema()

    def _conn(self) -> sqlite3.Connection:
        return sqlite3.connect(self.db_path)

    def _init_schema(self) -> None:
        with self._conn() as c:
            c.execute(
                "CREATE TABLE IF NOT EXISTS case_annotation ("
                "case_id TEXT PRIMARY KEY, notes TEXT DEFAULT '', "
                "tags TEXT DEFAULT '', updated_at TEXT)"
            )
            c.execute(
                "CREATE TABLE IF NOT EXISTS figure_caption ("
                "case_id TEXT, figure_id TEXT, caption TEXT DEFAULT '', "
                "updated_at TEXT, PRIMARY KEY (case_id, figure_id))"
            )

    def upsert_annotation(self, case_id: str, notes: str, tags: list[str]) -> None:
        tag_str = ",".join(_norm_tags(tags))
        with self._conn() as c:
            c.execute(
                "INSERT INTO case_annotation (case_id, notes, tags, updated_at) "
                "VALUES (?,?,?,?) ON CONFLICT(case_id) DO UPDATE SET "
                "notes=excluded.notes, tags=excluded.tags, updated_at=excluded.updated_at",
                (case_id, notes, tag_str, _now()),
            )

    def get_annotation(self, case_id: str) -> Annotation:
        with self._conn() as c:
            row = c.execute(
                "SELECT notes, tags FROM case_annotation WHERE case_id=?", (case_id,)
            ).fetchone()
        if not row:
            return Annotation(case_id=case_id)
        notes, tags = row
        tag_list = [t for t in (tags or "").split(",") if t]
        return Annotation(case_id=case_id, notes=notes or "", tags=tag_list)

    def upsert_caption(self, case_id: str, figure_id: str, caption: str) -> None:
        with self._conn() as c:
            c.execute(
                "INSERT INTO figure_caption (case_id, figure_id, caption, updated_at) "
                "VALUES (?,?,?,?) ON CONFLICT(case_id, figure_id) DO UPDATE SET "
                "caption=excluded.caption, updated_at=excluded.updated_at",
                (case_id, figure_id, caption, _now()),
            )

    def get_captions(self, case_id: str) -> dict[str, str]:
        with self._conn() as c:
            rows = c.execute(
                "SELECT figure_id, caption FROM figure_caption WHERE case_id=?", (case_id,)
            ).fetchall()
        return {fid: cap for fid, cap in rows}

    def all_tags(self) -> list[str]:
        with self._conn() as c:
            rows = c.execute("SELECT tags FROM case_annotation").fetchall()
        out: set[str] = set()
        for (tags,) in rows:
            out.update(t for t in (tags or "").split(",") if t)
        return sorted(out)
