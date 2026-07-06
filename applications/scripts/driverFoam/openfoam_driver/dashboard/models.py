from __future__ import annotations

from dataclasses import asdict, dataclass, field


@dataclass
class CaseCard:
    case_id: str
    family: str
    name: str
    goal: str | None = None
    solver: str | None = None
    ionic_model: str | None = None
    status: str = "not_run"            # ran | not_run
    last_run: str | None = None        # ISO8601
    regression: str = "none"           # available | none
    metrics: dict[str, str | float] = field(default_factory=dict)
    plots: list[str] = field(default_factory=list)
    artifacts: dict[str, str] = field(default_factory=dict)
    renders: list[str] = field(default_factory=list)
    is_3d: bool = False
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class CaseView:
    card: CaseCard
    user_notes: str = ""
    tags: list[str] = field(default_factory=list)
    captions: dict[str, str] = field(default_factory=dict)
