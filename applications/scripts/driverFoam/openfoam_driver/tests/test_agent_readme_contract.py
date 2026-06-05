from __future__ import annotations

import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[5]


LEGACY_AGENT_TOKENS = (
    "electroModel MonoDomainSolver",
    "MonoDomainSolverCoeffs",
    "PseudoECGDomain",
    "PseudoECGDomainCoeffs",
    "ecgModel PseudoECGDomain",
    "foamctl sim --tutorial ECG",
    "tutorials/ECG/",
)


def _agent_facing_readmes() -> list[Path]:
    candidates = [
        *sorted((REPO_ROOT / "tutorials" / "PATHOS").glob("*/README.md")),
        REPO_ROOT / "tutorials" / "manufacturedSolutions" / "bidomain" / "README.md",
        REPO_ROOT / "tutorials" / "manufacturedSolutions" / "bathBidomain" / "README.md",
        REPO_ROOT / "tutorials" / "manufacturedSolutions" / "monodomainPseudoECG" / "README.md",
    ]
    tracked = set(
        subprocess.run(
            ["git", "ls-files", "--", *(str(path.relative_to(REPO_ROOT)) for path in candidates)],
            cwd=REPO_ROOT,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        ).stdout.splitlines()
    )
    return [path for path in candidates if str(path.relative_to(REPO_ROOT)) in tracked]


def test_agent_facing_readmes_do_not_teach_legacy_selector_contracts() -> None:
    violations: list[str] = []
    for readme in _agent_facing_readmes():
        if not readme.exists():
            continue
        text = readme.read_text()
        for token in LEGACY_AGENT_TOKENS:
            if token in text:
                violations.append(f"{readme.relative_to(REPO_ROOT)}: {token}")

    assert violations == []


def test_documented_regression_scripts_exist() -> None:
    violations: list[str] = []
    for readme in _agent_facing_readmes():
        if not readme.exists():
            continue
        text = readme.read_text()
        if "./runRegressionTest.sh" in text and not (
            readme.parent / "runRegressionTest.sh"
        ).exists():
            violations.append(str(readme.relative_to(REPO_ROOT)))
        if "./regressionTest.sh" in text and not (
            readme.parent / "regressionTest.sh"
        ).exists():
            violations.append(str(readme.relative_to(REPO_ROOT)))

    assert violations == []
