from dataclasses import dataclass

from openfoam_driver.dict_entries import Phase


@dataclass(frozen=True)
class ValidationError:
    phase: Phase
    field: str
    message: str
    level: str  # "error" | "warning"
