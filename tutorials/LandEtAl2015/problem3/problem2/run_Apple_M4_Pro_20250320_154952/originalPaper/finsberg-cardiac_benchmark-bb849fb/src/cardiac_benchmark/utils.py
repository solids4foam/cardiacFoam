import json
import logging
from pathlib import Path
from typing import Any
from typing import Dict
from typing import Optional

import dolfin

logger = logging.getLogger(__name__)


class ConstantEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, dolfin.Constant):
            return float(obj)
        if isinstance(obj, Path):
            return obj.as_posix()
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


def _update_parameters(
    _par: Dict[str, Any],
    par: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    if par is None:
        par = {}
    for key, value in par.items():
        if key not in _par:
            logger.warning(f"Invalid key {key}")
            continue

        if isinstance(_par[key], dolfin.Constant):
            _par[key].assign(value)
        else:
            _par[key] = value
    return _par
