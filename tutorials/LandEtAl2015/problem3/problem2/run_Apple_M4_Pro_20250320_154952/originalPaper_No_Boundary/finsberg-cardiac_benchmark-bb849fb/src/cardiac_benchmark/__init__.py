try:
    # importlib.metadata is present in Python 3.8 and later
    import importlib.metadata as importlib_metadata
except ImportError:
    # use the shim package importlib-metadata pre-3.8
    import importlib_metadata as importlib_metadata  # type: ignore

from . import benchmark1
from . import benchmark2
from . import utils
from . import cli
from . import geometry
from . import material
from . import microstructure
from . import postprocess
from . import pressure_model
from . import activation_model
from . import problem
from . import solver

meta = importlib_metadata.metadata("cardiac-benchmark")
__version__ = meta["Version"]
__author__ = meta["Author"]
__license__ = meta["License"]
__email__ = meta["Author-email"]
__program_name__ = meta["Name"]


__all__ = [
    "benchmark1",
    "benchmark2",
    "utils",
    "cli",
    "geometry",
    "material",
    "microstructure",
    "postprocess",
    "pressure_model",
    "problem",
    "solver",
    "activation_model",
]
