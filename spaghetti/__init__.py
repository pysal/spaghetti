"""
# `spaghetti` --- Spatial Graphs: Networks, Topology, & Inference
"""

import contextlib
from importlib.metadata import PackageNotFoundError, version

from .network import (
    Network,
    PointPattern,
    SimulatedPointPattern,
    element_as_gdf,
    extract_component,
    regular_lattice,
    spanning_tree,
)

with contextlib.suppress(PackageNotFoundError):
    __version__ = version("spaghetti")
