"""
# `spaghetti` --- Spatial Graphs: Networks, Topology, & Inference
"""

from . import _version
from .network import (
    Network,
    PointPattern,
    SimulatedPointPattern,
    element_as_gdf,
    extract_component,
    regular_lattice,
    spanning_tree,
)

__version__ = _version.get_versions()["version"]
