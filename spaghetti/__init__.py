"""
# `spaghetti` --- Spatial Graphs: Networks, Topology, & Inference
"""

from .network import Network, PointPattern, SimulatedPointPattern
from .network import extract_component, spanning_tree
from .network import element_as_gdf, regular_lattice

from . import _version

__version__ = _version.get_versions()["version"]
