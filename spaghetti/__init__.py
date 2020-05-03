__version__ = "1.5.0.rc5"
"""
:mod:`spaghetti` --- Spatial Graphs: Networks, Topology, & Inference
====================================================================

"""
from .network import Network, PointPattern, SimulatedPointPattern
from .network import extract_component, spanning_tree
from .network import element_as_gdf, regular_lattice
