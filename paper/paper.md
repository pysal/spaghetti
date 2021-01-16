---
title: 'spaghetti: spatial network analysis in PySAL'
tags:
  - Python
  - PySAL
  - spatial networks
  - network analysis
authors:
  - name: James D. Gaboardi
    orcid: 0000-0002-4776-6826
    affiliation: 1
  - name: Sergio Rey
    orcid: 0000-0001-5857-9762
    affiliation: 2
  - name: Stefanie Lumnitz
    orcid: 0000-0002-7007-5812
    affiliation: 3
affiliations:
 - name: Pennsylvania State University
   index: 1
 - name: Center for Geospatial Sciences, University of California Riverside
   index: 2
 - name: Directorate of Earth Observation Programs, ESRIN, European Space Agency
   index: 3
date: 30 October 2020
bibliography: paper.bib
---

# Summary

The role spatial networks, such as streets, play on the human experience cannot be overstated. All of our daily activities fall along, or in close proximity to, roads, bike paths, and subway systems to name a few. Therefore, when performing spatial analysis in many cases considering network space, as opposed to Euclidean space, allows for a more precise representation of daily human action and movement patterns. For example, people generally cannot get to work by driving in a straight line directly from their home, but move along paths within networks. To this end, `spaghetti` (**spa**tial **g**rap**h**s: n**et**works, **t**opology, & **i**nference), a sub-module in the wider PySAL ecosystem, was developed to address network-centric research questions with a strong focus on spatial analysis [@pysal2007;@pysal2015;@gaboardi2018].

# Statement of Need

The most well-known network analysis package within the Python scientific stack is [NetworkX](https://networkx.github.io) [@hagberg2008], which can be used for modelling any type of complex network (e.g. social, spatial, etc.). [OSMnx](https://osmnx.readthedocs.io/en/stable/) [@boeing2017] is built on top of NetworkX and queries [OpenStreetMap](https://openstreetmap.org) for modelling street networks with resultant network objects returned within a `geopandas.GeoDataFrame` [@geopandas2020]. Another package, [`pandana`](https://github.com/UDST/pandana) [@foti2012generalized], is built on top of `pandas` [@mckinney-proc-scipy-2010;@reback2020pandas] with a focus on shortest path calculation and accessibility measures. Within the realm of Python, the functionality provided by  [`snkit`](https://github.com/tomalrussell/snkit) [@russell2019] is most comparable to `spaghetti`, though it's main purpose is the processing of raw line data into clean network objects. Outside of Python, [SANET](http://sanet.csis.u-tokyo.ac.jp) [@okabe2006a] is the most closely related project to `spaghetti`, however, it is not written in Python and provides a GUI plugin for GIS software such as QGIS. Moreover, SANET is not fully open source.

# Current Functionality

Considering the related projects in the [Statement of Need](#statement-of-need) detailed above, `spaghetti` fills a niche for not only the processing of spatial network objects, but also post-processing analysis. In other words, this package can be used to study the network *itself* or provide the foundation for studying network-based phenomena, such as crimes along city streets, all within a fully open-source environment. Considering this, the primary purpose of `spaghetti` is creating network objects: collections of vertices and arcs, and their topological relationships. The creation of a network object is realized through the following general steps:

 1. read in line data or create features (regular lattices)
 1. generate the network representation 
 1. extract contiguity weights (if desired)
 1. identify connected components (if desired)
 1. extract graph representation of the network (if desired)

After the creation of a base network object it can be manipulated, analyzed, and utilized as the input for subsequent modelling scenarios. The following are several such examples:

 * allocating observation point patterns to the network
 * calculating all neighbor distance matrices
    * point type A to point type A (auto)
    * point type A to point type B (cross)
 * simulating point patterns that can be used within the  [*K* function](https://pysal.org/spaghetti/generated/spaghetti.Network.html#spaghetti.Network.GlobalAutoK) for cluster analysis [osullivan_unwin2010.ch5;@okabe2012]
 * splitting the network into (nearly) uniform segments
 * extracting features as `geopandas.GeoDataFrame` objects:
    * network arcs, vertices and point patterns 
    * largest/longest components
    * shortest paths
    * minimum/maximum spanning trees

The following demonstrates several functionalities mentioned above, including feature creation, network instantiation, and feature extraction, along with a supplementary plot in \autoref{fig:gridweights}.

```python
import spaghetti
%matplotlib inline
# generate network
lattice = spaghetti.regular_lattice((0,0,3,3), 2, exterior=True)
ntw = spaghetti.Network(in_data=lattice)
# extract network elements
vertices_df, arcs_df = spaghetti.element_as_gdf(ntw, vertices=True, arcs=True)
# plot
base_kws = {"figsize":(12, 12), "lw":5, "color":"k", "zorder":0}
base = arcs_df.plot(**base_kws, alpha=.35)
node_kws, edge_kws = {"s":100, "zorder":2}, {"zorder":1}
w_kws = {"edge_kws":edge_kws, "node_kws":node_kws}
ntw.w_network.plot(arcs_df, indexed_on="id", ax=base, **w_kws)
vertices_df.plot(ax=base, fc="r", ec="k", markersize=50, zorder=2)
```

![A 4x4 regular lattice with network arcs in gray and vertices in red. Connectivity is demonstrated with `libpysal` spatial weights, which are plotted over the network in black [@libpysal2020].\label{fig:gridweights}](figs/spaghetti_network.png)

The overview presented here provides a high-level summary of functionality. More detailed examples and applications can be found in the *Tutorials* section of the `spaghetti` [documentation](https://pysal.org/spaghetti/tutorials.html).

# Planned Enhancements

As with any software project, there are always plans for further improvements and additional functionality. Three such major enhancements are described here. The first addition will likely be network partitioning through use of voronoi diagrams generated in network space. Network-constrained voronoi diagrams can be utilized as tools for analysis in and of themselves and can also be input for further analysis, such as the voronoi extension of the Network *K* function [@okabe2012]. Second, the current algorithm for allocating observations to a network within `spaghetti` allows for points to be snapped to a single location along the nearest network segment. While this is ideal for concrete observations, such as individual crime incidents, multiple network connections for abstract network events, such as census tract centroids, may be more appropriate [@gaboardi2020a]. Finally, the core functionality of `spaghetti` is nearly entirely written with pure Python data structures, which are excellent for code readability and initial development but generally suffer in terms of performance. There are currently several functions that can be utilized with an optional `geopandas` installation, however, further integration with the `pandas` stack has the potential to greatly improve performance.

# Concluding Remarks

Network-constrained spatial analysis is an important facet of scientific inquiry, especially within the social and geographic sciences [@Marshall2018]. Being able to perform this type of spatial analysis with a well-documented and tested open-source software package further facilitates fully reproducible and open science. With these motivations and core values, the `spaghetti` developers and wider PySAL team look forward to creating and supporting research into the future.

# Acknowledgements

Firstly, we would like to thank all the contributors to, and users of, this package. We would also like to acknowledge Jay Laura, who was the original lead developer of this package (`pysal.network`) prior to the introduction of the PySAL 2.0 ecosystem. The development of this package was partially supported by the [Atlanta Research Data Center](https://atlantardc.wordpress.com) and National Science Foundation Award [#1825768](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1825768).

# References
