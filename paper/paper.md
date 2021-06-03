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
date: 2 April 2021
bibliography: paper.bib
---


# Summary

The role spatial networks, such as streets, play on the human experience cannot be overstated. All of our daily activities fall along, or in close proximity to, roads, bike paths, and subway systems to name a few. Therefore, when performing spatial analysis in many cases considering network space, as opposed to Euclidean space, allows for a more precise representation of daily human action and movement patterns. For example, people generally cannot get to work by driving in a straight line directly from their home, but move along paths within networks. To this end, `spaghetti` (**spa**tial **g**rap**h**s: n**et**works, **t**opology, & **i**nference), a sub-module embedded in the wider [PySAL](https://pysal.org) ecosystem, was developed to address network-centric research questions with a strong focus on spatial analysis [@pysal2007;@pysal2015;@gaboardi2018].

Through `spaghetti`, first, network objects can be created and analysed from collections of line data by various means including reading in a shapefile or passing in a `geopandas.GeoDataFrame` at which time the line data are assigned network topology. Second, `spaghetti` provides computational tools to support statistical analysis of so-called network-based events along many different types of previously loaded networks. Network based-events or near-network observations are events that happen along spatial networks in our daily lives, i.e., locations of trees along footpaths, biking accidents along roads or locations of coffee shops along streets. As with `spaghetti.Network` objects, `spaghetti.PointPattern` objects can be [created from](https://pysal.org/spaghetti/generated/spaghetti.PointPattern.html#spaghetti.PointPattern) shapefiles, `geopandas.GeoDataFrame` objects or single `libpysal.cg.Point `objects. Near-network observations can then be snapped to nearest network segments enabling the calculation of observation distance matrices. Third, these observation distance matrices can be used both within `spaghetti` to perform clustering analysis or serve as input for other network-centric problems (e.g., optimal routing), or within the wider PySAL ecosystem to perform exploratory spatial analysis with [`esda`](https://pysal.org/esda/). Finally, `spaghetti`’s network elements (vertices and arcs) can also be extracted as `geopandas.GeoDataFrame` objects for visualization and integrated into further spatial statistical analysis within PySAL (e.g., `esda`).


# Related Work & Statement of Need

The most well-known network analysis package within the Python scientific stack is [NetworkX](https://networkx.github.io) [@hagberg2008], which can be used for modelling any type of complex network (e.g., social, spatial, etc.). [OSMnx](https://osmnx.readthedocs.io/en/stable/) [@boeing2017] is built on top of NetworkX and queries [OpenStreetMap](https://openstreetmap.org) for modelling street networks with resultant network objects returned within a `geopandas.GeoDataFrame` [@geopandas2021]. Another package, [`pandana`](https://github.com/UDST/pandana) [@foti2012generalized], is built on top of `pandas` [@mckinney-proc-scipy-2010;@reback2021pandas] with a focus on shortest path calculation and accessibility measures. Within the realm of Python, the functionality provided by  [`snkit`](https://github.com/tomalrussell/snkit) [@russell2019] is most comparable to `spaghetti`, though it's main purpose is the processing of raw line data into clean network objects. Outside of Python, [SANET](http://sanet.csis.u-tokyo.ac.jp) [@okabe2006a] is the most closely related project to `spaghetti`, however, it is not written in Python and provides a GUI plugin for GIS software such as QGIS. Moreover, SANET is not fully open source. While all the libraries above are important for network-based research, `spaghetti` was created and has evolved in line with the Python Spatial Analysis Library ecosystem for the specific purpose of utilizing the functionality of spatial weights in [`libpysal`](https://pysal.org/libpysal/) for generating network segment contiguity objects.


# Current Functionality

Considering the related projects in the [Related Work & Statement of Need section](# Related Work & Statement of Need) detailed above, `spaghetti` fills a niche for not only the processing of spatial network objects, but also post-processing analysis. In other words, this package can be used to study the network *itself* or provide the foundation for studying network-based phenomena, such as crimes along city streets, all within a fully open-source environment. Considering this, the primary purpose of `spaghetti` is creating network objects: collections of vertices and arcs, and their topological relationships. The creation of a network object is realized through the following general steps:

 1. read in line data or create features (regular lattices)
 1. generate the network representation 
 1. extract contiguity weights (if desired) as show in \autoref{fig:gridweights}
 1. identify connected components (if desired)
 1. extract graph representation of the network (if desired)

After the creation of a base network object it can be manipulated, analyzed, and utilized as the input for subsequent modelling scenarios. The following are several such examples:

 * allocating (snapping) observation point patterns to the network (see \autoref{fig:pointsnapmoran})
 * calculating all neighbor distance matrices
    * point type A to point type A (auto)
    * point type A to point type B (cross)
 * utilizing observation counts on network segments and network spatial weights within the [Moran’s *I*](https://pysal.org/spaghetti/generated/spaghetti.Network.html#spaghetti.Network.Moran) attribute to analyze global spatial autocorrelation [@moran:_cliff81;@esda:_2019] as seen in \autoref{fig:pointsnapmoran}
 * simulating point patterns that can be used within the  [*K* function](https://pysal.org/spaghetti/generated/spaghetti.Network.html#spaghetti.Network.GlobalAutoK) attribute for cluster analysis [@osullivan_unwin2010.ch5;@okabe2012]
 * splitting the network into (nearly) uniform segments
 * extracting features as `geopandas.GeoDataFrame` objects:
    * network arcs, vertices and point patterns 
    * largest/longest components
    * shortest paths
    * minimum/maximum spanning trees

The following two demonstrations show several functionalities mentioned above, including feature creation, network instantiation, network allocation, and feature extraction, along with supplementary plots in \autoref{fig:gridweights} and \autoref{fig:pointsnapmoran}.

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

![A 4x4 regular lattice with network arcs in gray and vertices in red. Connectivity is demonstrated with `libpysal` spatial weights, which are plotted over the network in black [@libpysal2020]. \label{fig:gridweights}](figs/spaghetti_network.png){ width=50% }

```python
import spaghetti, libpysal, matplotlib
# create a network from a line shapefile 
ntw = spaghetti.Network(in_data=libpysal.examples.get_path("streets.shp"))
# associate point observations with the network
pp_name = "schools"
pp_shp = libpysal.examples.get_path("%s.shp" % pp_name)
ntw.snapobservations(pp_shp, pp_name, attribute=True)
# calculation global spatial autocorrelation (Moran's I)
moran, yaxis = ntw.Moran(pp_name)
# extract network elements & observations
arcs_df = spaghetti.element_as_gdf(ntw, arcs=True)
schools = spaghetti.element_as_gdf(ntw, pp_name=pp_name)
schools_snapped = spaghetti.element_as_gdf(ntw, pp_name=pp_name, snapped=True)
# plot
base_kws = {"figsize":(7, 7), "lw":3, "color":"k", "zorder":0}
base = arcs_df.plot(**base_kws, alpha=.35)
schools.plot(ax=base, fc="b", ec="k", markersize=100, zorder=1, alpha=.5)
schools_snapped.plot(ax=base, fc="g", ec="k", markersize=50, zorder=2)
matplotlib.pyplot.title(f"Moran's $I$: {round(moran.I, 3)}", size="xx-large")
```

![Demonstrating the creation of a network and point pattern from shapefiles, followed by spatial autocorrelation analysis. A shapefile of school locations (blue) is read in and the points are snapped to the nearest network segments (green). A Moran's *I* statistic of -0.026 indicates near complete spatial randomness, though slightly dispersed. \label{fig:pointsnapmoran}](figs/spaghetti_pointpattern_moran.png){ width=50% }

The overview presented here provides a high-level summary of functionality. More detailed examples and applications can be found in the *Tutorials* section of the `spaghetti` [documentation](https://pysal.org/spaghetti/tutorials.html).


# Planned Enhancements

As with any software project, there are always plans for further improvements and additional functionality. Four such major enhancements are described here. The first addition will likely be network partitioning through use of voronoi diagrams generated in network space. Network-constrained voronoi diagrams can be utilized as tools for analysis in and of themselves and can also be input for further analysis, such as the voronoi extension of the Network *K* function [@okabe2012]. Second, the current algorithm for allocating observations to a network within `spaghetti` allows for points to be snapped to a single location along the nearest network segment. While this is ideal for concrete observations, such as individual crime incidents, multiple network connections for abstract network events, such as census tract centroids, may be more appropriate [@gaboardi2020a]. Third, the core functionality of `spaghetti` is nearly entirely written with pure Python data structures, which are excellent for code readability and initial development but generally suffer in terms of performance. There are currently several functions that can be utilized with an optional `geopandas` installation, however, further integration with the `pandas` stack has the potential to greatly improve performance. Finally, `spaghetti` developers will assess together with PySAL developers how to best support visualization and visual analysis targeted towards `spaghetti` network objects, implemented within visualization packages like [`splot`](https://splot.readthedocs.io/en/latest/?badge=latest) or [`mapclassify`](https://pysal.org/mapclassify/) and exposed as high level plotting functionality in `spaghetti`  [@splot:_Lumnitz2020].


# Concluding Remarks

Network-constrained spatial analysis is an important facet of scientific inquiry, especially within the social and geographic sciences [@Marshall2018]. Being able to perform this type of spatial analysis with a well-documented and tested open-source software package further facilitates fully reproducible and open science. With these motivations and core values, the `spaghetti` developers and wider PySAL team look forward to creating and supporting research into the future.


# Acknowledgements

Firstly, we would like to thank all the contributors to, and users of, this package. We would also like to acknowledge Jay Laura, who was the original lead developer of this package (`pysal.network`) prior to the introduction of the PySAL 2.0 ecosystem. The development of this package was partially supported by the [Atlanta Research Data Center](https://atlantardc.wordpress.com) and National Science Foundation Award [#1825768](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1825768).


# References



