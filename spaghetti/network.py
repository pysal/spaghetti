import contextlib
import copy
import os
import pickle
import warnings
from collections import OrderedDict, defaultdict
from itertools import islice

import esda
import numpy
from libpysal import cg, weights
from libpysal.common import requires

from . import util
from .analysis import GlobalAutoK

try:
    from libpysal import open as _open
except ImportError:
    import libpysal

    _open = libpysal.io.open


__all__ = ["Network", "PointPattern", "GlobalAutoK"]

SAME_SEGMENT = (-0.1, -0.1)


dep_msg = (
    "The next major release of pysal/spaghetti (2.0.0) will "
    "drop support for all ``libpysal.cg`` geometries. This change "
    "is a first step in refactoring ``spaghetti`` that is "
    "expected to result in dramatically reduced runtimes for "
    "network instantiation and operations. Users currently "
    "requiring network and point pattern input as ``libpysal.cg`` "
    "geometries should prepare for this simply by converting "
    "to ``shapely`` geometries."
)
warnings.warn(dep_msg, FutureWarning, stacklevel=1)


class Network:
    """Spatially-constrained network representation and analytical
    functionality. Naming conventions are as follows, (1) arcs and
    vertices for the full network object, and (2) edges and nodes for
    the simplified graph-theoretic object. The term 'link' is used to
    refer to a network arc or a graph edge.

    Parameters
    ----------
    in_data : {str, iterable, libpysal.cg.Chain, geopandas.GeoDataFrame}
        The input geographic data. Either (1) a path to a shapefile
        (str); (2) an iterable containing ``libpysal.cg.Chain``
        objects; (3) a single ``libpysal.cg.Chain``; or
        (4) a ``geopandas.GeoDataFrame``.
    vertex_sig : int
        Round the x and y coordinates of all vertices to ``vertex_sig``
        significant digits (combined significant digits on the left and
        right of the decimal place). Default is 11. Set to ``None`` for
        no rounding.
    unique_arcs : bool
        If ``True`` (default), keep only unique arcs (i.e., prune
        out any duplicated arcs). If ``False`` keep all segments.
    extractgraph : bool
        If ``True``, extract a graph-theoretic object with no degree 2
        nodes. Default is ``True``.
    w_components : bool
        Set to ``False`` to not record connected components from a
        ``libpysal.weights.W`` object. Default is ``True``.
    weightings : {dict, bool}
        If dict, lists of weightings for each arc. If bool,
        ``True`` flags ``self.arc_lengths`` as the weightings,
        ``False`` sets no weightings. Default is ``False``.
    weights_kws : dict
        Keyword arguments for ``libpysal.weights.W``. Default is ``dict()``.
    vertex_atol : {int, None}
        Precision for vertex absolute tolerance. Round vertex coordinates to
        ``vertex_atol`` decimal places. Default is ``None``. **ONLY** change
        the default when there are known issues with digitization.

    Attributes
    ----------
    adjacencylist : list
        List of lists storing vertex adjacency.
    vertex_coords : dict
        Keys are vertex IDs and values are :math:`(x,y)` coordinates of the vertices.
    vertex_list : list
        List of vertex IDs.
    vertices : dict
        Keys are tuples of vertex coords and values are the vertex ID.
    arcs : list
        List of arcs, where each arc is a sorted tuple
        of vertex IDs.
    arc_lengths : dict
        Keys are tuples of sorted vertex IDs representing an arc and
        values are the length.
    pointpatterns : dict
        Keys are a string name of the pattern and values are
        ``PointPattern`` class instances.
    distance_matrix : numpy.ndarray
        All network vertices (non-observations) distance matrix. Distances
        between vertices in disparate components are recorded as ``inf``
        by default.
    network_trees : dict
        Keys are the vertex IDs (``int``). Values are dictionaries
        with the keys being the IDs of the destination vertex
        and values being lists of vertices along the shortest path.
        If the destination vertex is a) the origin or b)
        unreachable (disparate component) it is listed as itself being the
        neighbor.
    edges : list
        Tuples of graph edge IDs.
    edge_lengths : dict
        Keys are the graph edge IDs (``tuple``). Values are the
        graph edge length (``float``).
    non_articulation_points : list
        All vertices with degree 2 that are not in an isolated
        island ring (loop) component.
    w_network : libpysal.weights.W
        Weights object created from the network arcs.
    network_n_components : int
        Count of connected components in the network.
    network_fully_connected : bool
        ``True`` if the network representation is a single connected
        component, otherwise ``False``.
    network_component_labels : numpy.ndarray
        Component labels for network arcs.
    network_component2arc : dict
        Lookup in the form {int: list} for arcs comprising network
        connected components keyed by component labels with arcs in
        a list as values.
    network_component_lengths : dict
        Length of each network component (keyed by component label).
    network_longest_component : int
        The ID of the longest component in the network. This is not
        necessarily equal to ``network_largest_component``.
    network_component_vertices : dict
        Lookup in the form {int: list} for vertices comprising network
        connected components keyed by component labels with vertices in
        a list as values.
    network_component_vertex_count : dict
        The number of vertices in each network component
        (keyed by component label).
    network_largest_component : int
        The ID of the largest component in the network. Within ``spaghetti``
        the largest component is the one with the most vertices. This is not
        necessarily equal to ``network_longest_component``.
    network_component_is_ring : dict
        Lookup in the form {int: bool} keyed by component labels with values
        as ``True`` if the component is a closed ring, otherwise ``False``.
    w_graph : libpysal.weights.W
        Weights object created from the graph edges.
    graph_n_components : int
        Count of connected components in the network.
    graph_fully_connected : bool
        ``True`` if the graph representation is a single connected
        component, otherwise ``False``.
    graph_component_labels : numpy.ndarray
        Component labels for graph edges.
    graph_component2edge : dict
        Lookup in the form {int: list} for edges comprising graph connected
        components keyed by component labels with edges in a list
        as values.
    graph_component_lengths : dict
        Length of each graph component (keyed by component label).
    graph_longest_component : int
        The ID of the longest component in the graph. This is not
        necessarily equal to ``graph_largest_component``.
    graph_component_vertices : dict
        Lookup in the form {int: list} for vertices comprising graph
        connected components keyed by component labels with vertices in
        a list as values.
    graph_component_vertex_count : dict
        The number of vertices in each graph component
        (keyed by component label).
    graph_largest_component : int
        The ID of the largest component in the graph. Within ``spaghetti``
        the largest component is the one with the most vertices. This is not
        necessarily equal to ``graph_longest_component``.
    graph_component_is_ring : dict
        Lookup in the form {int: bool} keyed by component labels with values as
        ``True`` if the component is a closed ring, otherwise ``False``.

    Notes
    -----

    **Important**: The core procedure for generating network representations is
    performed within the ``_extractnetwork()`` method. Here it is important to note
    that a ``spaghetti.Network`` instance is built up from the individual,
    constituent euclidean units of each line segment object. Therefore, the resulting
    network structure will generally have (1) more vertices and links than may expected,
    and, (2) many degree-2 vertices, which differs from a truly graph-theoretic object.
    This is demonstrated in the
    `Caveats Tutorial <https://pysal.org/spaghetti/notebooks/caveats.html>`_.

    See :cite:`Cliff1981`, :cite:`Tansel1983a`,
    :cite:`AhujaRavindraK`, :cite:`Labbe1995`,
    :cite:`Kuby2009`, :cite:`Barthelemy2011`,
    :cite:`daskin2013`, :cite:`Okabe2012`,
    :cite:`Ducruet2014`, :cite:`Weber2016`, for more in-depth discussion on
    spatial networks, graph theory, and location along networks.
    For related network-centric software see
    `Snkit <https://github.com/tomalrussell/snkit>`_ :cite:`tom_russell_2019_3379659`,
    `SANET <http://sanet.csis.u-tokyo.ac.jp>`_ :cite:`Okabe2006a`,
    `NetworkX <https://networkx.github.io>`_ :cite:`Hagberg2008`,
    `Pandana <http://udst.github.io/pandana/>`_ :cite:`Foti2012`,
    and `OSMnx <https://osmnx.readthedocs.io/en/stable/>`_ :cite:`Boeing2017`.

    Examples
    --------

    Create an instance of a network.

    >>> import spaghetti
    >>> from libpysal import examples
    >>> streets_file = examples.get_path("streets.shp")
    >>> ntw = spaghetti.Network(in_data=streets_file)

    Fetch the number connected components in the network.

    >>> ntw.network_n_components
    1

    Unique component labels in the network.

    >>> import numpy
    >>> list(numpy.unique(ntw.network_component_labels))
    [0]

    Show whether each component of the network is an isolated ring (or not).

    >>> ntw.network_component_is_ring
    {0: False}

    Show how many network arcs are associated with the component.

    >>> arcs = len(ntw.network_component2arc[ntw.network_component_labels[0]])
    >>> arcs
    303

    Do the same as above, but for the graph-theoretic representation
    of the network object.

    >>> ntw.graph_n_components
    1
    >>> list(numpy.unique(ntw.graph_component_labels))
    [0]
    >>> ntw.graph_component_is_ring
    {0: False}
    >>> edges = len(ntw.graph_component2edge[ntw.graph_component_labels[0]])
    >>> edges
    179

    The number of arcs in the network is always greater than or equal
    to the number of edges in the graph-theoretic representation.

    >>> arcs >= edges
    True

    Snap point observations to the network with attribute information.

    >>> crimes_file = examples.get_path("crimes.shp")
    >>> ntw.snapobservations(crimes_file, "crimes", attribute=True)

    And without attribute information.

    >>> schools_file = examples.get_path("schools.shp")
    >>> ntw.snapobservations(schools_file, "schools", attribute=False)

    Show the point patterns associated with the network.

    >>> ntw.pointpatterns.keys()
    dict_keys(['crimes', 'schools'])

    """

    def __init__(
        self,
        in_data=None,
        vertex_sig=11,
        unique_arcs=True,
        extractgraph=True,
        w_components=True,
        weightings=False,
        weights_kws=dict(),
        vertex_atol=None,
    ):
        # do this when creating a clean network instance from a
        # shapefile or a geopandas.GeoDataFrame, otherwise a shell
        # network instance is created (see `split_arcs()` method)
        if in_data is not None:
            # set parameters as attributes
            self.in_data = in_data
            self.vertex_sig = vertex_sig
            self.vertex_atol = vertex_atol
            self.unique_arcs = unique_arcs

            self.adjacencylist = defaultdict(list)
            self.vertices = {}

            # initialize network arcs and arc_lengths
            self.arcs = []
            self.arc_lengths = {}

            # initialize pointpatterns
            self.pointpatterns = {}

            # spatial representation of the network
            self._extractnetwork()
            self.arcs = sorted(self.arcs)
            self.vertex_coords = {v: k for k, v in self.vertices.items()}

            # extract connected components
            if w_components:
                as_graph = False
                network_weightings = False

                if weightings:
                    # set network arc weights to length if weights are
                    # desired, but no other input in given
                    weightings = self.arc_lengths
                    network_weightings = True

                # extract contiguity weights from libpysal
                self.w_network = self.contiguityweights(
                    graph=as_graph,
                    weightings=weightings,
                    weights_kws=weights_kws,
                )
                # identify connected components from the `w_network`
                self.identify_components(self.w_network, graph=as_graph)

            # extract the graph -- repeat similar as above
            # for extracting the network
            if extractgraph:
                self.extractgraph()

                if w_components:
                    as_graph = True

                    if network_weightings:
                        weightings = self.edge_lengths

                    self.w_graph = self.contiguityweights(
                        graph=as_graph,
                        weightings=weightings,
                        weights_kws=weights_kws,
                    )
                    self.identify_components(self.w_graph, graph=as_graph)

            # sorted list of vertex IDs
            self.vertex_list = sorted(self.vertices.values())

    def _round_sig(self, v):
        """Used internally to round the vertex to a set number of
        significant digits. If ``sig`` is set to 4, then the following
        are some possible results for a coordinate are as follows.
        (1) 0.0xxxx, (2) 0.xxxx, (3) x.xxx, (4) xx.xx,
        (5) xxx.x, (6) xxxx.0, (7) xxxx0.0

        Parameters
        ----------
        v : tuple
            Coordinate (x,y) of the vertex.

        """

        # set the number of significant digits
        sig = self.vertex_sig

        # simply return vertex (x,y) coordinates
        if sig is None:
            return v

        # for each coordinate in a coordinate pair
        # if the coordinate location is (0.0) simply return zero
        # else -- (1) take the absolute value of `val`; (2) take the
        # base 10 log for [1]; (3) take the floor of [2]; (4) convert
        # [3] into a negative integer; (5) add `sig - 1` to [4];
        # (6) round `val` by [5]
        out_v = [
            val
            if val == 0
            else round(val, -int(numpy.floor(numpy.log10(numpy.fabs(val)))) + (sig - 1))
            for val in v
        ]

        if self.vertex_atol:
            out_v = [round(v, self.vertex_atol) for v in out_v]

        return tuple(out_v)

    def identify_components(self, w, graph=False):
        """Identify connected component information from a
        ``libpysal.weights.W`` object

        Parameters
        ----------
        w : libpysal.weights.W
            Weights object created from the network segments (either
            raw or graph-theoretic).
        graph : bool
            Flag for a raw network (``False``) or graph-theoretic network
            (``True``). Default is ``False``.

        """

        # flag network (arcs) or graph (edges)
        if graph:
            links = self.edges
            obj_type = "graph_"
        else:
            links = self.arcs
            obj_type = "network_"

        # connected component count and labels
        n_components = w.n_components
        component_labels = w.component_labels

        # is the network a single, fully-connected component?
        fully_connected = bool(n_components == 1)

        # link to component lookup
        link2component = dict(zip(links, component_labels))

        # component ID lookups: links, lengths, vertices, vertex counts
        component2link = {}
        component_lengths = {}
        component_vertices = {}
        component_vertex_count = {}
        cp_labs_ = set(w.component_labels)
        l2c_ = link2component.items()
        for cpl in cp_labs_:
            component2link[cpl] = sorted([k for k, v in l2c_ if v == cpl])
            c2l_ = component2link[cpl]
            arclens_ = self.arc_lengths.items()
            component_lengths[cpl] = sum([v for k, v in arclens_ if k in c2l_])
            component_vertices[cpl] = {v for link in c2l_ for v in link}
            component_vertex_count[cpl] = len(component_vertices[cpl])

        # longest and largest components
        longest_component = max(component_lengths, key=component_lengths.get)
        largest_component = max(component_vertex_count, key=component_vertex_count.get)

        # component to ring lookup
        component_is_ring = {}
        adj_ = self.adjacencylist.items()
        for comp, verts in component_vertices.items():
            component_is_ring[comp] = False
            _2neighs = [len(neighs) == 2 for v, neighs in adj_ if v in verts]
            if all(_2neighs):
                component_is_ring[comp] = True

        # attribute label name depends on object type
        c2l_attr_name = "component2edge" if graph else "component2arc"

        # set all new variables into list
        extracted_attrs = [
            ["fully_connected", fully_connected],
            ["n_components", n_components],
            ["component_labels", component_labels],
            [c2l_attr_name, component2link],
            ["component_lengths", component_lengths],
            ["component_vertices", component_vertices],
            ["component_vertex_count", component_vertex_count],
            ["longest_component", longest_component],
            ["largest_component", largest_component],
            ["component_is_ring", component_is_ring],
        ]

        # iterate over list and set attribute with
        # either "network" or "graph" extension
        for attr_str, attr in extracted_attrs:
            setattr(self, obj_type + attr_str, attr)

    def _extractnetwork(self):
        """Used internally to extract a network."""

        # initialize vertex count
        vertex_count = 0

        # determine input network data type
        in_dtype = str(type(self.in_data)).split("'")[1]
        is_libpysal_chains = False
        supported_iterables = ["list", "tuple", "numpy.ndarray"]
        # type error message
        msg = "'{}' not supported for network instantiation."

        # set appropriate geometries
        if in_dtype == "str":
            shps = _open(self.in_data)
        elif in_dtype in supported_iterables:
            shps = self.in_data
            shp_type = str(type(shps[0])).split("'")[1]
            if shp_type == "libpysal.cg.shapes.Chain":
                is_libpysal_chains = True
            else:
                raise TypeError(msg.format(shp_type))
        elif in_dtype == "libpysal.cg.shapes.Chain":
            shps = [self.in_data]
            is_libpysal_chains = True
        elif in_dtype == "geopandas.geodataframe.GeoDataFrame":
            shps = self.in_data.geometry
        else:
            raise TypeError(msg.format(in_dtype))

        # iterate over each record of the network lines
        for shp in shps:
            # if the segments are native pysal geometries
            if is_libpysal_chains:
                vertices = shp.vertices
            else:
                # fetch all vertices between euclidean segments
                # in the line record -- these vertices are
                # coordinates in an (x, y) tuple.
                vertices = weights._contW_lists._get_verts(shp)

            # iterate over each vertex (v)
            for i, v in enumerate(vertices[:-1]):
                # -- For vertex 1
                # adjust precision -- this was originally
                # implemented to handle high-precision
                # network network vertices
                v = self._round_sig(v)

                # when the vertex already exists in lookup
                # set it as the current `vid`
                try:
                    vid = self.vertices[v]
                # when the vertex is not present in the lookup
                # add it and adjust vertex count
                except KeyError:
                    self.vertices[v] = vid = vertex_count
                    vertex_count += 1

                # -- For vertex 2
                # repeat the steps above for vertex 1
                v2 = self._round_sig(vertices[i + 1])
                try:
                    nvid = self.vertices[v2]
                except KeyError:
                    self.vertices[v2] = nvid = vertex_count
                    vertex_count += 1

                # records vertex 1 and vertex 2 adjacency
                self.adjacencylist[vid].append(nvid)
                self.adjacencylist[nvid].append(vid)

                # Sort the edges so that mono-directional
                # keys can be stored.
                arc_vertices = sorted([vid, nvid])
                arc = tuple(arc_vertices)

                # record the euclidean arc within the network
                self.arcs.append(arc)

                # record length
                length = util.compute_length(v, vertices[i + 1])
                self.arc_lengths[arc] = length

        if self.unique_arcs:
            # Remove duplicate edges and duplicate adjacent nodes.
            self.arcs = list(set(self.arcs))
            for k, v in self.adjacencylist.items():
                self.adjacencylist[k] = list(set(v))

    def extractgraph(self):
        """Using the existing network representation, create a
        graph-theoretic representation by removing all vertices with a
        neighbor incidence of two (non-articulation points). That is, we
        assume these vertices are bridges between vertices with higher
        or lower incidence.
        """

        # initialize edges and edge_lengths
        self.edges = []
        self.edge_lengths = {}

        # find all vertices with degree 2 that are not in an isolated
        # island ring (loop) component. These are non-articulation
        # points on the graph representation
        non_articulation_points = self._yield_napts()
        # retain non_articulation_points as an attribute
        self.non_articulation_points = list(non_articulation_points)

        # start with a copy of the spatial representation and
        # iteratively remove edges deemed to be segments
        self.edges = copy.deepcopy(self.arcs)
        self.edge_lengths = copy.deepcopy(self.arc_lengths)

        # mapping all the 'network arcs' contained within a single
        # 'graph represented' edge
        self.arcs_to_edges = {}

        # build up bridges "rooted" on the initial
        # non-articulation points
        bridge_roots = []

        # iterate over all vertices that are not contained within
        # isolated loops that have a degree of 2
        for s in non_articulation_points:
            # initialize bridge with an articulation point
            bridge = [s]

            # fetch all vertices adjacent to point `s`
            # that are also degree 2
            neighbors = self._yieldneighbor(s, non_articulation_points, bridge)
            while neighbors:
                # extract the current node in `neighbors`
                cnode = neighbors.pop()
                # remove it from `non_articulation_points`
                non_articulation_points.remove(cnode)
                # add it to bridge
                bridge.append(cnode)
                # fetch neighbors for the current node
                newneighbors = self._yieldneighbor(
                    cnode, non_articulation_points, bridge
                )
                # add the new neighbors back into `neighbors`
                neighbors += newneighbors

            # once all potential neighbors are exhausted add the
            # current bridge of non-articulation points to the
            # list of rooted bridges
            bridge_roots.append(bridge)

        # iterate over the list of newly created rooted bridges
        for bridge in bridge_roots:
            # if the vertex is only one non-articulation
            # point in the bridge
            if len(bridge) == 1:
                # that the singular element of the bridge
                n = self.adjacencylist[bridge[0]]
                # and create a new graph edge from it
                new_edge = tuple(sorted([n[0], n[1]]))

                # identify the arcs to be removed
                e1 = tuple(sorted([bridge[0], n[0]]))
                e2 = tuple(sorted([bridge[0], n[1]]))

                # remove the network arcs (spatial) from the
                # graph-theoretic representation
                self.edges.remove(e1)
                self.edges.remove(e2)

                # remove the former network arc lengths from the
                # graph edge lengths lookup
                length_e1 = self.edge_lengths[e1]
                length_e2 = self.edge_lengths[e2]
                self.edge_lengths.pop(e1, None)
                self.edge_lengths.pop(e2, None)

                # and add the new edge length in their place
                self.edge_lengths[new_edge] = length_e1 + length_e2

                # update the pointers
                self.arcs_to_edges[e1] = new_edge
                self.arcs_to_edges[e2] = new_edge

            # if there are more than one vertices in the bridge
            else:
                cumulative_length = 0
                start_end = {}

                # initialize a redundant set of bridge edges
                redundant = set()

                # iterate over the current bridge
                for b in bridge:
                    # iterate over each node in the bridge
                    for n in self.adjacencylist[b]:
                        # start the bridge with this node
                        if n not in bridge:
                            start_end[b] = n
                        # or create a redundant edge with the current
                        # node and `b`
                        else:
                            redundant.add(tuple(sorted([b, n])))

                # initialize a new graph edge
                new_edge = tuple(sorted(start_end.values()))

                # add start_end redundant edge
                for k, v in start_end.items():
                    redundant.add(tuple(sorted([k, v])))

                # remove all redundant network arcs while
                # adjusting the graph edge lengths lookup
                # and the edges_to_arcs lookup
                for r in redundant:
                    self.edges.remove(r)
                    cumulative_length += self.edge_lengths[r]
                    self.edge_lengths.pop(r, None)
                    self.arcs_to_edges[r] = new_edge

                # finally, add the new cumulative edge length
                self.edge_lengths[new_edge] = cumulative_length

            # add the updated graph edge
            self.edges.append(new_edge)

        # converted the graph edges into a sorted set to prune out
        # duplicate graph edges created during simplification
        self.edges = sorted(set(self.edges))

    def _yield_napts(self):
        """Find all nodes with degree 2 that are not in an isolated
        island ring (loop) component. These are non-articulation
        points on the graph representation.

        Returns
        -------
        napts : list
            Non-articulation points on a graph representation.

        """

        # non-articulation points
        napts = set()

        # network vertices remaining to evaluate
        unvisted = set(self.vertices.values())

        while unvisted:
            # iterate over each component
            for component_id, ring in self.network_component_is_ring.items():
                # evaluate for non-articulation points
                napts, unvisted = self._evaluate_napts(
                    napts, unvisted, component_id, ring
                )

        # convert set of non-articulation points into list
        napts = list(napts)

        return napts

    def _evaluate_napts(self, napts, unvisited, component_id, ring):
        """Evaluate one connected component in a network for
        non-articulation points (``napts``) and return an updated set of
        ``napts`` and unvisted vertices.

        Parameters
        ----------
        napts : set
            Non-articulation points (``napts``) in the network. The
            ``napts`` here do not include those within an isolated
            loop island.
        unvisited : set
            Vertices left to evaluate in the network.
        component_id : int
            ID for the network connected component for the
            current iteration of the algorithm.
        ring : bool
            Network component is isolated island loop ``True`` or
            not ``False``.

        Returns
        -------
        napts : set
            Updated ``napts`` object.
        unvisited : set
            Updated ``napts`` object.

        """

        # iterate over each `edge` of the `component`
        for component in self.network_component2arc[component_id]:
            # each `component` has two vertices
            for vertex in component:
                # if `component` is not an isolated island
                # and `vertex` has exactly 2 neighbors,
                # add `vertex` to `napts`
                if not ring and len(self.adjacencylist[vertex]) == 2:
                    napts.add(vertex)

                # remove `vertex` from `unvisited` if
                # it is still in the set else move along to
                # the next iteration
                with contextlib.suppress(KeyError):
                    unvisited.remove(vertex)

        return napts, unvisited

    def _yieldneighbor(self, vtx, arc_vertices, bridge):
        """Used internally, this method traverses a bridge arc
        to find the source and destination nodes.

        Parameters
        ----------
        vtx : int
            The vertex ID.
        arc_vertices : list
            All non-articulation points (``napts``) in the network.
            These are referred to as degree-2 vertices.
        bridge : list
            Inital bridge list containing only ``vtx``.

        Returns
        -------
        nodes : list
            Vertices to keep (articulation points). These elements are
            referred to as nodes.

        """

        # instantiate empty lis to fill with network articulation
        # points (nodes with a degree of 1 [endpoints] or greater
        # than 2 [intersections])
        nodes = []

        # get all nodes adjacent to `vtx` that are not in the
        # set of 'bridge' vertices
        for i in self.adjacencylist[vtx]:
            if i in arc_vertices and i not in bridge:
                nodes.append(i)

        return nodes

    def contiguityweights(
        self, graph=True, weightings=None, from_split=False, weights_kws=dict()
    ):
        """Create a contiguity-based ``libpysal.weights.W`` object.

        Parameters
        ----------
        graph : bool
            Controls whether the ``libpysal.weights.W`` is generated
            using the spatial representation (``False``) or the graph
            representation (``True``). Default is ``True``.
        weightings : {dict, None}
            Dictionary of lists of weightings for each arc/edge. Default is ``None``.
        from_split : bool
            Flag for whether the method is being called from within
            ``split_arcs()`` (``True``) or not (``False``). Default is ``False``.
        weights_kws : dict
            Keyword arguments for ``libpysal.weights.W``. Default is ``dict()``.

        Returns
        -------
         W : libpysal.weights.W
            A ``W`` representing the binary adjacency of the network.

        See also
        --------

        libpysal.weights.W

        Examples
        --------

        Instantiate a network.

        >>> import spaghetti
        >>> from libpysal import examples
        >>> import numpy
        >>> ntw = spaghetti.Network(examples.get_path("streets.shp"))

        Snap point observations to the network with attribute information.

        >>> ntw.snapobservations(
        ...     examples.get_path("crimes.shp"), "crimes", attribute=True
        ... )

        Find counts per network arc.

        >>> counts = ntw.count_per_link(
        ...     ntw.pointpatterns["crimes"].obs_to_arc, graph=False
        ... )
        >>> counts[(50, 165)]
        4

        Create a contiguity-based ``W`` object.

        >>> w = ntw.contiguityweights(graph=False)
        >>> w.n, w.n_components
        (303, 1)

        Notes
        -----

        See :cite:`pysal2007` for more details.

        """

        # instantiate OrderedDict to record network link
        # adjacency which will be keyed by the link ID (a tuple)
        # with values being lists of tuples (contiguous links)
        neighbors = OrderedDict()

        # flag network (arcs) or graph (edges)
        links = self.edges if graph else self.arcs

        # if weightings are desired instantiate a dictionary
        # other ignore weightings
        _weights = {} if weightings else None

        # iterate over all links until all possibilities
        # for network link adjacency are exhausted
        working = True
        while working:
            # for each network link (1)
            for key in links:
                # instantiate a slot in the OrderedDict
                neighbors[key] = []

                if weightings:
                    _weights[key] = []

                # for each network link (2)
                for neigh in links:
                    # skip if comparing link to itself
                    if key == neigh:
                        continue

                    # if link(1) and link(2) share any vertex
                    # update neighbors adjacency
                    if (
                        key[0] == neigh[0]
                        or key[0] == neigh[1]
                        or key[1] == neigh[0]
                        or key[1] == neigh[1]
                    ):
                        neighbors[key].append(neigh)

                        # and add weights if desired
                        if weightings:
                            _weights[key].append(weightings[neigh])

                    # break condition
                    # -- everything is sorted, so we know when we have
                    # stepped beyond a possible neighbor
                    if key[1] > neigh[1]:
                        working = False

            if len(links) == 1 or from_split:
                working = False

        # call libpysal for `W` instance
        weights_kws["weights"] = _weights
        w = weights.W(neighbors, **weights_kws)

        return w

    def distancebandweights(
        self, threshold, n_processes=1, gen_tree=False, weights_kws=dict()
    ):
        """Create distance-based weights.

        Parameters
        ----------
        threshold : float
            Distance threshold value.
        n_processes : {int, str}
            Specify the number of cores to utilize. Default is 1 core.
            Use ``"all"`` to request all available cores.
            Specify the exact number of cores with an integer.
        gen_tree : bool
            Rebuild shortest path with ``True``, or skip with ``False``.
            Default is ``False``.
        weights_kws : dict
            Keyword arguments for ``libpysal.weights.W``. Default is ``dict()``.

        Returns
        -------
        w : libpysal.weights.W
            A ``W`` object representing the binary adjacency of
            the network.

        Notes
        -----

        See :cite:`AnselinRey2014` and :cite:`rey_open_2015` for more details
        regarding spatial weights.

        See also
        --------

        libpysal.weights.W

        Examples
        --------

        Instantiate an instance of a network.

        >>> import spaghetti
        >>> from libpysal import examples
        >>> import warnings
        >>> streets_file = examples.get_path("streets.shp")
        >>> ntw = spaghetti.Network(in_data=streets_file)

        Create a contiguity-based ``W`` object based on network distance, ``500``
        US feet in this case.

        >>> w = ntw.distancebandweights(
        ...     threshold=500, weights_kws=dict(silence_warnings=True)
        ... )

        Show the number of units in the ``W`` object.

        >>> w.n
        230

        There are 7 components in the ``W`` object.

        >>> w.n_components
        7

        There are ``8`` units with ``3`` neighbors in the ``W`` object.

        >>> w.histogram[-1]
        (8, 3)

        """

        # if the a vertex-to-vertex network distance matrix is
        # not present in the `network.Network` object; calculate
        # one at this point
        if not hasattr(self, "distance_matrix"):
            self.full_distance_matrix(n_processes, gen_tree=gen_tree)

        # identify all network vertices which are within the
        # `threshold` parameter
        neighbor_query = numpy.where(self.distance_matrix < threshold)

        # create an instance for recording neighbors which
        # inserts a new key if not present in object
        neighbors = defaultdict(list)

        # iterate over neighbors within the `threshold`
        # and record all network vertices as neighbors
        # if the vertex is not being compared to itself
        for i, n in enumerate(neighbor_query[0]):
            neigh = neighbor_query[1][i]
            if n != neigh:
                neighbors[n].append(neigh)

        # call libpysal for `W` instance
        w = weights.W(neighbors, **weights_kws)

        return w

    def snapobservations(self, in_data, name, idvariable=None, attribute=False):
        """Snap a point pattern shapefile to a network object. The
        point pattern is stored in the ``network.pointpattern``
        attribute of the network object.

        Parameters
        ----------
        in_data : {geopandas.GeoDataFrame, str}
            The input geographic data. Either (1) a path to a
            shapefile (str); or (2) a ``geopandas.GeoDataFrame``.
        name : str
            Name to be assigned to the point dataset.
        idvariable : str
            Column name to be used as the ID variable.
        attribute : bool
            Defines whether attributes should be extracted. ``True`` for
            attribute extraction. ``False`` for no attribute extraction.
            Default is ``False``.

        Notes
        -----

        See :cite:`doi:10.1111/gean.12211` for a detailed discussion on
        the modeling consequences of snapping points to spatial networks.

        Examples
        --------

        Instantiate a network.

        >>> import spaghetti
        >>> from libpysal import examples
        >>> streets_file = examples.get_path("streets.shp")
        >>> ntw = spaghetti.Network(in_data=streets_file)

        Snap observations to the network.

        >>> pt_str = "crimes"
        >>> in_data = examples.get_path(pt_str+".shp")
        >>> ntw.snapobservations(in_data, pt_str, attribute=True)

        Isolate the number of points in the dataset.

        >>> ntw.pointpatterns[pt_str].npoints
        287

        """

        # create attribute of `pointpattern` but instantiating a
        # `network.PointPattern` class
        self.pointpatterns[name] = PointPattern(
            in_data=in_data, idvariable=idvariable, attribute=attribute
        )

        # allocate the point observations to the nework
        self._snap_to_link(self.pointpatterns[name])

    def compute_distance_to_vertices(self, x, y, arc):
        """Given an observation on a network arc, return the distance
        to the two vertices that bound that end.

        Parameters
        ----------
        x : float
            The x-coordinate of the snapped point.
        y : float
            The y-coordinate of the snapped point.
        arc : tuple
            The (vtx0, vtx1) representation of the network arc.

        Returns
        -------
        d1 : float
            The distance to vtx0. Always the vertex with the lesser ID.
        d2 : float
            The distance to vtx1. Always the vertex with the greater ID.

        """

        # distance to vertex 1
        d1 = util.compute_length((x, y), self.vertex_coords[arc[0]])

        # distance to vertex 2
        d2 = util.compute_length((x, y), self.vertex_coords[arc[1]])

        return d1, d2

    def compute_snap_dist(self, pattern, idx):
        """Given an observation snapped to a network arc, calculate the
        distance from the original location to the snapped location.

        Parameters
        ----------
        pattern : spaghetti.PointPattern
            The point pattern object.
        idx : int
            The point ID.

        Returns
        -------
        dist : float
            The euclidean distance from original location to the snapped
            location.

        """

        # set of original (x,y) point coordinates
        loc = pattern.points[idx]["coordinates"]

        # set of snapped (x,y) point coordinate
        snp = pattern.snapped_coordinates[idx]

        # distance from the original location to
        # the snapped location along the network
        dist = util.compute_length(loc, snp)

        return dist

    def _snap_to_link(self, pointpattern):
        """Used internally to snap point observations to network arcs.

        Parameters
        ----------
        pointpattern : spaghetti.PointPattern
            The point pattern object.

        Returns
        -------
        obs_to_arc : dict
            Dictionary with arcs as keys and lists of points as values.
        arc_to_obs : dict
            Dictionary with point IDs as keys and arc tuples as values.
        dist_to_vertex : dict
            Dictionary with point IDs as keys and values as dictionaries
            with keys for vertex IDs and values as distances from point
            to vertex.
        dist_snapped : dict
            Dictionary with point IDs as keys and distance from point
            to the network arc that it is snapped.

        """

        # instantiate observations snapped coordinates lookup
        pointpattern.snapped_coordinates = {}

        # record throw-away arcs (pysal.cg.Chain) enumerator
        arcs_ = []

        # snapped(point)-to-arc lookup
        s2a = {}

        # iterate over network arc IDs
        for arc in self.arcs:
            # record the start and end of the arc
            head = self.vertex_coords[arc[0]]
            tail = self.vertex_coords[arc[1]]

            # create a pysal.cg.Chain object of the arc
            # and add it to the arcs enumerator
            arcs_.append(util._chain_constr(None, [head, tail]))

            # add the arc into the snapped(point)-to-arc lookup
            s2a[(head, tail)] = arc

        # instantiate crosswalks
        points = {}  # point ID to coordinates lookup
        obs_to_arc = {}  # observations to arcs lookup
        dist_to_vertex = {}  # distance to vertices lookup
        dist_snapped = {}  # snapped distance lookup

        # fetch and records point coordinates keyed by ID
        for point_idx, point in pointpattern.points.items():
            points[point_idx] = point["coordinates"]

        # snap point observations to the network
        snapped = util.snap_points_to_links(points, arcs_)

        # record obs_to_arc, dist_to_vertex, and dist_snapped
        # -- iterate over the snapped observation points
        for point_idx, snap_info in snapped.items():
            # fetch the x and y coordinate
            x, y = snap_info[1].tolist()

            # look up the arc from snapped(point)-to-arc
            arc = s2a[tuple(snap_info[0])]

            # add the arc key to observations to arcs lookup
            if arc not in obs_to_arc:
                obs_to_arc[arc] = {}

            # add the (x,y) coordinates of the original observation
            # point location to the observations to arcs lookup
            obs_to_arc[arc][point_idx] = (x, y)

            # add the (x,y) coordinates of the snapped observation
            # point location to the snapped coordinates lookup
            pointpattern.snapped_coordinates[point_idx] = (x, y)

            # calculate the distance to the left and right vertex
            # along the network link from the snapped point location
            d1, d2 = self.compute_distance_to_vertices(x, y, arc)

            # record the distances in the distance to vertices lookup
            dist_to_vertex[point_idx] = {arc[0]: d1, arc[1]: d2}

            # record the snapped distance
            dist_snapped[point_idx] = self.compute_snap_dist(pointpattern, point_idx)

        # instantiate observations to network vertex lookup
        obs_to_vertex = defaultdict(list)

        # iterate over the observations to arcs lookup
        for k, v in obs_to_arc.items():
            # record the left and right vertex ids
            keys = v.keys()
            obs_to_vertex[k[0]] = keys
            obs_to_vertex[k[1]] = keys

        # iterate over components and assign observations
        component_to_obs = {}
        for comp, _arcids in self.network_component2arc.items():
            component_to_obs[comp] = []
            for lk, odict in obs_to_arc.items():
                if lk in _arcids:
                    component_to_obs[comp].extend(list(odict.keys()))

        # set crosswalks as attributes of the `pointpattern` class
        pointpattern.obs_to_arc = obs_to_arc
        pointpattern.component_to_obs = component_to_obs
        pointpattern.dist_to_vertex = dist_to_vertex
        pointpattern.dist_snapped = dist_snapped
        pointpattern.obs_to_vertex = list(obs_to_vertex)

    def count_per_link(self, obs_on, graph=False):
        """Compute the counts per arc or edge (link).

        Parameters
        ----------
        obs_on : dict
            Dictionary of observations on the network.
            Either in the form ``{(<LINK>):{<POINT_ID>:(<COORDS>)}}``
            or ``{<LINK>:[(<COORD>),(<COORD>)]}``.
        graph : bool
            Count observations on graph edges (``True``) or
            network arcs (``False``). Default is ``False``.

        Returns
        -------
        counts : dict
            Counts per network link in the form ``{(<LINK>):<COUNT>}``.

        Examples
        --------

        Note that this passes the ``obs_to_arc`` or ``obs_to_edge`` attribute
        of a point pattern snapped to the network.

        >>> import spaghetti
        >>> from libpysal import examples
        >>> ntw = spaghetti.Network(examples.get_path("streets.shp"))

        Snap observations to the network.

        >>> ntw.snapobservations(
        ...     examples.get_path("crimes.shp"), "crimes", attribute=True
        ... )

        >>> counts = ntw.count_per_link(
        ...     ntw.pointpatterns["crimes"].obs_to_arc, graph=False
        ... )
        >>> counts[(140, 142)]
        10

        >>> s = sum([v for v in list(counts.values())])
        >>> s
        287

        """

        # instantiate observation counts by link lookup
        counts = {}

        # graph-theoretic object of nodes and edges
        if graph:
            # iterate the links-to-observations lookup
            for key, observations in obs_on.items():
                # isolate observation count for the link
                cnt = len(observations)

                # extract link (edges) key
                if key in self.arcs_to_edges:
                    key = self.arcs_to_edges[key]

                # either add to current count or a dictionary
                # entry or create new dictionary entry
                try:
                    counts[key] += cnt
                except KeyError:
                    counts[key] = cnt

        # network object of arcs and vertices
        else:
            # simplified version of the above process
            for key in obs_on:
                counts[key] = len(obs_on[key])

        return counts

    def _newpoint_coords(self, arc, distance):
        """Used internally to compute new point coordinates during snapping."""

        # extract coordinates for vertex 1 of arc
        x1 = self.vertex_coords[arc[0]][0]
        y1 = self.vertex_coords[arc[0]][1]

        # extract coordinates for vertex 2 of arc
        x2 = self.vertex_coords[arc[1]][0]
        y2 = self.vertex_coords[arc[1]][1]

        # if the network arc is vertical set the (x) coordinate
        # and proceed to calculating the (y) coordinate
        if x1 == x2:
            x0 = x1

            # if the vertical direction is positive from
            # vertex 1 to vertex 2 on the euclidean plane
            if y1 < y2:
                y0 = y1 + distance

            # if the vertical direction is negative from
            # vertex 1 to vertex 2 on the euclidean plane
            # -- this shouldn't happen due to vertex sorting in
            # -- self._extractnetwork() and self.extractgraph()
            elif y1 > y2:
                y0 = y2 + distance

            # otherwise the link is zero-length
            # -- this should never happen
            else:
                y0 = y1

            return x0, y0

        # calculate the slope of the arc, `m`
        m = (y2 - y1) / (x2 - x1)

        # if the horizontal direction is negative from
        # vertex 1 to vertex 2 on the euclidean plane
        if x1 > x2:
            x0 = x1 - distance / numpy.sqrt(1 + m**2)

        # if the horizontal direction is positive from
        # vertex 1 to vertex 2 on the euclidean plane
        elif x1 < x2:
            x0 = x1 + distance / numpy.sqrt(1 + m**2)

        # calculate the (y) coordinate
        y0 = m * (x0 - x1) + y1

        # the new (x,y) coordinates for the snapped observation
        return x0, y0

    def simulate_observations(self, count, distribution="uniform"):
        """Generate a simulated point pattern on the network.

        Parameters
        ----------
        count : int
            The number of points to create.
        distribution : str
            A distribution of random points. Currently, the only
            supported distribution is uniform.

        Returns
        -------
        random_pts : dict
            Keys are the edge tuple. Values are lists of new point coordinates.

        See also
        --------

        numpy.random.Generator.uniform

        Examples
        --------

        Instantiate a network.

        >>> import spaghetti
        >>> from libpysal import examples
        >>> ntw = spaghetti.Network(examples.get_path("streets.shp"))

        Snap observations to the network.

        >>> ntw.snapobservations(
        ...     examples.get_path("crimes.shp"), "crimes", attribute=True
        ... )

        Isolate the number of points in the dataset.

        >>> npts = ntw.pointpatterns["crimes"].npoints
        >>> npts
        287

        Simulate ``npts`` number of points along the network
        in a `uniform` distribution.

        >>> sim = ntw.simulate_observations(npts)
        >>> isinstance(sim, spaghetti.network.SimulatedPointPattern)
        True
        >>> sim.npoints
        287

        """

        # instantiate an empty `SimulatedPointPattern()`
        simpts = SimulatedPointPattern()

        # record throw-away arcs enumerator
        arcs_ = []

        # create array and fill each entry as length of network arc
        lengths = numpy.zeros(len(self.arc_lengths))
        for i, key in enumerate(self.arc_lengths.keys()):
            arcs_.append(key)
            lengths[i] = self.arc_lengths[key]

        # cumulative network length
        stops = numpy.cumsum(lengths)
        cumlen = stops[-1]

        # create lengths with a uniform distribution
        if distribution.lower() == "uniform":
            nrandompts = numpy.random.uniform(0, cumlen, size=(count,))
        else:
            raise RuntimeError(f"{distribution} distribution not currently supported.")

        # iterate over random distances created above
        for i, r in enumerate(nrandompts):
            # take the first element of the index array (arc ID) where the
            # random distance is greater than that of its value in `stops`
            idx = numpy.where(r < stops)[0][0]

            # assign the simulated point to the arc
            assignment_arc = arcs_[idx]

            # calculate and set the distance from the arc start
            distance_from_start = stops[idx] - r

            # populate the coordinates dict
            x0, y0 = self._newpoint_coords(assignment_arc, distance_from_start)

            # record the snapped coordinates and associated vertices
            simpts.snapped_coordinates[i] = (x0, y0)
            simpts.obs_to_vertex[assignment_arc[0]].append(i)
            simpts.obs_to_vertex[assignment_arc[1]].append(i)

            # calculate and set the distance from the arc end
            distance_from_end = self.arc_lengths[arcs_[idx]] - distance_from_start

            # populate the distances to vertices
            simpts.dist_to_vertex[i] = {
                assignment_arc[0]: distance_from_start,
                assignment_arc[1]: distance_from_end,
            }

            # set snapped coordinates and point count attributes
            simpts.points = simpts.snapped_coordinates
            simpts.npoints = len(simpts.points)

        return simpts

    def enum_links_vertex(self, v0):
        """Returns the arcs (links) adjacent to vertices.

        Parameters
        ----------
        v0 : int
            The vertex ID.

        Returns
        -------
        links : list
            List of tuple arcs adjacent to the vertex.

        Examples
        --------

        Create an instance of a network.

        >>> import spaghetti
        >>> from libpysal import examples
        >>> ntw = spaghetti.Network(examples.get_path("streets.shp"))

        Enumerate the links/arcs that are adjacent to vertex ``24``.

        >>> ntw.enum_links_vertex(24)
        [(24, 48), (24, 25), (24, 26)]

        """

        # instantiate links list
        links = []

        neighbor_vertices = self.adjacencylist[v0]

        # enumerate links associated with the current vertex
        for n in neighbor_vertices:
            links.append(tuple(sorted([n, v0])))

        return links

    def full_distance_matrix(self, n_processes, gen_tree=False):
        """All vertex-to-vertex distances on a network. This method
        is called from within ``allneighbordistances()``,
        ``nearestneighbordistances()``, and ``distancebandweights()``.

        Parameters
        ----------
        n_processes : int
            Specify the number of cores to utilize. Default is 1 core.
            Use ``"all"`` to request all available cores.
            Specify the exact number of cores with an integer.
        gen_tree : bool
            Rebuild shortest path ``True``, or skip ``False``.
            Default is ``False``.

        Notes
        -----

        Based on :cite:`Dijkstra1959a` and :cite:`doi:10.1002/9781119967101.ch3`.

        """

        # create an empty matrix which will store shortest path distance
        nvtx = len(self.vertex_list)
        self.distance_matrix = numpy.empty((nvtx, nvtx))

        # create `network_trees` attribute that stores
        # all network path trees (if desired)
        self.network_trees = {}

        # single-core processing
        if n_processes == 1:
            # iterate over each network vertex
            for vtx in self.vertex_list:
                # calculate the shortest path and preceding
                # vertices for traversal route
                distance, pred = util.dijkstra(self, vtx)
                pred = numpy.array(pred)

                # generate the shortest path tree
                tree = util.generatetree(pred) if gen_tree else None

                # populate distances and paths
                self.distance_matrix[vtx] = distance
                self.network_trees[vtx] = tree

        # multiprocessing
        else:
            # set up multiprocessing schema
            import multiprocessing as mp
            from itertools import repeat

            cores = mp.cpu_count() if n_processes == "all" else n_processes

            with mp.Pool(processes=cores) as p:
                # calculate the shortest path and preceding
                # vertices for traversal route by mapping each process
                distance_pred = p.map(
                    util.dijkstra_mp, zip(repeat(self), self.vertex_list)
                )

                # set range of iterations
                iterations = range(len(distance_pred))

                # fill shortest paths
                distance = [distance_pred[itr][0] for itr in iterations]

                # fill preceding vertices
                pred = numpy.array([distance_pred[itr][1] for itr in iterations])

                # iterate of network vertices and generate
                # the shortest path tree for each
                for vtx in self.vertex_list:
                    tree = util.generatetree(pred[vtx]) if gen_tree else None

                    # populate distances and paths
                    self.distance_matrix[vtx] = distance[vtx]
                    self.network_trees[vtx] = tree

    def allneighbordistances(
        self,
        sourcepattern,
        destpattern=None,
        fill_diagonal=None,
        n_processes=1,
        gen_tree=False,
        snap_dist=False,
    ):
        """Compute either all distances between :math:`i` and :math:`j` in a
        single point pattern or all distances between each :math:`i` from a
        source pattern and all :math:`j` from a destination pattern.

        Parameters
        ----------
        sourcepattern : {str, spaghetti.PointPattern}
            The key of a point pattern snapped to the network or
            the full ``spaghetti.PointPattern`` object.
        destpattern : str
            (Optional) The key of a point pattern snapped to the network
            or the full ``spaghetti.PointPattern`` object.
        fill_diagonal : {float, int}
            (Optional) Fill the diagonal of the cost matrix. Default is
            ``None`` and will populate the diagonal with ``numpy.nan``.
            Do not declare a ``destpattern`` for a custom
            ``fill_diagonal``.
        n_processes : {int, str}
            Specify the number of cores to utilize. Default is 1 core.
            Use ``"all"`` to request all available cores.
            Specify the exact number of cores with an integer.
        gen_tree : bool
            Rebuild shortest path ``True``, or skip ``False``.
            Default is ``False``.
        snap_dist : bool
            Flag as ``True`` to include the distance from the original
            location to the snapped location along the network. Default
            is ``False``.

        Returns
        -------
        nearest : numpy.ndarray
            An array of shape (n,m) storing distances between all
            source and destination points.
        tree_nearest : dict
            Nearest network node to point pattern vertex shortest
            path lookup. The values of the dictionary are a tuple
            of the nearest source vertex and the nearest destination
            vertex to query the lookup tree. If two observations are
            snapped to the same network arc a flag of -.1 is set for
            both the source and destination network vertex
            indicating the same arc is used while also raising an
            ``IndexError`` when rebuilding the path.

        Examples
        --------

        Create a network instance.

        >>> import spaghetti
        >>> from libpysal import examples
        >>> import numpy
        >>> ntw = spaghetti.Network(examples.get_path("streets.shp"))

        Snap observations to the network.

        >>> ntw.snapobservations(
        ...     examples.get_path("crimes.shp"), "crimes", attribute=True
        ... )

        Calculate all distances between observations in the ``crimes`` dataset.

        >>> s2s_dist = ntw.allneighbordistances("crimes")

        If calculating a ``type-a`` to ``type-a`` distance matrix
        the distance between an observation and itself is ``nan`` and
        the distance between one observation and another will be positive value.

        >>> s2s_dist[0,0], s2s_dist[1,0]
        (nan, 3105.189475447081)

        If calculating a ``type-a`` to ``type-b`` distance matrix
        the distance between all observations will likely be positive
        values, may be zero (or approximately zero), but will never be negative.

        >>> ntw.snapobservations(
        ...     examples.get_path("schools.shp"), "schools", attribute=False
        ... )
        >>> s2d_dist = ntw.allneighbordistances("crimes", destpattern="schools")
        >>> numpy.round((s2d_dist[0,0], s2d_dist[1,0]), 5)
        array([4520.72354, 6340.42297])


        Shortest paths can also be reconstructed when desired by
        setting the ``gen_tree`` keyword argument to ``True``. Here
        it is shown that the shortest path between school ``6`` and
        school ``7`` flows along network arcs through network
        vertices ``173`` and ``64``. The ``ntw.network_trees`` attribute
        may then be queried for the network elements comprising that path.

        >>> d2d_dist, tree = ntw.allneighbordistances("schools", gen_tree=True)
        >>> tree[(6, 7)]
        (173, 64)

        """

        # calculate the network vertex to vertex distance matrix
        # if it is not already an attribute
        if not hasattr(self, "distance_matrix"):
            self.full_distance_matrix(n_processes, gen_tree=gen_tree)

        # set the source and destination observation point patterns
        if type(sourcepattern) is str:
            sourcepattern = self.pointpatterns[sourcepattern]
            if destpattern:
                destpattern = self.pointpatterns[destpattern]

        # source pattern setup
        # set local copy of source pattern index
        src_indices = list(sourcepattern.points.keys())
        # set local copy of source distance to vertex lookup
        src_d2v = copy.deepcopy(sourcepattern.dist_to_vertex)
        # source point count
        nsource_pts = len(src_indices)
        # create source point to network vertex lookup
        src_vertices = {}
        for s in src_indices:
            v1, v2 = src_d2v[s].keys()
            src_vertices[s] = (v1, v2)

        # destination pattern setup
        # if only a source pattern is specified, also set it as
        # the destination pattern
        symmetric = False
        if destpattern is None:
            symmetric = True
            destpattern = sourcepattern
        # set local copy of destination pattern index
        dest_indices = list(destpattern.points.keys())
        # set local copy of destination distance to vertex lookup
        dst_d2v = copy.deepcopy(destpattern.dist_to_vertex)
        # destination point count
        ndest_pts = len(dest_indices)
        # create `deepcopy` of destination points to
        # consider for searching
        dest_searchpts = copy.deepcopy(dest_indices)
        # create destination point to network vertex lookup
        dest_vertices = {}
        for s in dest_indices:
            v1, v2 = dst_d2v[s].keys()
            dest_vertices[s] = (v1, v2)

        # add snapping distance to each pointpattern
        if snap_dist:
            # declare both point patterns and both
            # distance to vertex lookup in single lists
            patterns = [sourcepattern, destpattern]
            dist_copies = [src_d2v, dst_d2v]
            # iterate over each point pattern
            for elm, pp in enumerate(patterns):
                # extract associated vertex distances
                for pidx, dists_dict in dist_copies[elm].items():
                    # add snapped distance to each point
                    for vidx, vdist in dists_dict.items():
                        dists_dict[vidx] = vdist + pp.dist_snapped[pidx]

        # output setup
        # create empty source x destination array
        # and fill with infinity values
        nearest = numpy.empty((nsource_pts, ndest_pts))
        nearest[:] = numpy.inf
        # create empty dictionary to store path trees
        tree_nearest = {}

        # iterate over each point in sources
        for p1 in src_indices:
            # get the source vertices and dist to source vertices
            source1, source2 = src_vertices[p1]
            set1 = set(src_vertices[p1])

            # distance from source vertex1 to point and
            # distance from source vertex2 to point
            sdist1, sdist2 = src_d2v[p1].values()

            if symmetric:
                # only compute the upper triangle if symmetric
                dest_searchpts.remove(p1)

            # iterate over each point remaining in destinations
            for p2 in dest_searchpts:
                # get the destination vertices and
                # dist to destination vertices
                dest1, dest2 = dest_vertices[p2]
                set2 = set(dest_vertices[p2])

                # when the observations are snapped to the same arc
                if set1 == set2:
                    # calculate only the length between points along
                    # that arc
                    x1, y1 = sourcepattern.snapped_coordinates[p1]
                    x2, y2 = destpattern.snapped_coordinates[p2]

                    computed_length = util.compute_length((x1, y1), (x2, y2))
                    nearest[p1, p2] = computed_length

                    # set the nearest network vertices to a flag of -.1
                    # indicating the same arc is used while also raising
                    # and indexing error when rebuilding the path
                    tree_nearest[p1, p2] = SAME_SEGMENT

                # otherwise lookup distance between the source and
                # destination vertex
                else:
                    # distance from destination vertex1 to point and
                    # distance from destination vertex2 to point
                    ddist1, ddist2 = dst_d2v[p2].values()

                    # set the four possible combinations of
                    # source to destination shortest path traversal
                    d11 = self.distance_matrix[source1][dest1]
                    d21 = self.distance_matrix[source2][dest1]
                    d12 = self.distance_matrix[source1][dest2]
                    d22 = self.distance_matrix[source2][dest2]

                    # find the shortest distance from the path passing
                    # through each of the two origin vertices to the
                    # first destination vertex
                    sd_1 = d11 + sdist1
                    sd_21 = d21 + sdist2
                    sp_combo1 = source1, dest1
                    if sd_1 > sd_21:
                        sd_1 = sd_21
                        sp_combo1 = source2, dest1

                    # now add the point to vertex1 distance on
                    # the destination arc
                    len_1 = sd_1 + ddist1

                    # repeat the prior but now for the paths entering
                    # at the second vertex of the second arc
                    sd_2 = d12 + sdist1
                    sd_22 = d22 + sdist2
                    sp_combo2 = source1, dest2
                    if sd_2 > sd_22:
                        sd_2 = sd_22
                        sp_combo2 = source2, dest2
                    len_2 = sd_2 + ddist2

                    # now find the shortest distance path between point
                    # 1 on arc 1 and point 2 on arc 2, and assign
                    sp_12 = len_1
                    s_vertex, d_vertex = sp_combo1
                    if len_1 > len_2:
                        sp_12 = len_2
                        s_vertex, d_vertex = sp_combo2

                    # set distance and path tree
                    nearest[p1, p2] = sp_12
                    tree_nearest[p1, p2] = (s_vertex, d_vertex)

                if symmetric:
                    # mirror the upper and lower triangle
                    # when symmetric
                    nearest[p2, p1] = nearest[p1, p2]

        # populate the main diagonal when symmetric
        if symmetric:
            # fill the matrix diagonal with NaN values is no fill
            # value is specified
            if fill_diagonal is None:
                numpy.fill_diagonal(nearest, numpy.nan)

            # otherwise fill with specified value
            else:
                numpy.fill_diagonal(nearest, fill_diagonal)

        # if the nearest path tree is desired return it along
        # with the cost matrix
        if gen_tree:
            return nearest, tree_nearest

        else:
            return nearest

    def nearestneighbordistances(
        self,
        sourcepattern,
        destpattern=None,
        n_processes=1,
        gen_tree=False,
        all_dists=None,
        snap_dist=False,
        keep_zero_dist=True,
    ):
        """Compute the interpattern nearest neighbor distances or the
        intrapattern nearest neighbor distances between a source
        pattern and a destination pattern.

        Parameters
        ----------
        sourcepattern : str
            The key of a point pattern snapped to the network.
        destpattern : str
            (Optional) The key of a point pattern snapped to the
            network.
        n_processes : {int, str}
            Specify the number of cores to utilize. Default is 1 core.
            Use ``"all"`` to request all available cores.
            Specify the exact number of cores with an integer.
        gen_tree : bool
            Rebuild shortest path ``True``, or skip ``False``.
            Default is ``False``.
        all_dists : numpy.ndarray
            An array of shape :math:`(n,n)` storing distances between all
            points.
        snap_dist : bool
            Flag as ``True`` to include the distance from the original
            location to the snapped location along the network. Default
            is ``False``.
        keep_zero_dist : bool
            Include zero values in minimum distance ``True`` or exclude
            ``False``. Default is ``True``. If the source pattern is the
            same as the destination pattern the diagonal is filled with
            ``numpy.nan``.

        Returns
        -------
        nearest : dict
            Nearest neighbor distances keyed by the source point ID with
            the value as as tuple of lists containing
            nearest destination point ID(s) and distance.

        Examples
        --------

        Instantiate a network.

        >>> import spaghetti
        >>> from libpysal import examples
        >>> ntw = spaghetti.Network(examples.get_path("streets.shp"))

        Snap observations to the network.

        >>> ntw.snapobservations(examples.get_path("crimes.shp"), "crimes")

        Fetch nearest neighbor distances while (potentially)
        keeping neighbors that have been geocoded directly on top of
        each other. Here it is demonstrated that observation ``11``
        has two neighbors (``18`` and ``19``) at an exactly equal distance.
        However, observation ``18`` is shown to have only one neighbor
        (``18``) with no distance between them.

        >>> nn = ntw.nearestneighbordistances("crimes", keep_zero_dist=True)
        >>> nn[11], nn[18]
        (([18, 19], 165.33982412719126), ([19], 0.0))

        This may be remedied by setting the ``keep_zero_dist`` keyword
        argument to ``False``. With this parameter set, observation ``11``
        still has the same neighbor/distance values, but
        observation ``18`` now has a single nearest neighbor (``11``)
        with a non-zero, postive distance.

        >>> nn = ntw.nearestneighbordistances("crimes", keep_zero_dist=False)
        >>> nn[11], nn[18]
        (([18, 19], 165.33982412719126), ([11], 165.33982412719126))

        There are valid reasons for both retaining or masking zero distance
        neighbors. When conducting analysis, thought must be given as to
        which model more accurately represents the specific scenario.

        """

        # raise exception is the specified point pattern does not exist
        if sourcepattern not in self.pointpatterns.keys():
            raise KeyError(f"Available point patterns are {self.pointpatterns.keys()}")

        # calculate the network vertex to vertex distance matrix
        # if it is not already an attribute
        if not hasattr(self, "distance_matrix"):
            self.full_distance_matrix(n_processes, gen_tree=gen_tree)

        # determine if the source and destination patterns are equal
        symmetric = sourcepattern != destpattern

        # (for source-to-source patterns) if zero-distance neighbors are
        # desired, keep the diagonal as NaN and take the minimum
        # distance neighbor(s), which may include zero distance
        # neighors.
        fill_diagonal = None
        if not keep_zero_dist and symmetric:
            # (for source-to-source patterns) if zero-distance neighbors
            # should be ignored, convert the diagonal to 0.0 and take
            # the minimum distance neighbor(s) that is/are not 0.0
            # distance.
            fill_diagonal = 0.0

        # set the source and destination observation point patterns
        sourcepattern = self.pointpatterns[sourcepattern]
        if destpattern:
            destpattern = self.pointpatterns[destpattern]

        # if the full source to destination is not calculated,
        # do that at this time
        if all_dists is None:
            all_dists = self.allneighbordistances(
                sourcepattern,
                destpattern=destpattern,
                fill_diagonal=fill_diagonal,
                n_processes=n_processes,
                gen_tree=gen_tree,
                snap_dist=snap_dist,
            )

        # create empty nearest neighbors lookup
        nearest = {}

        # iterate over each source point
        for source_index in sourcepattern.points:
            # this considers all zero-distance neighbors
            if keep_zero_dist and symmetric:
                val = numpy.nanmin(all_dists[source_index, :])

            # this does not consider zero-distance neighbors
            else:
                val = numpy.min(
                    all_dists[source_index, :][
                        numpy.nonzero(all_dists[source_index, :])
                    ]
                )

            # nearest destination (may be more than one if
            # observations are equal distances away)
            dest_idxs = numpy.where(all_dists[source_index, :] == val)[0].tolist()

            # set nearest destination point(s) and distance
            nearest[source_index] = (dest_idxs, val)

        return nearest

    def shortest_paths(self, tree, pp_orig, pp_dest=None):
        """Return the shortest paths between observation points as
        ``libpysal.cg.Chain`` objects.

        Parameters
        ----------
        tree : dict
            See ``tree_nearest`` in
            ``spaghetti.Network.allneighbordistances()``.
        pp_orig : str
            Origin point pattern for shortest paths.
            See ``name`` in ``spaghetti.Network.snapobservations()``.
        pp_dest : str
            Destination point pattern for shortest paths.
            See ``name`` in ``spaghetti.Network.snapobservations()``.
            Defaults ``pp_orig`` if not declared.

        Returns
        -------
        paths : list
            The shortest paths between observations as geometric objects.
            Each element of the list is a list where the first element
            is an origin-destination pair tuple and the second
            element is a ``libpysal.cg.Chain``.

        Raises
        ------
        AttributeError
            This exception is raised when an attempt to extract shortest
            path geometries is being made that but the ``network_trees``
            attribute does not exist within the network object.

        Examples
        --------

        Instantiate a network.

        >>> import spaghetti
        >>> from libpysal import examples
        >>> ntw = spaghetti.Network(examples.get_path("streets.shp"))

        Snap observations to the network.

        >>> ntw.snapobservations(examples.get_path("schools.shp"), "schools")

        Create shortest path trees between observations.

        >>> _, tree = ntw.allneighbordistances("schools", gen_tree=True)

        Generate geometric objects from trees.

        >>> paths = ntw.shortest_paths(tree, "schools")

        Extract the first path, which is between observations
        ``0`` and ``1``.

        >>> path = paths[0]
        >>> path[0]
        (0, 1)

        The are ``n`` vertices in the path between observations
        ``0`` and ``1``.

        >>> n = len(path[1].vertices)
        >>> n
        10

        """

        # build the network trees object if it is not already an attribute
        if not hasattr(self, "network_trees"):
            msg = (
                "The 'network_trees' attribute has not been created. "
                "Rerun 'spaghetti.Network.allneighbordistances()' "
                "with the 'gen_tree' parameter set to 'True'."
            )
            raise AttributeError(msg)

        # isolate network attributes
        pp_orig = self.pointpatterns[pp_orig]
        pp_dest = self.pointpatterns[pp_dest] if pp_dest else pp_orig
        vtx_coords = self.vertex_coords
        net_trees = self.network_trees

        # instantiate a list to store paths
        paths = []

        # iterate over each path in the tree
        for idx, ((obs0, obs1), (v0, v1)) in enumerate(tree.items()):
            # if the observations share the same segment
            # create a partial segment path
            if (v0, v1) == SAME_SEGMENT:
                # isolate the snapped coordinates and put in a list
                partial_segment_verts = [
                    cg.Point(pp_orig.snapped_coordinates[obs0]),
                    cg.Point(pp_dest.snapped_coordinates[obs1]),
                ]
                path = partial_segment_verts

            else:
                # source and destination network vertices
                svtx, dvtx = tree[obs0, obs1]

                # path passes through these nodes
                # (source and destination inclusive)
                thru_nodes = net_trees[svtx][dvtx][::-1] + [dvtx]

                # full-length network segments along path
                full_segs_path = []
                iter_limit = len(thru_nodes) - 1
                for _idx, item in enumerate(islice(thru_nodes, iter_limit)):
                    full_segs_path.append((item, thru_nodes[_idx + 1]))

                # create copy of arc paths dataframe
                full_segments = []
                for fsp in full_segs_path:
                    full_segments.append(util._chain_constr(vtx_coords, fsp))

                # unpack the vertices containers
                segm_verts = [v for fs in full_segments for v in fs.vertices]

                # remove duplicate vertices
                for idx, v in enumerate(segm_verts):
                    try:
                        if v == segm_verts[idx + 1]:
                            segm_verts.remove(v)
                    except IndexError as e:
                        if e.args[0] == "list index out of range":
                            continue
                        else:
                            raise

                # partial-length network segments along path
                partial_segment_verts = [
                    cg.Point(pp_orig.snapped_coordinates[obs0]),
                    cg.Point(pp_dest.snapped_coordinates[obs1]),
                ]

                # combine the full and partial segments into a single list
                first_vtx, last_vtx = partial_segment_verts
                path = [first_vtx] + segm_verts + [last_vtx]

            # populate the ``paths`` dataframe
            paths.append([(obs0, obs1), util._chain_constr(None, path)])

        return paths

    def split_arcs(
        self, split_param, split_by="distance", w_components=True, weights_kws=dict()
    ):
        """Split all network arcs at either a fixed distance or fixed count.

        Parameters
        ----------
        split_param : {int, float}
            Either the number of desired resultant split arcs or
            the distance at which arcs are split.
        split_by : str
            Either ``'distance'`` or ``'count'``. Default is ``'distance'``.
        w_components : bool
            Set to ``False`` to not record connected components from a
            ``libpysal.weights.W`` object. Default is ``True``.
        weights_kws : dict
            Keyword arguments for ``libpysal.weights.W``. Default is ``dict()``.

        Returns
        -------
        split_network : spaghetti.Network
            A newly instantiated ``spaghetti.Network`` object.

        Examples
        --------

        Instantiate a network.

        >>> import spaghetti
        >>> from libpysal import examples
        >>> ntw = spaghetti.Network(examples.get_path("streets.shp"))

        Split the network into a segments of 200 distance units in length
        (US feet in this case).
        This will include "remainder" segments unless the network is
        comprised of arcs with lengths exactly divisible by ``distance``.

        >>> n200 = ntw.split_arcs(200.0)
        >>> len(n200.arcs)
        688

        The number of arcs within the new object can be accessed via the
        weights object, as well. These counts will be equal.

        >>> len(n200.arcs) == n200.w_network.n
        True

        Neighboring arcs can also be queried through the weight object.

        >>> n200.w_network.neighbors[72,392]
        [(71, 72), (72, 252), (72, 391), (392, 393)]

        Network arcs can also be split by a specified number of divisions with
        the ``split_by`` keyword set to ``'count'``, which is ``'distance'`` by
        default. For example, each arc can be split into 2 equal parts.

        >>> n2 = ntw.split_arcs(2, split_by="count")
        >>> len(n2.arcs)
        606

        """

        def int_coord(c):
            """convert coordinates for integers if possible
            e.g., (1.0, 0.5) --> (1, 0.5)
            """
            return int(c) if (type(c) == float and c.is_integer()) else c

        # catch invalid split types
        _split_by = split_by.lower()
        valid_split_types = ["distance", "count"]
        if _split_by not in valid_split_types:
            msg = (
                f"'{split_by}' is not a valid value for 'split_by'. "
                f"Valid arguments include: {valid_split_types}."
            )
            raise ValueError(msg)

        # catch invalid count params
        if _split_by == "count":
            if split_param <= 1:
                msg = (
                    "Splitting arcs by 1 or less is not possible. "
                    f"Currently 'split_param' is set to {split_param}."
                )
                raise ValueError(msg)
            split_integer = int(split_param)
            if split_param != split_integer:
                msg = (
                    "Network arcs must split by an integer. "
                    f"Currently 'split_param' is set to {split_param}."
                )
                raise TypeError(msg)

        # create new shell network instance
        split_network = Network()

        # duplicate input network attributes
        split_network.adjacencylist = copy.deepcopy(self.adjacencylist)
        split_network.arc_lengths = copy.deepcopy(self.arc_lengths)
        split_network.arcs = copy.deepcopy(self.arcs)
        split_network.vertex_coords = copy.deepcopy(self.vertex_coords)
        split_network.vertex_list = copy.deepcopy(self.vertex_list)
        split_network.vertices = copy.deepcopy(self.vertices)
        split_network.pointpatterns = copy.deepcopy(self.pointpatterns)
        split_network.in_data = self.in_data

        # set vertex ID to start iterations
        current_vertex_id = max(self.vertices.values())

        # instantiate sets for newly created network arcs and
        # input network arcs to remove
        new_arcs = set()
        remove_arcs = set()

        # iterate over all network arcs
        for arc in split_network.arcs:
            # fetch network arc length
            length = split_network.arc_lengths[arc]

            # set initial segmentation interval
            if _split_by == "distance":
                interval = split_param
            else:
                interval = length / float(split_param)

            # initialize arc new arc length at zero
            totallength = 0

            # initialize the current vertex and ending vertex
            currentstart, end_vertex = arc[0], arc[1]

            # determine direction of arc vertices
            csx, csy = split_network.vertex_coords[currentstart]
            evx, evy = split_network.vertex_coords[end_vertex]
            if csy > evy and csx == evx:
                currentstart, end_vertex = end_vertex, currentstart

            # if the arc will be split remove the current
            # arc from the adjacency list
            if interval < length:
                # remove old arc adjacency information
                split_network.adjacencylist[currentstart].remove(end_vertex)
                split_network.adjacencylist[end_vertex].remove(currentstart)

                # remove old arc length information
                split_network.arc_lengths.pop(arc, None)

                # add old arc to set of arcs to remove
                remove_arcs.add(arc)

            # if the arc will not be split, do nothing and continue
            else:
                continue

            # traverse the length of the arc
            while totallength < length:
                # once an arc can not be split further
                if totallength + interval >= length:
                    # record the ending vertex
                    currentstop = end_vertex
                    # set the length remainder
                    interval = length - totallength
                    # full old length reached
                    totallength = length

                else:
                    # set the current vertex ID
                    current_vertex_id += 1
                    # set the current stopping ID
                    currentstop = current_vertex_id
                    # add the interval distance to the traversed length
                    totallength += interval

                    # compute the new vertex coordinate
                    newx, newy = self._newpoint_coords(arc, totallength)
                    new_vertex = (int_coord(newx), int_coord(newy))

                    # update the vertex and coordinate info if needed
                    if new_vertex not in split_network.vertices.keys():
                        split_network.vertices[new_vertex] = currentstop
                        split_network.vertex_coords[currentstop] = new_vertex
                        split_network.vertex_list.append(currentstop)
                    else:
                        # retrieve vertex ID if coordinate already exists
                        current_vertex_id -= 1
                        currentstop = split_network.vertices[new_vertex]

                # update the new network adjacency list
                split_network.adjacencylist[currentstart].append(currentstop)
                split_network.adjacencylist[currentstop].append(currentstart)

                # add the new arc to the arc dictionary
                # iterating over this so we need to add after iterating
                _new_arc = tuple(sorted([currentstart, currentstop]))
                new_arcs.add(_new_arc)

                # set the length of the arc
                split_network.arc_lengths[_new_arc] = interval

                # increment the starting vertex to the stopping vertex
                currentstart = currentstop

        # add the newly created arcs to the network and remove the old arcs
        split_network.arcs = set(split_network.arcs)
        split_network.arcs.update(new_arcs)
        split_network.arcs.difference_update(remove_arcs)
        split_network.arcs = sorted(split_network.arcs)

        # extract connected components
        if w_components:
            # extract contiguity weights from libpysal
            split_network.w_network = split_network.contiguityweights(
                graph=False, from_split=True, weights_kws=weights_kws
            )
            # identify connected components from the `w_network`
            split_network.identify_components(split_network.w_network, graph=False)

        # update the snapped point pattern
        for instance in split_network.pointpatterns.values():
            split_network._snap_to_link(instance)

        return split_network

    def GlobalAutoK(  # noqa N802
        self,
        pointpattern,
        nsteps=10,
        permutations=99,
        threshold=0.5,
        distribution="uniform",
        upperbound=None,
    ):
        r"""Compute a global auto :math:`K`-function based on a network constrained
        cost matrix through
        `Monte Carlo simulation <https://en.wikipedia.org/wiki/Monte_Carlo_method>`_
        according to the formulation adapted from
        :cite:`doi:10.1002/9780470549094.ch5`. See the **Notes**
        section for further description.

        Parameters
        ----------
        pointpattern : spaghetti.PointPattern
            A ``spaghetti`` point pattern object.
        nsteps : int
            The number of steps at which the count of the nearest
            neighbors is computed. Default is ``10``.
        permutations : int
            The number of permutations to perform. Default is ``99``.
        threshold : float
            The level at which significance is computed.
            (0.5 would be 97.5% and 2.5%). Default is ``0.5``.
        distribution : str
            The distribution from which random points are sampled.
            Currently, the only supported distribution is ``'uniform'``.
        upperbound : float
            The upper bound at which the :math:`K`-function is computed.
            Defaults to the maximum observed nearest neighbor distance.

        Returns
        -------
        GlobalAutoK : spaghetti.analysis.GlobalAutoK
            The global auto :math:`K`-function class instance.

        Notes
        -----

        The :math:`K`-function can be formulated as:

        .. math::

           \displaystyle K(r)=\frac{\sum^n_{i=1} \#[\hat{A} \in D(a_i, r)]}{n\lambda},

        where $n$ is the set cardinality of :math:`A`, :math:`\hat{A}` is the
        subset of observations in :math:`A` that are within :math:`D` units of
        distance from :math:`a_i` (each single observation in :math:`A`), and :math:`r`
        is the range of distance values over which the :math:`K`-function is
        calculated. The :math:`\lambda` term is the intensity of observations
        along the network, calculated as:

        .. math::

           \displaystyle \lambda = \frac{n}{\big|N_{arcs}\big|},

        where :math:`\big|N_{arcs}\big|` is the summed length of network arcs.
        The global auto :math:`K`-function measures overall clustering in one set of
        observations by comparing all intra-set distances over a range of
        distance buffers :math:`D \in r`. The :math:`K`-function improves upon
        nearest-neighbor distance measures through the analysis of all neighbor
        distances. For an explanation on how to interpret the results of the
        :math:`K`-function see the Network Spatial Dependence tutorial
        `here <https://pysal.org/spaghetti/notebooks/network-spatial-dependence.html>`_.

        For original implementation see :cite:`Ripley1976`
        and :cite:`Ripley1977`.
        For further Network-`K` formulations see
        :cite:`doi:10.1111/j.1538-4632.2001.tb00448.x`,
        :cite:`doi:10.1002/9781119967101.ch6`, and
        :cite:`Baddeley2020`.

        See also
        --------

        pointpats.K

        Examples
        --------

        Create a network instance.

        >>> import spaghetti
        >>> from libpysal import examples
        >>> ntw = spaghetti.Network(in_data=examples.get_path("streets.shp"))

        Snap observation points onto the network.

        >>> pt_str = "schools"
        >>> in_data = examples.get_path(pt_str+".shp")
        >>> ntw.snapobservations(in_data, pt_str, attribute=True)
        >>> schools = ntw.pointpatterns[pt_str]

        Compute a :math:`K`-function from school observations
        with ``99`` ``permutations`` at ``10`` intervals.

        >>> kres = ntw.GlobalAutoK(schools, permutations=99, nsteps=10)
        >>> kres.lowerenvelope.shape[0]
        10

        """

        # call analysis.GlobalAutoK
        return GlobalAutoK(
            self,
            pointpattern,
            nsteps=nsteps,
            permutations=permutations,
            threshold=threshold,
            distribution=distribution,
            upperbound=upperbound,
        )

    def Moran(self, pp_name, permutations=999, graph=False):  # noqa N802
        """Calculate a Moran's *I* statistic on a set of observations
        based on network arcs. The Moran’s *I* test statistic allows
        for the inference of how clustered (or dispersed) a dataset is
        while considering both attribute values and spatial relationships.
        A value of closer to +1 indicates absolute clustering while a
        value of closer to -1 indicates absolute dispersion. Complete
        spatial randomness takes the value of 0. See the ``esda``
        `documentation <https://pysal.org/esda/generated/esda.Moran.html#esda.Moran>`_
        for in-depth descriptions and tutorials.

        Parameters
        ----------
        pp_name : str
            The name of the point pattern in question.
        permutations : int
            The number of permutations to perform. Default is ``999``.
        graph : bool
            Perform the Moran calculation on the graph `W` object
            (``True``). Default is ``False``, which performs the
            Moran calculation on the network `W` object.

        Returns
        -------
        moran : esda.Moran
            A Moran's *I* statistic object results.
        y : list
            The y-axis (counts).

        Examples
        --------

        Create a network instance.

        >>> import spaghetti
        >>> from libpysal import examples
        >>> ntw = spaghetti.Network(in_data=examples.get_path("streets.shp"))

        Snap observation points onto the network.

        >>> crimes = "crimes"
        >>> in_data = examples.get_path(crimes+".shp")
        >>> ntw.snapobservations(in_data, crimes, attribute=True)

        Compute a Moran's :math:`I` from crime observations.

        >>> moran_res, _ = ntw.Moran(crimes)
        >>> round(moran_res.I, 6)
        0.005193

        Notes
        -----

        See :cite:`moran:_cliff81` and :cite:`esda:_2019` for more details.

        """

        # set proper weights attribute
        w = self.w_graph if graph else self.w_network

        # Compute the counts
        pointpat = self.pointpatterns[pp_name]
        counts = self.count_per_link(pointpat.obs_to_arc, graph=graph)

        # Build the y vector
        y = [counts[i] if i in counts else 0.0 for i in w.neighbors]

        # Moran's I
        moran = esda.moran.Moran(y, w, permutations=permutations)

        return moran, y

    def savenetwork(self, filename):
        """Save a network to disk as a binary file.

        Parameters
        ----------
        filename : str
            The filename where the network should be saved. This should
            be a full path or it will be saved in the current directory.

        Examples
        --------

        Create a network instance.

        >>> import spaghetti
        >>> from libpysal import examples
        >>> ntw = spaghetti.Network(examples.get_path("streets.shp"))

        Save out the network instance.

        >>> ntw.savenetwork("mynetwork.pkl")

        """

        with _open(filename, "wb") as networkout:
            pickle.dump(self, networkout, protocol=2)

    @staticmethod
    def loadnetwork(filename):
        """Load a network from a binary file saved on disk.

        Parameters
        ----------
        filename : str
            The filename where the network is saved.

        Returns
        -------
        self : spaghetti.Network
            A pre-computed ``spaghetti`` network object.

        """

        with _open(filename, "rb") as networkin:
            self = pickle.load(networkin)

        return self


def extract_component(net, component_id, weightings=None):
    """Extract a single component from a network object.

    Parameters
    ----------
    net : spaghetti.Network
        Full network object.
    component_id : int
        The ID of the desired network component.
    weightings : {dict, bool}
        See the ``weightings`` keyword argument in ``spaghetti.Network``.

    Returns
    -------
    cnet : spaghetti.Network
        The pruned network containing the component specified in
        ``component_id``.

    Notes
    -----

    Point patterns are not reassigned when extracting a component. Therefore,
    component extraction should be performed prior to snapping any point
    sets onto the network. Also, if the ``spaghetti.Network`` object
    has ``distance_matrix`` or ``network_trees`` attributes, they are
    deleted and must be computed again on the single component.

    Examples
    --------

    Instantiate a network object.

    >>> from libpysal import examples
    >>> import spaghetti
    >>> snow_net = examples.get_path("Soho_Network.shp")
    >>> ntw = spaghetti.Network(
    ...     in_data=snow_net,
    ...     extractgraph=False,
    ...     weights_kws=dict(silence_warnings=True)
    ... )

    The network is not fully connected.

    >>> ntw.network_fully_connected
    False

    Examine the number of network components.

    >>> ntw.network_n_components
    45

    Extract the longest component.

    >>> longest = spaghetti.extract_component(ntw, ntw.network_longest_component)
    >>> longest.network_n_components
    1
    >>> longest.network_component_lengths
    {0: 13508.169276875526}

    """

    def _reassign(attr, cid):
        """Helper for reassigning attributes."""

        # set for each attribute(s)
        if attr == "_fully_connected":
            _val = [True for objt in obj_type]
            attr = [objt + attr for objt in obj_type]
        elif attr == "_n_components":
            _val = [1 for objt in obj_type]
            attr = [objt + attr for objt in obj_type]
        elif attr in ["_longest_component", "_largest_component"]:
            _val = [cid for objt in obj_type]
            attr = [objt + attr for objt in obj_type]
        elif attr == "vertex_list":
            # reassigns vertex list + network, graph component vertices
            supp = [objt + "_component_vertices" for objt in obj_type]
            _val = [getattr(cnet, supp[0])[cid]]
            _val += [{cid: getattr(cnet, s)[cid]} for s in supp]
            attr = [attr] + supp
        elif attr == "vertex_coords":
            # reassigns both vertex_coords and vertices
            supp = getattr(cnet, "vertex_list")
            _val = [{k: v for k, v in getattr(cnet, attr).items() if k in supp}]
            _val += [{v: k for k, v in _val[0].items()}]
            attr = [attr, "vertices"]
        elif attr == "_component_vertex_count":
            # reassigns both network and graph _component_vertex_count
            supp = len(getattr(cnet, "vertex_list"))
            _val = [{cid: supp} for objt in obj_type]
            attr = [objt + attr for objt in obj_type]
        elif attr == "adjacencylist":
            supp_adj = copy.deepcopy(list(getattr(cnet, attr).keys()))
            supp_vtx = getattr(cnet, "vertex_list")
            supp_rmv = [v for v in supp_adj if v not in supp_vtx]
            [getattr(cnet, attr).pop(s) for s in supp_rmv]
            return
        elif attr == "_component_is_ring":
            # reassigns both network and graph _component_is_ring
            supp = [getattr(cnet, objt + attr) for objt in obj_type]
            _val = [{cid: s[cid]} for s in supp]
            attr = [objt + attr for objt in obj_type]
        elif attr == "non_articulation_points":
            supp_vtx = getattr(cnet, "vertex_list")
            _val = [[s for s in getattr(cnet, attr) if s in supp_vtx]]
            attr = [attr]
        elif attr == "_component2":
            # reassigns both network and graph _component2 attributes
            supp = [_n + "_component2" + _a]
            if hasgraph:
                supp += [_g + "_component2" + _e]
            _val = [{cid: getattr(cnet, s)[cid]} for s in supp]
            attr = supp
        elif attr == "arcs":
            # reassigns both arcs and edges
            c2 = "_component2"
            supp = [_n + c2 + _a]
            if hasgraph:
                supp += [_g + c2 + _e]
            _val = [getattr(cnet, s)[cid] for s in supp]
            attr = [attr]
            if hasgraph:
                attr += ["edges"]
        elif attr == "_component_labels":
            # reassigns both network and graph _component_labels
            supp = [len(getattr(cnet, o + "s")) for o in obj]
            _val = [numpy.array([cid] * s) for s in supp]
            attr = [objt + attr for objt in obj_type]
        elif attr == "_component_lengths":
            # reassigns both network and graph _component_lengths
            supp = [objt + attr for objt in obj_type]
            _val = [{cid: getattr(cnet, s)[cid]} for s in supp]
            attr = supp
        elif attr == "_lengths":
            # reassigns both arc and edge _lengths
            supp_name = [o + attr for o in obj]
            supp_lens = [getattr(cnet, s) for s in supp_name]
            supp_link = [getattr(cnet, o + "s") for o in obj]
            supp_ll = list(zip(supp_lens, supp_link))
            _val = [{k: v for k, v in l1.items() if k in l2} for l1, l2 in supp_ll]
            attr = supp_name

        # reassign attributes
        for a, av in zip(attr, _val):
            setattr(cnet, a, av)

    # provide warning (for now) if the network contains a point pattern
    if getattr(net, "pointpatterns"):
        msg = (
            "There is a least one point pattern associated with the network."
            " Component extraction should be performed prior to snapping"
            " point patterns to the network object; failing to do so may"
            " lead to unexpected results."
        )
        warnings.warn(msg, stacklevel=2)
    # provide warning (for now) if the network contains a point pattern
    dm, nt = "distance_matrix", "network_trees"
    if hasattr(net, dm) or hasattr(net, nt):
        msg = (
            f"Either one or both ({dm}, {nt}) attributes"
            " are present and will be deleted. These must be"
            " recalculated following component extraction."
        )
        warnings.warn(msg, stacklevel=2)
        for attr in [dm, nt]:
            if hasattr(net, attr):
                delattr(net, attr)

    # make initial copy of the network
    cnet = copy.deepcopy(net)

    # set labels
    _n, _a, _g, _e = "network", "arc", "graph", "edge"
    obj_type = [_n]
    obj = [_a]
    hasgraph = False
    if hasattr(cnet, "w_graph"):
        obj_type += [_g]
        obj += [_e]
        hasgraph = True

    # attributes to reassign
    update_attributes = [
        "_fully_connected",
        "_n_components",
        "_longest_component",
        "_largest_component",
        "vertex_list",
        "vertex_coords",
        "_component_vertex_count",
        "adjacencylist",
        "_component_is_ring",
        "_component2",
        "arcs",
        "_component_lengths",
        "_lengths",
        "_component_labels",
    ]
    if hasgraph:
        update_attributes.append("non_articulation_points")

    # reassign attributes
    for attribute in update_attributes:
        _reassign(attribute, component_id)

    # recreate spatial weights
    cnet.w_network = cnet.contiguityweights(graph=False, weightings=weightings)
    if hasgraph:
        cnet.w_graph = cnet.contiguityweights(graph=True, weightings=weightings)

    return cnet


def spanning_tree(net, method="sort", maximum=False, silence_warnings=True):
    """Extract a minimum or maximum spanning tree from a network.

    Parameters
    ----------
    net : spaghetti.Network
        Instance of a network object.
    method : str
        Method for determining spanning tree. Currently, the only
        supported method is 'sort', which sorts the network arcs
        by length prior to building intermediary networks and checking
        for cycles within the tree/subtrees. Future methods may
        include linear programming approachs, etc.
    maximum : bool
        When ``True`` a maximum spanning tree is created. When ``False``
        a minimum spanning tree is created. Default is ``False``.
    silence_warnings : bool
        Warn if there is more than one connected component. Default is
        ``False`` due to the nature of constructing a minimum
        spanning tree.

    Returns
    -------
    net : spaghetti.Network
        Pruned instance of the network object.

    Notes
    -----

    For in-depth background and details see
    :cite:`GrahamHell_1985`,
    :cite:`AhujaRavindraK`, and
    :cite:`Okabe2012`.

    See also
    --------

    networkx.algorithms.tree.mst
    scipy.sparse.csgraph.minimum_spanning_tree

    Examples
    --------

    Create a network instance.

    >>> from libpysal import cg
    >>> import spaghetti
    >>> p00 = cg.Point((0,0))
    >>> lines = [cg.Chain([p00, cg.Point((0,3)), cg.Point((4,0)), p00])]
    >>> ntw = spaghetti.Network(in_data=lines)

    Extract the minimum spanning tree.

    >>> minst_net = spaghetti.spanning_tree(ntw)
    >>> min_len = sum(minst_net.arc_lengths.values())
    >>> min_len
    7.0

    Extract the maximum spanning tree.

    >>> maxst_net = spaghetti.spanning_tree(ntw, maximum=True)
    >>> max_len = sum(maxst_net.arc_lengths.values())
    >>> max_len
    9.0

    >>> max_len > min_len
    True

    """

    # (un)silence warning
    weights_kws = {"silence_warnings": silence_warnings}
    # do not extract graph object while testing for cycles
    net_kws = {"extractgraph": False, "weights_kws": weights_kws}

    # if the network has no cycles, it is already a spanning tree
    if util.network_has_cycle(net.adjacencylist):
        if method.lower() == "sort":
            spanning_tree = mst_weighted_sort(net, maximum, net_kws)
        else:
            msg = f"'{method}' not a valid method for minimum spanning tree creation."
            raise ValueError(msg)

        # instantiate the spanning tree as a network object
        net = Network(in_data=spanning_tree, weights_kws=weights_kws)

    return net


def mst_weighted_sort(net, maximum, net_kws):
    """Extract a minimum or maximum spanning tree from a network used
    the length-weighted sort method.

    Parameters
    ----------
    net : spaghetti.Network
        See ``spanning_tree()``.
    maximum : bool
        See ``spanning_tree()``.
    net_kws : dict
        Keywords arguments for instaniating a ``spaghetti.Network``.

    Returns
    -------
    spanning_tree : list
        All networks arcs that are members of the spanning tree.

    Notes
    -----

    This function is based on the method found in Chapter 3
    Section 4.3 of :cite:`Okabe2012`.

    """

    # network arcs dictionary sorted by arc length
    sort_kws = {"key": net.arc_lengths.get, "reverse": maximum}
    sorted_lengths = sorted(net.arc_lengths, **sort_kws)

    # the spanning tree is initially empty
    spanning_tree = []

    # iterate over each lengths of network arc
    while sorted_lengths:
        _arc = sorted_lengths.pop(0)
        # make a spatial representation of an arc
        chain_rep = util.chain_constr(net.vertex_coords, [_arc])
        # current set of network arcs as libpysal.cg.Chain
        _chains = spanning_tree + chain_rep
        # current network iteration
        _ntw = Network(in_data=_chains, **net_kws)
        # determine if the network contains a cycle
        if not util.network_has_cycle(_ntw.adjacencylist):
            # If no cycle is present, add the arc to the spanning tree
            spanning_tree.extend(chain_rep)

    return spanning_tree


@requires("geopandas", "shapely")
def element_as_gdf(
    net,
    vertices=False,
    arcs=False,
    pp_name=None,
    snapped=False,
    routes=None,
    id_col="id",
    geom_col=None,
):
    """Return a ``geopandas.GeoDataFrame`` of network elements. This can be
    (a) the vertices of a network; (b) the arcs of a network; (c) both the
    vertices and arcs of the network; (d) the raw point pattern associated
    with the network; (e) the snapped point pattern of (d); or (f) the
    shortest path routes between point observations.

    Parameters
    ----------
    net : spaghetti.Network
        A `spaghetti` network object.
    vertices : bool
        Extract the network vertices (``True``). Default is ``False``.
    arcs : bool
        Extract the network arcs (``True``). Default is ``False``.
    pp_name : str
        Name of the ``network.PointPattern`` to extract.
        Default is ``None``.
    snapped : bool
        If extracting a ``network.PointPattern``, set to ``True`` for
        snapped point locations along the network. Default is ``False``.
    routes : dict
        See ``paths`` from ``spaghetti.Network.shortest_paths``.
        Default is ``None``.
    id_col : str
        ``geopandas.GeoDataFrame`` column name for IDs. Default is ``"id"``.
        When extracting routes this creates an (origin, destination) tuple.
    geom_col : str
        Deprecated and will be removed in the minor release.
        ``geopandas.GeoDataFrame`` column name for IDs. Default is ``"id"``.
        When extracting routes this creates an (origin, destination) tuple.

    Raises
    ------
    KeyError
        In order to extract a ``network.PointPattern`` it must already
        be a part of the network object. This exception is raised
        when a ``network.PointPattern`` is being extracted that does
        not exist within the network object.

    Returns
    -------
    points : geopandas.GeoDataFrame
        Network point elements (either vertices or ``network.PointPattern``
        points) as a ``geopandas.GeoDataFrame`` of ``shapely.geometry.Point``
        objects with an ``"id"`` column and ``"geometry""`` column.
        If the network object has a ``network_component_vertices`` attribute,
        then component labels are also added in a column.
    lines : geopandas.GeoDataFrame
        Network arc elements as a ``geopandas.GeoDataFrame`` of
        ``shapely.geometry.LineString`` objects with an ``"id"``
        column and ``"geometry"`` column. If the network object has
        a ``network_component_labels`` attribute, then component labels
        are also added in a column.
    paths : geopandas.GeoDataFrame
        Shortest path routes along network arc elements as a
        ``geopandas.GeoDataFrame`` of ``shapely.geometry.LineString``
        objects with an ``"id"`` (see ``spaghetti.Network.shortest_paths()``)
        column and ``"geometry"`` column.

    Notes
    -----

    When both network vertices and arcs are desired, the variable
    declaration must be in the order: <vertices>, <arcs>.
    This function requires ``geopandas``.

    See also
    --------

    geopandas.GeoDataFrame

    Examples
    --------

    Instantiate a network object.

    >>> import spaghetti
    >>> from libpysal import examples
    >>> ntw = spaghetti.Network(examples.get_path("streets.shp"))

    Extract the network elements (vertices and arcs) as
    ``geopandas.GeoDataFrame`` objects.

    >>> vertices_df, arcs_df = spaghetti.element_as_gdf(
    ...     ntw, vertices=True, arcs=True
    ... )

    Examine the first vertex. It is a member of the component labeled ``0``.

    >>> vertices_df.loc[0]
    id                                            0
    geometry      POINT (728368.04762 877125.89535)
    comp_label                                    0
    Name: 0, dtype: object

    Calculate the total length of the network.

    >>> arcs_df.geometry.length.sum()
    104414.09200823458

    """

    # see GH#722
    if geom_col:
        dep_msg = (
            "The ``geom_col`` keyword argument is deprecated and will "
            "be dropped in the next minor release of pysal/spaghetti (1.8.0) "
            "in favor of the default 'geometry' name. Users can rename "
            "the geometry column following processing, if desired."
        )
        warnings.warn(dep_msg, FutureWarning, stacklevel=2)

    # shortest path routes between observations
    if routes:
        paths = util._routes_as_gdf(routes, id_col)
        # see GH#722
        if geom_col:
            paths.rename_geometry(geom_col, inplace=True)
        return paths

    # need vertices place holder to create network segment LineStrings
    # even if only network edges are desired.
    vertices_for_arcs = False
    if arcs and not vertices:
        vertices_for_arcs = True

    # vertices/nodes/points
    if vertices or vertices_for_arcs or pp_name:
        points = util._points_as_gdf(
            net,
            vertices,
            vertices_for_arcs,
            pp_name,
            snapped,
            id_col=id_col,
        )

        # return points geodataframe if arcs not specified or
        # if extracting `PointPattern` points
        if not arcs or pp_name:
            # see GH#722
            if geom_col:
                points.rename_geometry(geom_col, inplace=True)
            return points

    # arcs
    arcs = util._arcs_as_gdf(net, points, id_col=id_col)

    if vertices_for_arcs:
        # see GH#722
        if geom_col:
            arcs.rename_geometry(geom_col, inplace=True)
        return arcs

    else:
        # see GH#722
        if geom_col:
            points.rename_geometry(geom_col, inplace=True)
            arcs.rename_geometry(geom_col, inplace=True)
        return points, arcs


def regular_lattice(bounds, nh, nv=None, exterior=False):
    """Generate a regular lattice of line segments (``libpysal.cg.Chain``
    `objects <https://pysal.org/libpysal/generated/libpysal.cg.Chain.html>`_).

    Parameters
    ----------
    bounds : {tuple, list}
        Area bounds in the form - <minx,miny,maxx,maxy>.
    nh : int
        The number of internal horizontal lines of the lattice.
    nv : int
        The number of internal vertical lines of the lattice. Defaults to
        ``nh`` if left as None.
    exterior : bool
        Flag for including the outer bounding box segments. Default is False.

    Returns
    -------
    lattice : list
        The ``libpysal.cg.Chain`` objects forming a regular lattice.

    Notes
    -----

    The ``nh`` and ``nv`` parameters do not include the external
    line segments. For example, setting ``nh=3, nv=2, exterior=True``
    will result in 5 horizontal line sets and 4 vertical line sets.

    Examples
    --------

    Create a 5x5 regular lattice with an exterior

    >>> import spaghetti
    >>> lattice = spaghetti.regular_lattice((0,0,4,4), 3, exterior=True)
    >>> lattice[0].vertices
    [(0.0, 0.0), (1.0, 0.0)]

    Create a 5x5 regular lattice without an exterior

    >>> lattice = spaghetti.regular_lattice((0,0,5,5), 3, exterior=False)
    >>> lattice[-1].vertices
    [(3.75, 3.75), (3.75, 5.0)]

    Create a 7x9 regular lattice with an exterior from the
    bounds of ``streets.shp``.

    >>> path = libpysal.examples.get_path("streets.shp")
    >>> shp = libpysal.io.open(path)
    >>> lattice = spaghetti.regular_lattice(shp.bbox, 5, nv=7, exterior=True)
    >>> lattice[0].vertices
    [(723414.3683108028, 875929.0396895551), (724286.1381211297, 875929.0396895551)]

    """

    # check for bounds validity
    if len(bounds) != 4:
        err_msg = (
            f"The 'bounds' parameter is {len(bounds)} elements "
            "but should be exactly 4 - <minx,miny,maxx,maxy>."
        )
        raise RuntimeError(err_msg)

    # check for bounds validity
    if not nv:
        nv = nh
    try:
        nh, nv = int(nh), int(nv)
    except TypeError as err:
        err_msg = (
            f"The 'nh' and 'nv' parameters ({type(nh)}, {type(nv)}) "
            "could not be converted to integers."
        )
        raise TypeError(err_msg) from err

    # bounding box line lengths
    len_h, len_v = bounds[2] - bounds[0], bounds[3] - bounds[1]

    # horizontal and vertical increments
    incr_h, incr_v = len_h / float(nh + 1), len_v / float(nv + 1)

    # define the horizontal and vertical space
    space_h = [incr_h * slot for slot in range(nv + 2)]
    space_v = [incr_v * slot for slot in range(nh + 2)]

    # create vertical and horizontal lines
    lines_h = util.build_chains(space_h, space_v, exterior, bounds)
    lines_v = util.build_chains(space_h, space_v, exterior, bounds, h=False)

    # combine into one list
    lattice = lines_h + lines_v

    return lattice


class PointPattern:
    """A stub point pattern class used to store a point pattern.

    Note from the original author of ``pysal.network``:
    This class is monkey patched with network specific attributes when the
    points are snapped to a network. In the future this class may be
    replaced with a generic point pattern class.

    Parameters
    ----------
    in_data : {str, list, tuple, libpysal.cg.Point, geopandas.GeoDataFrame}
        The input geographic data. Either (1) a path to a shapefile
        (str); (2) an iterable containing ``libpysal.cg.Point``
        objects; (3) a single ``libpysal.cg.Point``; or
        (4) a ``geopandas.GeoDataFrame``.
    idvariable : str
        Field in the shapefile to use as an ID variable.
    attribute :  bool
        A flag to indicate whether all attributes are tagged to this
        class (``True``) or excluded (``False``). Default is ``False``.

    Attributes
    ----------
    points : dict
        Keys are the point IDs (int). Values are the :math:`(x,y)`
        coordinates (tuple).
    npoints : int
        The number of points.
    obs_to_arc : dict
        Keys are arc IDs (tuple). Values are snapped point information
        (``dict``).  Within the snapped point information (``dict``)
        keys are observation IDs (``int``), and values are snapped
        coordinates.
    obs_to_vertex : list
       List of incident network vertices to snapped observation points
       converted from a ``default_dict``. Originally in the form of
       paired left/right nearest network vertices {netvtx1: obs_id1,
       netvtx2: obs_id1, netvtx1: obs_id2... netvtx1: obs_idn}, then
       simplified to a list in the form
       [netvtx1, netvtx2, netvtx1, netvtx2, ...].
    dist_to_vertex : dict
        Keys are observations IDs (``int``). Values are distance lookup
        (``dict``). Within distance lookup (``dict``) keys are the two
        incident vertices of the arc and values are distance to each of
        those arcs.
    snapped_coordinates : dict
        Keys are the point IDs (int). Values are the snapped :math:`(x,y)`
        coordinates (tuple).
    snap_dist : bool
            Flag as ``True`` to include the distance from the original
            location to the snapped location along the network. Default
            is ``False``.

    """

    def __init__(self, in_data=None, idvariable=None, attribute=False):
        # initialize points dictionary and counter
        self.points = {}
        self.npoints = 0

        # determine input point data type
        in_dtype = str(type(in_data)).split("'")[1]
        # flag for points from a shapefile
        from_shp = False
        # flag for points as libpysal.cg.Point objects
        is_libpysal_points = False
        supported_iterables = ["list", "tuple"]
        # type error message
        msg = "'{}' not supported for point pattern instantiation."

        # set appropriate geometries
        if in_dtype == "str":
            from_shp = True
        elif in_dtype in supported_iterables:
            dtype = str(type(in_data[0])).split("'")[1]
            if dtype == "libpysal.cg.shapes.Point":
                is_libpysal_points = True
            else:
                raise TypeError(msg.format(dtype))
        elif in_dtype == "libpysal.cg.shapes.Point":
            in_data = [in_data]
            is_libpysal_points = True
        elif in_dtype == "geopandas.geodataframe.GeoDataFrame":
            from_shp = False
        else:
            raise TypeError(msg.format(str(in_dtype)))

        # either set native point ID from dataset or create new IDs
        if idvariable and not is_libpysal_points:
            ids = weights.util.get_ids(in_data, idvariable)
        else:
            ids = None

        # extract the point geometries
        if not is_libpysal_points:
            if from_shp:
                pts = _open(in_data)
            else:
                pts_objs = list(in_data.geometry)
                pts = [cg.shapes.Point((p.x, p.y)) for p in pts_objs]
        else:
            pts = in_data

        # fetch attributes if requested
        if attribute and not is_libpysal_points:
            # open the database file if data is from shapefile
            if from_shp:
                dbname = os.path.splitext(in_data)[0] + ".dbf"
                db = _open(dbname)

            # if data is from a GeoDataFrame, drop the geometry column
            # and declare attribute values as a list of lists
            else:
                db = in_data.drop(in_data.geometry.name, axis=1).values.tolist()
                db = [[d] for d in db]
        else:
            db = None

        # iterate over all points
        for i, pt in enumerate(pts):
            # IDs, attributes
            if ids and db is not None:
                self.points[ids[i]] = {"coordinates": pt, "properties": db[i]}

            # IDs, no attributes
            elif ids and db is None:
                self.points[ids[i]] = {"coordinates": pt, "properties": None}

            # no IDs, attributes
            elif not ids and db is not None:
                self.points[i] = {"coordinates": pt, "properties": db[i]}

            # no IDs, no attributes
            else:
                self.points[i] = {"coordinates": pt, "properties": None}

        # close the shapefile and database file
        # if the input data is a .shp
        if from_shp:
            pts.close()
            if db:
                db.close()

        # record number of points
        self.npoints = len(self.points.keys())


class SimulatedPointPattern:
    """Note from the original author of ``pysal.network``:
    Struct style class to mirror the ``PointPattern`` class.
    If the ``PointPattern`` class has methods, it might make
    sense to make this a child of that class. This class is not intended
    to be used by the external user.

    Attributes
    ----------
    npoints : int
        The number of points.
    obs_to_arc : dict
        Keys are arc IDs (tuple). Values are snapped point information
        (dict).  Within the snapped point information (dict)
        keys are observation IDs (int), and values are snapped
        coordinates.
    obs_to_vertex : list
       List of incident network vertices to snapped observation points
       converted from a default_dict. Originally in the form of
       paired left/right nearest network vertices {netvtx1: obs_id1,
       netvtx2: obs_id1, netvtx1: obs_id2... netvtx1: obs_idn}, then
       simplified to a list in the form
       [netvtx1, netvtx2, netvtx1, netvtx2, ...].
    dist_to_vertex : dict
        Keys are observations IDs (int). Values are distance lookup
        (dict). Within distance lookup (dict) keys are the two
        incident vertices of the arc and values are distance to each of
        those arcs.
    snapped_coordinates : dict
        Keys are the point IDs (int). Values are the snapped :math:`(x,y)`
        coordinates (tuple).
    snap_dist : bool
            Flag as ``True`` to include the distance from the original
            location to the snapped location along the network. Default
            is ``False``.

    """

    def __init__(self):
        # duplicate post-snapping PointPattern class structure
        self.npoints = 0
        self.obs_to_arc = {}
        self.obs_to_vertex = defaultdict(list)
        self.dist_to_vertex = {}
        self.snapped_coordinates = {}
