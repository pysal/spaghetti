from warnings import warn

import numpy
from libpysal import cg
from libpysal.common import requires
from rtree import Rtree

try:
    import geopandas
    import pandas
    import shapely
    from shapely.geometry import LineString
except ImportError:
    warn(
        "geopandas/shapely not available. Some functionality will be disabled.",
        UserWarning,
        stacklevel=1,
    )


def compute_length(v0, v1):
    """Compute the euclidean distance between two points.

    Parameters
    ----------
    v0 : tuple
        Coordinate sequence in the form x,y.
    vq : tuple
        Coordinate sequence in the form x,y.

    Returns
    -------
    euc_dist : float
        Euclidean distance.

    Examples
    --------

    >>> import spaghetti
    >>> point1, point2 = (0,0), (1,1)
    >>> spaghetti.util.compute_length(point1, point2)
    1.4142135623730951

    """

    euc_dist = cg.standalone.get_points_dist(v0, v1)

    return euc_dist


def get_neighbor_distances(ntw, v0, link):
    """Get distances to the nearest vertex neighbors along
    connecting arcs.

    Parameters
    ----------
    ntw : spaghetti.Network
        A spaghetti network object.
    v0 : int
        The vertex ID.
    link : dict
        The key is a tuple (start vertex, end vertex); value is ``float``.
        Cost per arc to travel, e.g. distance.

    Returns
    -------
    neighbors : dict
       The key is an integer (vertex ID); value is ``float`` (distance).

    Examples
    --------

    >>> import spaghetti
    >>> from libpysal import examples
    >>> ntw = spaghetti.Network(examples.get_path("streets.shp"))
    >>> neighs = spaghetti.util.get_neighbor_distances(ntw, 0, ntw.arc_lengths)
    >>> numpy.round(neighs[1], 10)
    102.6235345344

    """

    # fetch links associated with vertices
    arcs = ntw.enum_links_vertex(v0)

    # create neighbor distance lookup
    neighbors = {}

    # iterate over each associated link
    for arc in arcs:
        # set distance from vertex1 to vertex2 (link length)
        if arc[0] != v0:
            neighbors[arc[0]] = link[arc]
        else:
            neighbors[arc[1]] = link[arc]

    return neighbors


def generatetree(pred):
    """Rebuild the shortest path from root origin to destination.

    Parameters
    ----------
    pred : list
        List of preceding vertices for route traversal.

    Returns
    -------
    tree : dict
        The key is the root origin; value is the root origin to destination.

    Examples
    --------

    >>> import spaghetti
    >>> from libpysal import examples
    >>> ntw = spaghetti.Network(examples.get_path("streets.shp"))
    >>> distance, pred = spaghetti.util.dijkstra(ntw, 0)
    >>> tree = spaghetti.util.generatetree(pred)
    >>> tree[3]
    [23, 22, 20, 19, 170, 2, 0]

    """

    # instantiate tree lookup
    tree = {}

    # iterate over the list of predecessor vertices
    for i, p in enumerate(pred):
        # if the route begins/ends with itself set the
        # root vertex and continue to next iteration
        if p == -1:
            # tree keyed by root vertex with root vertex as path
            tree[i] = [i]
            continue

        # set the initial vertex `p` as `idx`
        idx = p
        # and add it as the first vertex in the path
        path = [idx]

        # iterate through the path until back to home vertex
        while idx >= 0:
            # set the next vertex on the path
            next_vertex = pred[idx]
            # and redeclare the current `idx`
            idx = next_vertex

            # add the vertex to path while not at home vertex
            if idx >= 0:
                path.append(next_vertex)

        # tree keyed by root vertex with network vertices as path
        tree[i] = path

    return tree


def dijkstra(ntw, v0, initial_dist=numpy.inf):
    """Compute the shortest path between a start vertex and
    all other vertices in an origin-destination matrix.

    Parameters
    ----------
    ntw :  spaghetti.Network
        A spaghetti network object.
    v0 : int
        Start vertex ID.
    initial_dist : float
        Integer break point to stop iteration and return n neighbors.
        Default is ``numpy.inf``.

    Returns
    -------
    distance : list
        List of distances from vertex to all other vertices.
    pred : list
        List of preceeding vertices for traversal route.

    Notes
    -----

    Based on :cite:`Dijkstra1959a`.

    Examples
    --------

    >>> import spaghetti
    >>> from libpysal import examples
    >>> ntw = spaghetti.Network(examples.get_path("streets.shp"))
    >>> distance, pred = spaghetti.util.dijkstra(ntw, 0)
    >>> round(distance[196], 4)
    5505.6682
    >>> pred[196]
    133

    """

    # cost per arc to travel, e.g. distance
    cost = ntw.arc_lengths

    # initialize travel costs as `inf` for all distances
    distance = [initial_dist for x in ntw.vertex_list]

    # label distance to self as 0
    distance[ntw.vertex_list.index(v0)] = 0

    # instantiate set of unvisited vertices
    unvisited = {v0}

    # initially label as predecessor vertices with -1 as path
    pred = [-1 for x in ntw.vertex_list]

    # iterate over `unvisited` until all vertices have been visited
    while len(unvisited) > 0:
        # get vertex with the lowest value from distance
        dist = initial_dist

        for vertex in unvisited:
            if distance[vertex] < dist:
                dist = distance[vertex]
                current = vertex

        # remove that vertex from the set
        unvisited.remove(current)

        # get the neighbors (and costs) to the current vertex
        neighbors = get_neighbor_distances(ntw, current, cost)

        # iterate over neighbors to find least cost along path
        for v1, indiv_cost in neighbors.items():
            # if the labeled cost is greater than
            # the currently calculated cost
            if distance[v1] > distance[current] + indiv_cost:
                # relabel to the currently calculated cost
                distance[v1] = distance[current] + indiv_cost

                # set the current vertex as a predecessor on the path
                pred[v1] = current

                # add the neighbor vertex to `unvisted`
                unvisited.add(v1)

    # cast preceding vertices list as an array of integers
    pred = numpy.array(pred, dtype=int)

    return distance, pred


def dijkstra_mp(ntw_vertex):
    """Compute the shortest path between a start vertex and all other
    vertices in the matrix utilizing multiple cores upon request.

    Parameters
    ----------
    ntw_vertex : tuple
        Tuple of arguments to pass into ``dijkstra()`` as
        (1) ``ntw`` - ``spaghetti.Network object``;
        (2) ``vertex`` - int (start node ID)

    Returns
    -------
    distance : list
        List of distances from vertex to all other vertices.
    pred : list
        List of preceeding vertices for traversal route.

    Notes
    -----

    Based on :cite:`Dijkstra1959a`.

    Examples
    --------

    >>> import spaghetti
    >>> from libpysal import examples
    >>> ntw = spaghetti.Network(examples.get_path("streets.shp"))
    >>> distance, pred = spaghetti.util.dijkstra_mp((ntw, 0))
    >>> round(distance[196], 4)
    5505.6682
    >>> pred[196]
    133

    """

    # unpack network object and source vertex
    ntw, vertex = ntw_vertex

    # calculate shortest path distances and predecessor vertices
    distance, pred = dijkstra(ntw, vertex)

    return distance, pred


def squared_distance_point_link(point, link):
    """Find the squared distance between a point and a link.

    Parameters
    ----------
    point : tuple
        Point coordinates (x,y).
    link : list
        List of 2 point coordinate tuples [(x0, y0), (x1, y1)].

    Returns
    -------
    sqd : float
        The distance squared between the point and edge.
    nearp : numpy.ndarray
        An array of (xb, yb); the nearest point on the edge.

    Examples
    --------

    >>> import spaghetti
    >>> point, link = (1,1), ((0,0), (2,0))
    >>> spaghetti.util.squared_distance_point_link(point, link)
    (1.0, array([1., 0.]))

    """

    # cast vertices comprising the network link as an array
    p0, p1 = (numpy.array(p) for p in link)

    # cast the observation point as an array
    p = numpy.array(point)

    # subtract point 0 coords from point 1
    v = p1 - p0
    # subtract point 0 coords from the observation coords
    w = p - p0

    # if the point 0 vertex is the closest point along the link
    c1 = numpy.dot(w, v)
    if c1 <= 0.0:
        sqd = numpy.dot(w.T, w)
        nearp = p0

        return sqd, nearp

    # if the point 1 vertex is the closest point along the link
    c2 = numpy.dot(v, v)
    if c2 <= c1:
        dp1 = p - p1
        sqd = numpy.dot(dp1.T, dp1)
        nearp = p1

        return sqd, nearp

    # otherwise the closest point along the link lies between p0 and p1
    b = c1 / c2
    bv = numpy.dot(b, v)
    pb = p0 + bv
    d2 = p - pb
    sqd = numpy.dot(d2, d2)
    nearp = pb

    return sqd, nearp


def snap_points_to_links(points, links):
    """Place points onto closest link in a set of links (arc/edges).

    Parameters
    ----------
    points : dict
        Point ID as key and (x,y) coordinate as value.
    links : list
        Elements are of type ``libpysal.cg.shapes.Chain``
        ** Note ** each element is a link represented as a chain with
        *one head and one tail vertex* in other words one link only.

    Returns
    -------
    point2link : dict
        Key [point ID (see points in arguments)]; value [a 2-tuple
        ((head, tail), point) where (head, tail) is the target link,
        and point is the snapped location on the link.

    Examples
    --------

    >>> import spaghetti
    >>> from libpysal.cg.shapes import Point, Chain
    >>> points = {0: Point((1,1))}
    >>> link = [Chain([Point((0,0)), Point((2,0))])]
    >>> spaghetti.util.snap_points_to_links(points, link)
    {0: ([(0.0, 0.0), (2.0, 0.0)], array([1., 0.]))}

    """

    # instantiate an rtree
    rtree = Rtree()
    # set the smallest possible float epsilon on machine
    SMALL = numpy.finfo(float).eps

    # initialize network vertex to link lookup
    vertex_2_link = {}

    # iterate over network links
    for i, link in enumerate(links):
        # extract network link (x,y) vertex coordinates
        head, tail = link.vertices
        x0, y0 = head
        x1, y1 = tail

        if (x0, y0) not in vertex_2_link:
            vertex_2_link[(x0, y0)] = []

        if (x1, y1) not in vertex_2_link:
            vertex_2_link[(x1, y1)] = []

        vertex_2_link[(x0, y0)].append(link)
        vertex_2_link[(x1, y1)].append(link)

        # minimally increase the bounding box exterior
        bx0, by0, bx1, by1 = link.bounding_box
        bx0 -= SMALL
        by0 -= SMALL
        bx1 += SMALL
        by1 += SMALL

        # insert the network link and its associated
        # rectangle into the rtree
        rtree.insert(i, (bx0, by0, bx1, by1), obj=link)

    # build a KDtree on link vertices
    kdtree = cg.KDTree(list(vertex_2_link.keys()))

    point2link = {}

    for pt_idx, point in points.items():
        # first, find nearest neighbor link vertices for the point
        dmin, vertex = kdtree.query(point, k=1)
        vertex = tuple(kdtree.data[vertex])
        closest = vertex_2_link[vertex][0].vertices

        # Use this link as the candidate closest link:  closest
        # Use the distance as the distance to beat:     dmin
        point2link[pt_idx] = (closest, numpy.array(vertex))
        x0 = point[0] - dmin
        y0 = point[1] - dmin
        x1 = point[0] + dmin
        y1 = point[1] + dmin

        # Find all links with bounding boxes that intersect
        # a query rectangle centered on the point with sides
        # of length dmin * dmin
        rtree_lookup = rtree.intersection([x0, y0, x1, y1], objects=True)
        candidates = [cand.object.vertices for cand in rtree_lookup]

        # Sorting the candidate ensures reproducible results from OS to OS.
        # See:
        #   https://github.com/pysal/spaghetti/pull/595
        #   https://github.com/pysal/spaghetti/issues/598
        #   https://github.com/pysal/spaghetti/pull/599
        candidates.sort(reverse=True)
        dmin += SMALL
        dmin2 = dmin * dmin

        # of the candidate arcs, find the nearest to the query point
        for candidate in candidates:
            dist2cand, nearp = squared_distance_point_link(point, candidate)
            if dist2cand <= dmin2:
                closest = candidate
                dmin2 = dist2cand
                point2link[pt_idx] = (closest, nearp)

    return point2link


def network_has_cycle(adjacency):
    """Searches for a cycle in the complete network/graph.

    Parameters
    ----------
    adjacency : spaghetti.Network.adjacencylist
        Vertex adjacency relationships.

    Returns
    -------
    network_cycle_found : bool
        ``True`` for a cycle being found in the network/graph,
        otherwise ``False``.

    """

    def tree_has_cycle(_parent, _v):
        """Searches for a cycle in the subtree.

        Parameters
        ----------
        _parent : int
            Root vertex for the subnetwork/graph.
        _v : int
            Current vertex index of in the complete network.

        Returns
        -------
        subtree_cycle_found : bool
            Current recursion found a cycle in the subtree.

        """

        # Set current cycle tag as False
        subtree_cycle_found = False

        # Label the current network vertex as seen
        seen[_v] = True

        # Perform recursion for all adjacent network/graph vertices
        for rv in adjacency[_v]:
            # If vertex already seen, skip it
            if not seen[rv]:
                # Perform recursion down the depth-first search tree
                if tree_has_cycle(_v, rv):
                    subtree_cycle_found = True
                    break

            # If an adjacent vertex has not been seen and it is not the
            # parent of current vertex, then a cycle is present
            elif _parent != rv:
                subtree_cycle_found = True
                break

        return subtree_cycle_found

    # set initial cycle tag as False
    network_cycle_found = False

    # Label all network/graph vertices as not seen
    vids = list(adjacency.keys())
    seen = {vid: False for vid in vids}

    # Perform depth-first search recursion to isolate cycles
    for v in vids:
        # If vertex already seen, skip it; or recurse down the depth-first search tree
        if not seen[v] and tree_has_cycle(-1, v):
            network_cycle_found = True
            break

    return network_cycle_found


def chain_constr(vcoords, arcs):
    """Create the spatial representation of a network arc.

    Parameters
    ----------
    vcoords : dict
        Vertex to coordinate lookup (see ``spaghetti.Network.vertex_coords``).
    arcs : list
        Arcs represented as start and end vertices.

    Returns
    -------
    spatial_reps : list
        Spatial representations of arcs - ``libpysal.cg.Chain`` objects.

    """
    spatial_reps = [_chain_constr(vcoords, vs) for vs in arcs]
    return spatial_reps


def _chain_constr(_vcoords, _vs):
    """Construct a libpysal.cg.Chain object.

    Parameters
    ----------
    _vcoords : {dict, None}
        See ``vcoords`` in ``get_chains()``.
    _vs : tuple
        Start and end vertex IDs of arc.

    Returns
    -------
    libpysal.cg.Chain
        Spatial representation of the arc.

    """

    return cg.Chain([cg.Point(_vcoords[v]) for v in _vs] if _vcoords else _vs)


def build_chains(space_h, space_v, exterior, bounds, h=True):
    """Generate line segments for a lattice.

    Parameters
    ----------
    space_h : list
        Horizontal spacing.
    space_v : list
        Vertical spacing.
    exterior : bool
        Flag for including the outer bounding box segments.
    bounds : list
        Area bounds in the form - <minx,miny,maxx,maxy>.
    h : bool
        Generate horizontal line segments.
        Default is ``True``. ``False`` generates vertical segments.

    Returns
    -------
    chains : list
        All horizontal or vertical line segments in the lattice.

    """

    # Initialize starting and ending indices
    start_h, end_h, start_v, end_v = 0, len(space_h), 0, len(space_v)

    # set inital index track back to 0
    minus_y, minus_x = 0, 0

    if h:  # start track back at 1 for horizontal lines
        minus_x = 1
        if not exterior:  # do not include borders
            start_v += 1
            end_v -= 1

    else:  # start track back at 1 for vertical lines
        minus_y = 1
        if not exterior:  # do not include borders
            start_h += 1
            end_h -= 1

    # Create empty line list and fill
    chains = []

    # for element in the horizontal index
    for plus_h in range(start_h, end_h):
        # for element in the vertical index
        for plus_v in range(start_v, end_v):
            # ignore if a -1 index
            if plus_h - minus_x == -1 or plus_v - minus_y == -1:
                continue
            else:
                # Point 1 (start point + previous slot in
                #          horizontal or vertical space index)
                p1x = bounds[0] + space_h[plus_h - minus_x]
                p1y = bounds[1] + space_v[plus_v - minus_y]
                p1 = cg.Point((p1x, p1y))

                # Point 2 (start point + current slot in
                #          horizontal or vertical space index)
                p2x = bounds[0] + space_h[plus_h]
                p2y = bounds[1] + space_v[plus_v]
                p2 = cg.Point((p2x, p2y))

                # libpysal.cg.Chain
                chains.append(_chain_constr(None, [p1, p2]))

    return chains


@requires("geopandas", "shapely")
def _points_as_gdf(net, vertices, vertices_for_arcs, pp_name, snapped, id_col=None):
    """Internal function for returning a point ``geopandas.GeoDataFrame``
    called from within ``spaghetti.element_as_gdf()``.

    Parameters
    ----------
    vertices_for_arcs : bool
        Flag for points being an object returned (``False``) or for merely
        creating network arcs (``True``). Set from within the parent
        function (``spaghetti.element_as_gdf()``).

    Raises
    ------

    KeyError
        In order to extract a ``network.PointPattern`` it must already
        be a part of the network object. This exception is raised
        when a ``network.PointPattern`` is being extracted that does not
        exist within the network object.

    Returns
    -------
    points : geopandas.GeoDataFrame
        Network point elements (either vertices or ``network.PointPattern``
        points) as a simple ``geopandas.GeoDataFrame`` of
        ``shapely.geometry.Point`` objects with an ``"id"`` column and
        ``"geometry"`` column.

    Notes
    -----

    1. See ``spaghetti.element_as_gdf()`` for description of arguments.
    2. This function requires ``geopandas``.

    """

    # vertices / nodes
    if vertices or vertices_for_arcs:
        pts_dict = net.vertex_coords

    if pp_name:
        try:
            pp = net.pointpatterns[pp_name]
        except KeyError as err:
            err_msg = f"Available point patterns are {net.pointpatterns.keys()}"
            raise KeyError(err_msg) from err

        # raw point pattern
        if not snapped:
            pp_pts = pp.points
            n_pp_pts = range(len(pp_pts))
            pts_dict = {point: pp_pts[point]["coordinates"] for point in n_pp_pts}

        # snapped point pattern
        else:
            pts_dict = pp.snapped_coordinates

    # instantiate geopandas.GeoDataFrame
    points = geopandas.GeoDataFrame(
        pts_dict.keys(),
        columns=[id_col],
        geometry=shapely.points(numpy.asarray(list(pts_dict.values()))),
    )

    # additional columns
    if not pp_name:
        ncv_tag = "network_component_vertices"
        if hasattr(net, ncv_tag):
            ncv = getattr(net, ncv_tag)
            ncv_map = {v: k for k, verts in ncv.items() for v in verts}
            points["comp_label"] = points[id_col].map(ncv_map)
    if pp_name:
        c2o_tag = "component_to_obs"
        if hasattr(pp, c2o_tag):
            c2o = getattr(pp, c2o_tag)
            o2c_map = {o: c for c, obs in c2o.items() for o in obs}
            points["comp_label"] = points[id_col].map(o2c_map)

    return points


@requires("geopandas", "shapely")
def _arcs_as_gdf(net, points, id_col=None):
    """Internal function for returning an arc ``geopandas.GeoDataFrame``
    called from within ``spaghetti.element_as_gdf()``.

    Returns
    -------
    arcs : geopandas.GeoDataFrame
        Network arc elements as a ``geopandas.GeoDataFrame`` of
        ``shapely.geometry.LineString`` objects with an ``"id"``
        column and ``geometry`` column.

    Notes
    -----

    1. See ``spaghetti.element_as_gdf()`` for description of arguments.
    2. This function requires ``geopandas``.

    """

    def _line_coords(loc):
        return (
            (points.loc[loc[0]].geometry.x, points.loc[loc[0]].geometry.y),
            (points.loc[loc[1]].geometry.x, points.loc[loc[1]].geometry.y),
        )

    # instantiate GeoDataFrame
    arcs = pandas.DataFrame(zip(sorted(net.arcs)), columns=[id_col])
    arcs = arcs.set_geometry(
        shapely.linestrings(arcs[id_col].map(_line_coords).values.tolist())
    )

    # additional columns
    if hasattr(net, "network_component_labels"):
        arcs["comp_label"] = net.network_component_labels

    return arcs


@requires("geopandas", "shapely")
def _routes_as_gdf(paths, id_col):
    """Internal function for returning a shortest paths
    ``geopandas.GeoDataFrame`` called from within
    ``spaghetti.element_as_gdf()``.

    Returns
    -------
    paths : geopandas.GeoDataFrame
        Network shortest paths as a ``geopandas.GeoDataFrame`` of
        ``shapely.geometry.LineString`` objects with an ``"O"`` (origin),
        ``D`` (destination), and ``geometry`` column. An additional
        column storing the ID as a tuple is available.

    Notes
    -----

    1. See ``spaghetti.element_as_gdf()`` for description of arguments.
    2. This function requires ``geopandas``.

    """

    # instantiate as a geodataframe
    paths = dict(paths)
    ids, geoms = zip(paths.keys()), [LineString(g.vertices) for g in paths.values()]
    paths = geopandas.GeoDataFrame(ids, columns=[id_col], geometry=geoms)

    return paths
