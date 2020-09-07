from libpysal import cg, examples, io
import numpy
import unittest
import copy

try:
    import geopandas

    GEOPANDAS_EXTINCT = False
except ImportError:
    GEOPANDAS_EXTINCT = True

# empirical data ---------------------------------------------------------------
# network shapefile
STREETS = examples.get_path("streets.shp")
# observations
schools = "schools"
SCHOOLS = examples.get_path(schools + ".shp")
crimes = "crimes"
CRIMES = examples.get_path(crimes + ".shp")

# native pysal geometries ------------------------------------------------------
P00 = cg.Point((0, 0))
P03 = cg.Point((0, 3))
P030001 = cg.Point((0, 3.0001))
P01 = cg.Point((0, 1))
P10 = cg.Point((1, 0))
P11 = cg.Point((1, 1))
P12 = cg.Point((1, 2))
P21 = cg.Point([2, 1])
P22 = cg.Point((2, 2))
P33 = cg.Point((3, 3))
P34 = cg.Point((3, 4))
P40 = cg.Point((4, 0))
P400010 = cg.Point((4.0001, 0))
P11P22_CHAIN = cg.Chain([P11, P22])

P0505 = cg.Point([0.5, 0.5])
P052 = cg.Point([0.5, 2.0])
P0525 = cg.Point([0.5, 2.5])
P2505 = cg.Point([2.5, 0.5])
P2525 = cg.Point([2.5, 2.5])

P31 = cg.Point([3, 1])
P325125 = cg.Point([3.25, 1.25])
P3375125 = cg.Point([3.375, 1.25])
P35125 = cg.Point([3.5, 1.25])
P3751 = cg.Point([3.75, 1])
P35075 = cg.Point([3.5, 0.75])
P3375075 = cg.Point([3.375, 0.75])
P325075 = cg.Point([3.25, 0.75])

RING = [
    cg.Chain([P31, P325125]),
    cg.Chain([P325125, P3375125]),
    cg.Chain([P3375125, P35125]),
    cg.Chain([P35125, P3751]),
    cg.Chain([P3751, P35075]),
    cg.Chain([P35075, P3375075]),
    cg.Chain([P3375075, P325075]),
    cg.Chain([P325075, P31]),
]
EXTENSION = [cg.Chain([P12, P22, P21])]

GOOD_TRIANGLE = [
    cg.Chain([P00, P03]),
    cg.Chain([P03, P40]),
    cg.Chain([P40, P00]),
]

BAD_TRIANGLE = [
    cg.Chain([P00, P03]),
    cg.Chain([P030001, P400010]),
    cg.Chain([P40, P00]),
]

synth_obs = "synth_obs"
points1 = "points1"
points2 = "points2"


# -------------------------------------------------------------------------------
class TestNetwork(unittest.TestCase):
    def setUp(self):
        # empirical network instantiated from shapefile
        self.ntw_shp = self.spaghetti.Network(in_data=STREETS, weightings=True)
        self.n_known_shp_arcs, self.n_known_shp_vertices = 303, 230

        # native pysal geometries
        self.chains_from_shp = [
            cg.Chain([cg.Point(self.ntw_shp.vertex_coords[vertex]) for vertex in arc])
            for arc in self.ntw_shp.arcs
        ]

        # lattice and ring + extension
        bounds, h, v = (0, 0, 2, 2), 1, 1
        self.lattice = self.spaghetti.regular_lattice(bounds, h, nv=v)
        self.lines = self.lattice + RING + EXTENSION
        self.ntw_from_lattice_ring = self.spaghetti.Network(in_data=self.lines)
        self.ntw_from_lattice_ring.snapobservations([P0505, P052], "points")

        # Pythagorean Triple
        self.triangle = self.spaghetti.Network(in_data=GOOD_TRIANGLE)

    def tearDown(self):
        pass

    def test_network_data_read(self):
        # shp test against known
        self.assertEqual(len(self.ntw_shp.arcs), self.n_known_shp_arcs)
        self.assertEqual(len(self.ntw_shp.vertices), self.n_known_shp_vertices)

        arc_lengths = self.ntw_shp.arc_lengths.values()
        self.assertAlmostEqual(sum(arc_lengths), 104414.0920159, places=5)

        self.assertIn(0, self.ntw_shp.adjacencylist[1])
        self.assertIn(0, self.ntw_shp.adjacencylist[2])
        self.assertNotIn(0, self.ntw_shp.adjacencylist[3])

    def test_network_from_libpysal_chains(self):
        known_components = self.ntw_shp.network_n_components
        known_length = sum(self.ntw_shp.arc_lengths.values())
        # network instantiated from libpysal.cg.Chain objects
        for dtype in (list, tuple, numpy.array):
            ntw_data = dtype(self.chains_from_shp)
            self.ntw_from_chains = self.spaghetti.Network(in_data=ntw_data)
            observed_components = self.ntw_from_chains.network_n_components
            observed_length = sum(self.ntw_from_chains.arc_lengths.values())
            self.assertEqual(observed_components, known_components)
            self.assertAlmostEqual(observed_length, known_length, places=3)

    def test_network_from_single_libpysal_chain(self):
        # network instantiated from a single libpysal.cg.Chain
        self.ntw_chain_out = self.spaghetti.Network(in_data=P11P22_CHAIN)
        known_edges = self.ntw_chain_out.edges
        fname = "test_network.pkl"
        # test save and load network
        self.ntw_chain_out.savenetwork(fname)
        self.ntw_chain_in = self.spaghetti.Network.loadnetwork(fname)
        observed_arcs = self.ntw_chain_in.arcs
        self.assertEqual(observed_arcs, known_edges)

    def test_network_from_vertical_libpysal_chains(self):
        vert_up = cg.Chain([P0505, P052])
        self.ntw_up_chain = self.spaghetti.Network(in_data=vert_up)
        self.assertEqual(len(self.ntw_up_chain.arcs), len(vert_up.segments))

        vert_down = cg.Chain([P052, P0505])
        self.ntw_down_chain = self.spaghetti.Network(in_data=vert_down)
        self.assertEqual(len(self.ntw_down_chain.arcs), len(vert_down.segments))

    def test_network_failures(self):
        # try instantiating network with single point
        with self.assertRaises(TypeError):
            self.spaghetti.Network(in_data=P11)
        # try instantiating network with list of single point
        with self.assertRaises(TypeError):
            self.spaghetti.Network(in_data=[P11])

    @unittest.skipIf(GEOPANDAS_EXTINCT, "Missing Geopandas")
    def test_network_from_geopandas(self):
        # network instantiated from geodataframe
        gdf = geopandas.read_file(STREETS)
        self.ntw_gdf = self.spaghetti.Network(in_data=gdf, w_components=True)

        # gdf test against known
        self.assertEqual(len(self.ntw_gdf.arcs), self.n_known_shp_arcs)
        self.assertEqual(len(self.ntw_gdf.vertices), self.n_known_shp_vertices)

        # shp against gdf
        self.assertEqual(len(self.ntw_shp.arcs), len(self.ntw_gdf.arcs))
        self.assertEqual(len(self.ntw_shp.vertices), len(self.ntw_gdf.vertices))

    def test_round_sig(self):
        # round to 2 significant digits test
        x_round2, y_round2 = 1200, 1900
        self.ntw_shp.vertex_sig = 2
        obs_xy_round2 = self.ntw_shp._round_sig((1215, 1865))
        self.assertEqual(obs_xy_round2, (x_round2, y_round2))

        # round to no significant digits test
        x_roundNone, y_roundNone = 1215, 1865
        self.ntw_shp.vertex_sig = None
        obs_xy_roundNone = self.ntw_shp._round_sig((1215, 1865))
        self.assertEqual(obs_xy_roundNone, (x_roundNone, y_roundNone))

    def test_vertex_atol(self):
        known_components = 1
        ntw_triangle = self.spaghetti.Network(in_data=BAD_TRIANGLE, vertex_atol=2)
        observed_components = ntw_triangle.network_n_components
        self.assertEqual(observed_components, known_components)

    def test_contiguity_weights(self):
        known_network_histo = [(2, 35), (3, 89), (4, 105), (5, 61), (6, 13)]
        observed_network_histo = self.ntw_shp.w_network.histogram
        self.assertEqual(known_network_histo, observed_network_histo)

        known_graph_histo = [(2, 2), (3, 2), (4, 47), (5, 80), (6, 48)]
        observed_graph_histo = self.ntw_shp.w_graph.histogram
        self.assertEqual(observed_graph_histo, known_graph_histo)

    def test_components(self):
        known_network_arc = (225, 226)
        observed_network_arc = self.ntw_shp.network_component2arc[0][-1]
        self.assertEqual(observed_network_arc, known_network_arc)

        known_graph_edge = (207, 208)
        observed_graph_edge = self.ntw_shp.graph_component2edge[0][-1]
        self.assertEqual(observed_graph_edge, known_graph_edge)

    def test_connected_components(self):
        ## test warnings
        ntw = copy.deepcopy(self.ntw_from_lattice_ring)

        # observed values
        observed_connected = ntw.network_fully_connected
        # known values
        known_connected = False
        self.assertEqual(observed_connected, known_connected)

        # observed values
        observed_component_vertices = ntw.network_component_vertices
        # known values
        known_component_vertices = {
            0: [0, 1, 2, 3, 4, 13],
            1: [5, 6, 7, 8, 9, 10, 11, 12],
        }
        self.assertEqual(observed_component_vertices, known_component_vertices)

        # observed values
        observed_network_vtx = ntw.network_component_vertex_count
        observed_graph_vtx = ntw.graph_component_vertex_count
        # known values
        known_network_vtx = {0: 6, 1: 8}
        known_graph_vtx = {0: 3, 1: 8}
        self.assertEqual(observed_network_vtx, known_network_vtx)
        self.assertEqual(observed_graph_vtx, known_graph_vtx)

        # observed values
        observed_edge_lengths = ntw.edge_lengths[(0, 1)]
        # known values
        known_edge_lengths = 1.0
        self.assertEqual(observed_edge_lengths, known_edge_lengths)

        # observed values
        observed_largest_net = ntw.network_largest_component
        observed_longest_graph = ntw.graph_longest_component
        # known values
        known_largest = 1
        known_longest = 0
        self.assertEqual(observed_largest_net, known_largest)
        self.assertEqual(observed_longest_graph, known_longest)

        # observed values
        observed_lengths = ntw.network_component_lengths
        # known values
        known_lengths = {0: 6.0, 1: 1.914213562373095}
        self.assertEqual(observed_lengths, known_lengths)

    def test_distance_band_weights(self):
        w = self.ntw_shp.distancebandweights(threshold=500)
        self.assertEqual(w.n, 230)
        self.assertEqual(
            w.histogram,
            [(1, 22), (2, 58), (3, 63), (4, 40), (5, 36), (6, 3), (7, 5), (8, 3)],
        )

    def test_split_arcs_200(self):
        n200 = self.ntw_shp.split_arcs(200.0)
        self.assertEqual(len(n200.arcs), 688)

    def test_enum_links_vertex(self):
        coincident = self.ntw_shp.enum_links_vertex(24)
        self.assertIn((24, 48), coincident)

    def test_shortest_paths(self):
        # symmetric point pattern
        known_vertices = 10
        self.ntw_shp.snapobservations(SCHOOLS, schools)
        _, tree = self.ntw_shp.allneighbordistances(schools, gen_tree=True)
        observed_paths = self.ntw_shp.shortest_paths(tree, schools)
        observed_vertices = len(observed_paths[0][1].vertices)
        self.assertEqual(observed_vertices, known_vertices)

        # asymmetric point pattern
        bounds, h, v = (0, 0, 3, 3), 2, 2
        lattice = self.spaghetti.regular_lattice(bounds, h, nv=v, exterior=False)
        ntw = self.spaghetti.Network(in_data=lattice)

        POINTS1 = [P0505, P2525]
        POINTS2 = [P0525, P2505, cg.Point((0.75, 0.6))]
        ntw.snapobservations(POINTS1, points1)
        ntw.snapobservations(POINTS2, points2)
        _, tree = ntw.allneighbordistances(points1, points2, gen_tree=True)
        observed_paths = ntw.shortest_paths(tree, points1, pp_dest=points2)

        # observed values
        observed_vertices1 = observed_paths[2][1].vertices
        observed_vertices2 = len(observed_paths[3][1].vertices)
        # known values
        known_vertices1 = [(1.0, 0.5), (1.0, 0.6)]
        known_vertices2 = 4
        self.assertEqual(observed_vertices1, observed_vertices1)
        self.assertEqual(observed_vertices2, known_vertices2)

        # test error
        with self.assertRaises(AttributeError):
            lattice = self.spaghetti.regular_lattice((0, 0, 4, 4), 4)
            ntw = self.spaghetti.Network(in_data=lattice)
            paths = ntw.shortest_paths([], synth_obs)

    def test_extract_component(self):
        ntw = copy.deepcopy(self.ntw_from_lattice_ring)
        s2s, tree = ntw.allneighbordistances("points", gen_tree=True)

        # test longest component
        longest = self.spaghetti.extract_component(ntw, ntw.network_longest_component)
        # observed values
        observed_napts = longest.non_articulation_points
        # known values
        known_napts = [2, 4, 13]
        self.assertEqual(observed_napts, known_napts)

        # test largest component
        largest = self.spaghetti.extract_component(ntw, ntw.network_largest_component)
        # observed values
        observed_arcs, observed_edges = largest.arcs, largest.edges
        # known values
        known_arcs = [
            (5, 6),
            (5, 12),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10),
            (10, 11),
            (11, 12),
        ]
        known_edges = known_arcs
        self.assertEqual(observed_arcs, known_arcs)
        self.assertEqual(observed_arcs, known_edges)
        self.assertEqual(observed_edges, known_arcs)
        self.assertEqual(observed_edges, known_edges)

    def test_spanning_tree(self):
        # minimum
        known_len = 7.0
        mst = self.spaghetti.spanning_tree(
            self.triangle, method="sort", maximum=False, silence_warnings=True
        )
        observed_len = sum(mst.arc_lengths.values())
        self.assertEqual(observed_len, known_len)

        # maximum
        known_len = 9.0
        mst = self.spaghetti.spanning_tree(
            self.triangle, method="sort", maximum=True, silence_warnings=True
        )
        observed_len = sum(mst.arc_lengths.values())
        self.assertEqual(observed_len, known_len)

        # method error
        with self.assertRaises(ValueError):
            self.spaghetti.spanning_tree(self.triangle, method="tors")

    @unittest.skipIf(GEOPANDAS_EXTINCT, "Missing Geopandas")
    def test_element_as_gdf(self):
        # extract both vertices and arcs
        vertices, arcs = self.spaghetti.element_as_gdf(
            self.ntw_shp, vertices=True, arcs=True
        )
        # test arcs
        known_vertex_wkt = "POINT (728368.04762 877125.89535)"
        observed_vertex = vertices.loc[(vertices["id"] == 0), "geometry"].squeeze()
        observed_vertex_wkt = observed_vertex.wkt
        self.assertEqual(observed_vertex_wkt, known_vertex_wkt)
        # test arcs
        known_arc_wkt = (
            "LINESTRING (728368.04762 877125.89535, 728368.13931 877023.27186)"
        )
        observed_arc = arcs.loc[(arcs["id"] == (0, 1)), "geometry"].squeeze()
        observed_arc_wkt = observed_arc.wkt
        self.assertEqual(observed_arc_wkt, known_arc_wkt)

        # extract only arcs
        arcs = self.spaghetti.element_as_gdf(self.ntw_shp, arcs=True)
        observed_arc = arcs.loc[(arcs["id"] == (0, 1)), "geometry"].squeeze()
        observed_arc_wkt = observed_arc.wkt
        self.assertEqual(observed_arc_wkt, known_arc_wkt)

        # extract symmetric routes
        known_length, bounds, h, v = 2.6, (0, 0, 3, 3), 2, 2
        lattice = self.spaghetti.regular_lattice(bounds, h, nv=v, exterior=False)
        ntw = self.spaghetti.Network(in_data=lattice)
        SYNTH_OBS = [cg.Point([0.2, 1.3]), cg.Point([0.2, 1.7]), cg.Point([2.8, 1.5])]
        ntw.snapobservations(SYNTH_OBS, synth_obs)
        _, tree = ntw.allneighbordistances(synth_obs, gen_tree=True)
        paths = ntw.shortest_paths(tree, synth_obs)
        paths_gdf = self.spaghetti.element_as_gdf(ntw, routes=paths)
        observed_length = paths_gdf.loc[0, "geometry"].length
        self.assertEqual(observed_length, known_length)

        # extract asymmetric routes
        known_origins, bounds, h, v = 2, (0, 0, 3, 3), 2, 2
        lattice = self.spaghetti.regular_lattice(bounds, h, nv=v, exterior=False)
        ntw = self.spaghetti.Network(in_data=lattice)
        POINTS1 = [P0505, P2525]
        POINTS2 = [P0525, P2505]
        ntw.snapobservations(POINTS1, points1)
        ntw.snapobservations(POINTS2, points2)
        _, tree = ntw.allneighbordistances(points1, points2, gen_tree=True)
        paths = ntw.shortest_paths(tree, points1, pp_dest=points2)
        paths_gdf = self.spaghetti.element_as_gdf(ntw, routes=paths)
        observed_origins = paths_gdf["O"].nunique()
        self.assertEqual(observed_origins, known_origins)

    def test_regular_lattice(self):
        # 4x4 regular lattice with the exterior
        known = [P00, P10]
        bounds = (0, 0, 3, 3)
        lattice = self.spaghetti.regular_lattice(bounds, 2, nv=2, exterior=True)
        observed = lattice[0].vertices
        self.assertEqual(observed, known)

        # 5x5 regular lattice without the exterior
        known = [P33, P34]
        bounds = (0, 0, 4, 4)
        lattice = self.spaghetti.regular_lattice(bounds, 3, exterior=False)
        observed = lattice[-1].vertices
        self.assertEqual(observed, known)

        # 7x9 regular lattice from shapefile bounds
        known_vertices = [
            (723414.3683108028, 875929.0396895551),
            (724286.1381211297, 875929.0396895551),
        ]
        shp = io.open(STREETS)
        lattice = self.spaghetti.regular_lattice(shp.bbox, 5, nv=7, exterior=True)
        observed_vertices = lattice[0].vertices
        for observed, known in zip(observed_vertices, known_vertices):
            self.assertEqual((observed[0], observed[1]), known)

        # test for Type Error
        with self.assertRaises(TypeError):
            self.spaghetti.regular_lattice(bounds, [[4]])

        # test for Runtime Error
        with self.assertRaises(RuntimeError):
            self.spaghetti.regular_lattice((0, 0, 1), 1)


# -------------------------------------------------------------------------------
class TestNetworkPointPattern(unittest.TestCase):
    def setUp(self):
        self.ntw = self.spaghetti.Network(in_data=STREETS)
        self.obs = [schools, crimes]
        self.OBS = [SCHOOLS, CRIMES]
        self.idxs = ["pp1", "pp2"]
        iterator = zip(self.obs, self.OBS, self.idxs)
        for (obs, OBS, idx) in iterator:
            self.ntw.snapobservations(OBS, obs, attribute=True)
            setattr(self, idx, self.ntw.pointpatterns[obs])

        self.known_pp1_npoints = 8

    def tearDown(self):
        pass

    def test_pp_from_libpysal_points(self):
        # known
        cpp = self.ntw.pointpatterns[crimes]
        known_snapped = set(cpp.snapped_coordinates.values())
        cg_crimes = "cg_%s" % crimes
        # points from pysal geometries
        points = [cg.Point(cpp.points[i]["coordinates"]) for i in cpp.points]
        for dtype in (list, tuple):
            point_data = dtype(points)
            self.ntw.snapobservations(point_data, cg_crimes, attribute=True)
            observed = self.ntw.pointpatterns[cg_crimes]
            observed_snapped = set(observed.snapped_coordinates.values())
            self.assertEqual(observed_snapped, known_snapped)

    def test_pp_from_single_libpysal_point(self):
        # network instantiated from a single libpysal.cg.Chain
        known_dist = 1.4142135623730951
        self.ntw_from_chain = self.spaghetti.Network(in_data=P11P22_CHAIN)
        self.ntw_from_chain.snapobservations(P00, synth_obs)
        snap_dist = self.ntw_from_chain.pointpatterns[synth_obs].dist_snapped[0]
        self.assertAlmostEqual(snap_dist, known_dist, places=10)

        # network instantiated from a single vertical (up) libpysal.cg.Chain
        chain = cg.Chain([P11, P12])
        known_dist = 1.0
        self.ntw_from_chain = self.spaghetti.Network(in_data=chain)
        self.ntw_from_chain.snapobservations(cg.Point((0, 1.5)), synth_obs)
        snap_dist = self.ntw_from_chain.pointpatterns[synth_obs].dist_snapped[0]
        self.assertEqual(snap_dist, known_dist)

        # network instantiated from a single vertical (down) libpysal.cg.Chain
        chain = cg.Chain([cg.Point((5, 5)), cg.Point((5, 4))])
        known_dist = 1.5
        self.ntw_from_chain = self.spaghetti.Network(in_data=chain)
        self.ntw_from_chain.snapobservations(cg.Point((6.5, 4.5)), synth_obs)
        snap_dist = self.ntw_from_chain.pointpatterns[synth_obs].dist_snapped[0]
        self.assertEqual(snap_dist, known_dist)

    def test_pp_failures(self):
        # network instantiated from a single libpysal.cg.Chain
        self.ntw_from_chain = self.spaghetti.Network(in_data=P11P22_CHAIN)
        # try snapping chain
        with self.assertRaises(TypeError):
            self.ntw_from_chain.snapobservations(P11P22_CHAIN, "chain")
        # try snapping list of chain
        with self.assertRaises(TypeError):
            self.ntw_from_chain.snapobservations([P11P22_CHAIN], "chain")

    @unittest.skipIf(GEOPANDAS_EXTINCT, "Missing Geopandas")
    def test_pp_from_geopandas(self):
        idxs = ["gdf_%s" % pp for pp in self.idxs]
        iterator = zip(self.obs, self.OBS, idxs)
        for (obs, OBS, idx) in iterator:
            OBS = geopandas.read_file(OBS)
            kwargs = {"attribute": True}
            if obs == crimes:
                kwargs.update({"idvariable": "POLYID"})
            self.ntw.snapobservations(OBS, obs, **kwargs)
            setattr(self, idx, self.ntw.pointpatterns[obs])

        self.assertEqual(self.pp1.npoints, self.gdf_pp1.npoints)
        self.assertEqual(self.pp2.npoints, self.gdf_pp2.npoints)

    def test_split_arcs_1000(self):
        n1000 = self.ntw.split_arcs(1000.0)
        self.assertEqual(len(n1000.arcs), 303)

    def test_add_point_pattern(self):
        self.assertEqual(self.pp1.npoints, self.known_pp1_npoints)
        self.assertIn("properties", self.pp1.points[0])
        self.assertIn([1], self.pp1.points[0]["properties"])

    def test_count_per_link_network(self):
        counts = self.ntw.count_per_link(self.pp1.obs_to_arc, graph=False)
        meancounts = sum(counts.values()) / float(len(counts.keys()))
        self.assertAlmostEqual(meancounts, 1.0, places=5)

    def test_count_per_edge_graph(self):
        counts = self.ntw.count_per_link(self.pp1.obs_to_arc, graph=True)
        meancounts = sum(counts.values()) / float(len(counts.keys()))
        self.assertAlmostEqual(meancounts, 1.0, places=5)

    def test_simulate_uniform_observations(self):
        sim = self.ntw.simulate_observations(self.known_pp1_npoints)
        self.assertEqual(self.known_pp1_npoints, sim.npoints)

    def test_simulate_unsupported_distribution_observations(self):
        with self.assertRaises(RuntimeError):
            self.ntw.simulate_observations(1, distribution="mrofinu")

    def test_all_neighbor_distances(self):
        matrix1, tree = self.ntw.allneighbordistances(schools, gen_tree=True)
        known_mtx_val = 17682.436988
        known_tree_val = (173, 64)

        self.assertAlmostEqual(numpy.nansum(matrix1[0]), known_mtx_val, places=4)
        self.assertEqual(tree[(6, 7)], known_tree_val)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        matrix2 = self.ntw.allneighbordistances(schools, fill_diagonal=0.0)
        observed = matrix2.diagonal()
        known = numpy.zeros(matrix2.shape[0])
        self.assertEqual(observed.all(), known.all())
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        matrix3 = self.ntw.allneighbordistances(schools, snap_dist=True)
        known_mtx_val = 3218.2597894
        observed_mtx_val = matrix3
        self.assertAlmostEqual(observed_mtx_val[0, 1], known_mtx_val, places=4)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        matrix4 = self.ntw.allneighbordistances(schools, fill_diagonal=0.0)
        observed = matrix4.diagonal()
        known = numpy.zeros(matrix4.shape[0])
        self.assertEqual(observed.all(), known.all())
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        matrix5, tree = self.ntw.allneighbordistances(crimes, gen_tree=True)
        known_mtx_val = 1484112.694526529
        known_tree_val = (-0.1, -0.1)

        self.assertAlmostEqual(numpy.nansum(matrix5[0]), known_mtx_val, places=4)
        self.assertEqual(tree[(18, 19)], known_tree_val)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

    def test_all_neighbor_distances_multiproccessing(self):
        matrix1, tree = self.ntw.allneighbordistances(
            schools, fill_diagonal=0.0, n_processes="all", gen_tree=True
        )
        known_mtx1_val = 17682.436988
        known_tree_val = (173, 64)

        observed = matrix1.diagonal()
        known = numpy.zeros(matrix1.shape[0])
        self.assertEqual(observed.all(), known.all())
        self.assertAlmostEqual(numpy.nansum(matrix1[0]), known_mtx1_val, places=4)
        self.assertEqual(tree[(6, 7)], known_tree_val)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        matrix2 = self.ntw.allneighbordistances(schools, n_processes=2)
        known_mtx2_val = 17682.436988
        self.assertAlmostEqual(numpy.nansum(matrix2[0]), known_mtx2_val, places=4)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        matrix3, tree = self.ntw.allneighbordistances(
            schools, fill_diagonal=0.0, n_processes=2, gen_tree=True
        )
        known_mtx3_val = 17682.436988
        known_tree_val = (173, 64)

        self.assertAlmostEqual(numpy.nansum(matrix3[0]), known_mtx3_val, places=4)
        self.assertEqual(tree[(6, 7)], known_tree_val)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

    def test_nearest_neighbor_distances(self):
        # general test
        with self.assertRaises(KeyError):
            self.ntw.nearestneighbordistances("i_should_not_exist")
        nnd1 = self.ntw.nearestneighbordistances(schools)
        nnd2 = self.ntw.nearestneighbordistances(schools, destpattern=schools)
        nndv1 = numpy.array(list(nnd1.values()), dtype=object)[:, 1].astype(float)
        nndv2 = numpy.array(list(nnd2.values()), dtype=object)[:, 1].astype(float)
        numpy.testing.assert_array_almost_equal_nulp(nndv1, nndv2)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        # nearest neighbor keeping zero test
        known_zero = ([19], 0.0)[0]
        nn_c = self.ntw.nearestneighbordistances(crimes, keep_zero_dist=True)
        self.assertEqual(nn_c[18][0], known_zero)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        # nearest neighbor omitting zero test
        known_nonzero = ([11], 165.33982412719126)[1]
        nn_c = self.ntw.nearestneighbordistances(crimes, keep_zero_dist=False)
        self.assertAlmostEqual(nn_c[18][1], known_nonzero, places=4)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        # nearest neighbor with snap distance
        known_neigh = ([3], 402.5219673922477)[1]
        nn_c = self.ntw.nearestneighbordistances(crimes, snap_dist=True)
        self.assertAlmostEqual(nn_c[0][1], known_neigh, places=4)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

    @unittest.skipIf(GEOPANDAS_EXTINCT, "Missing Geopandas")
    def test_element_as_gdf(self):
        obs = self.spaghetti.element_as_gdf(self.ntw, pp_name=schools)
        snap_obs = self.spaghetti.element_as_gdf(
            self.ntw, pp_name=schools, snapped=True
        )

        known_dist = 205.65961300587043
        observed_point = obs.loc[(obs["id"] == 0), "geometry"].squeeze()
        snap_point = snap_obs.loc[(snap_obs["id"] == 0), "geometry"].squeeze()
        observed_dist = observed_point.distance(snap_point)
        self.assertAlmostEqual(observed_dist, known_dist, places=8)

        with self.assertRaises(KeyError):
            self.spaghetti.element_as_gdf(self.ntw, pp_name="i_should_not_exist")


# -------------------------------------------------------------------------------
class TestNetworkAnalysis(unittest.TestCase):
    def setUp(self):
        bounds, h, v = (0, 0, 3, 3), 2, 2
        lattice = self.spaghetti.regular_lattice(bounds, h, nv=v, exterior=True)
        self.ntw = self.spaghetti.Network(in_data=lattice)
        chains = [
            cg.Chain(
                [
                    cg.Point(self.ntw.vertex_coords[p1]),
                    cg.Point(self.ntw.vertex_coords[p2]),
                ]
            )
            for (p1, p2) in self.ntw.arcs
        ]
        midpoints = []
        for chain in chains:
            (v1x, v1y), (v2x, v2y) = chain.vertices
            mid = cg.Point(((v1x + v2x) / 2.0, (v1y + v2y) / 2.0))
            midpoints.append(mid)
        self.mids = "mids"
        self.ntw.snapobservations(midpoints, self.mids)
        npts = self.ntw.pointpatterns[self.mids].npoints
        self.test_permutations = 99
        self.test_steps = 10

    def tearDown(self):
        pass

    def test_global_auto_k_uniform(self):
        known_lowerenvelope = numpy.array(
            [
                0.0,
                0.1875,
                0.60416667,
                1.22916667,
                2.10416667,
                2.85416667,
                3.70833333,
                4.70833333,
                5.1875,
                5.41666667,
            ]
        )
        numpy.random.seed(0)
        obtained = self.ntw.GlobalAutoK(
            self.ntw.pointpatterns[self.mids],
            permutations=self.test_permutations,
            nsteps=self.test_steps,
            distribution="uniform",
        )
        self.assertEqual(obtained.lowerenvelope.shape[0], self.test_steps)
        numpy.testing.assert_allclose(obtained.lowerenvelope, known_lowerenvelope)

    def test_global_auto_k_unsupported_distribution(self):
        with self.assertRaises(RuntimeError):
            self.ntw.GlobalAutoK(
                self.ntw.pointpatterns[self.mids],
                permutations=self.test_permutations,
                nsteps=self.test_steps,
                distribution="mrofinu",
            )


# -------------------------------------------------------------------------------
class TestNetworkUtils(unittest.TestCase):
    def setUp(self):
        self.ntw = self.spaghetti.Network(in_data=STREETS)
        self.P00, self.P01 = P00, P01
        self.P10, self.P11 = P10, P11

    def tearDown(self):
        pass

    def test_compute_length(self):
        point1, point2 = (0, 0), (1, 1)
        known_length = 1.4142135623730951
        observed_length = self.util.compute_length(point1, point2)
        self.assertAlmostEqual(observed_length, known_length, places=4)

    def test_get_neighbor_distances(self):
        known_neighs = {1: 102.62353453439829, 2: 660.000001049743}
        observed_neighs = self.util.get_neighbor_distances(
            self.ntw, 0, self.ntw.arc_lengths
        )
        for k in known_neighs.keys():
            self.assertAlmostEqual(observed_neighs[k], known_neighs[k], places=4)

    def test_generate_tree(self):
        known_path = [23, 22, 20, 19, 170, 2, 0]
        distance, pred = self.util.dijkstra(self.ntw, 0)
        tree = self.util.generatetree(pred)
        self.assertEqual(tree[3], known_path)

    def test_dijkstra(self):
        distance, pred = self.util.dijkstra(self.ntw, 0)
        self.assertAlmostEqual(distance[196], 5505.668247, places=4)
        self.assertEqual(pred[196], 133)

    def test_dijkstra_mp(self):
        distance, pred = self.util.dijkstra_mp((self.ntw, 0))
        self.assertAlmostEqual(distance[196], 5505.668247, places=4)
        self.assertEqual(pred[196], 133)

    def test_chain_constr(self):
        known_len = 1.4142135623730951
        chain = self.util.chain_constr({0: (0.0, 0.0), 1: (1.0, 1.0)}, [(0, 1)])[0]
        self.assertAlmostEqual(chain.len, known_len, places=10)

    def test_build_chains(self):
        # 1x1 cross (regular lattice without the exterior)
        v_known = [self.P10, self.P11]
        h_known = [self.P01, self.P11]
        space_h = space_v = [0.0, 1.0, 2.0]
        exterior, bounds = False, (0, 0, 2, 2)
        _vl = self.util.build_chains(space_h, space_v, exterior, bounds, h=False)
        _hl = self.util.build_chains(space_h, space_v, exterior, bounds, h=True)
        v_observed = _vl[0].vertices
        h_observed = _hl[0].vertices
        self.assertEqual(v_observed, v_known)
        self.assertEqual(h_observed, h_known)

        # 3x3 cross (regular lattice with the exterior)
        v_known = [self.P00, self.P01]
        h_known = [self.P00, self.P10]
        exterior, bounds = True, (0, 0, 2, 2)
        _vl = self.util.build_chains(space_h, space_v, exterior, bounds, h=False)
        _hl = self.util.build_chains(space_h, space_v, exterior, bounds, h=True)
        v_observed = _vl[0].vertices
        h_observed = _hl[0].vertices
        self.assertEqual(v_observed, v_known)
        self.assertEqual(h_observed, h_known)

    def test_squared_distance_point_link(self):
        point, link = (1, 1), ((0, 0), (2, 0))
        sqrd_nearp = self.util.squared_distance_point_link(point, link)
        self.assertEqual(sqrd_nearp[0], 1.0)
        self.assertEqual(sqrd_nearp[1].all(), numpy.array([1.0, 0.0]).all())

    def test_snap_points_to_links(self):
        points = {0: P11}
        links = [cg.shapes.Chain([P00, cg.shapes.Point((2, 0))])]
        snapped = self.util.snap_points_to_links(points, links)
        known_coords = [xy._Point__loc for xy in snapped[0][0]]
        self.assertEqual(known_coords, [(0.0, 0.0), (2.0, 0.0)])
        self.assertEqual(snapped[0][1].all(), numpy.array([1.0, 0.0]).all())
