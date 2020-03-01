from libpysal import cg, examples, io
import numpy
import unittest
import copy

try:
    import geopandas

    GEOPANDAS_EXTINCT = False
except ImportError:
    GEOPANDAS_EXTINCT = True


# empirical network shapefile
STREETS = examples.get_path("streets.shp")

# native pysal geometries ------------------------------------------------------
# single libpysal.cg.Chain
P00 = cg.Point((0, 0))
P11 = cg.Point((1, 1))
P22 = cg.Point((2, 2))
P11P22_CHAIN = cg.Chain([P11, P22])


class TestNetwork(unittest.TestCase):
    def setUp(self):

        # empirical network instantiated from shapefile
        self.ntw_shp = self.spaghetti.Network(in_data=STREETS, weightings=True)
        self.n_known_shp_arcs, self.n_known_shp_vertices = 303, 230

        # native pysal geometries
        self.chains_from_shp = chains = [
            cg.Chain([cg.Point(self.ntw_shp.vertex_coords[vertex]) for vertex in arc])
            for arc in self.ntw_shp.arcs
        ]

        # synthetic network instantiated from pysal geometries

        # lattice and ring + extension
        bounds = (0, 0, 2, 2)
        h, v = 1, 1
        self.lattice = self.spaghetti.regular_lattice(bounds, h, nv=v, exterior=False)
        self.ring = [
            cg.Chain([cg.Point([3, 1]), cg.Point([3.25, 1.25])]),
            cg.Chain([cg.Point([3.25, 1.25]), cg.Point([3.375, 1.25])]),
            cg.Chain([cg.Point([3.375, 1.25]), cg.Point([3.5, 1.25])]),
            cg.Chain([cg.Point([3.5, 1.25]), cg.Point([3.75, 1])]),
            cg.Chain([cg.Point([3.75, 1]), cg.Point([3.5, 0.75])]),
            cg.Chain([cg.Point([3.5, 0.75]), cg.Point([3.375, 0.75])]),
            cg.Chain([cg.Point([3.375, 0.75]), cg.Point([3.25, 0.75])]),
            cg.Chain([cg.Point([3.25, 0.75]), cg.Point([3, 1])]),
        ]
        self.extension = [
            cg.Chain([cg.Point([1, 2]), cg.Point([2, 2]), cg.Point([2, 1])])
        ]

        self.lines = self.lattice + self.ring + self.extension
        self.ntw_from_lattice_ring = self.spaghetti.Network(in_data=self.lines)

        self.p0505 = cg.Point([0.5, 0.5])
        self.p052 = cg.Point([0.5, 2.0])
        self.p0525 = cg.Point([0.5, 2.5])

        self.ntw_from_lattice_ring.snapobservations(
            [self.p0505, self.p052], "p", attribute=False
        )

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
            self.assertEqual(
                self.ntw_from_chains.network_n_components, known_components
            )
            self.assertAlmostEqual(
                sum(self.ntw_from_chains.arc_lengths.values()), known_length, places=3
            )

    def test_network_from_single_libpysal_chain(self):
        ## network instantiated from a single libpysal.cg.Chain
        self.ntw_chain_out = self.spaghetti.Network(in_data=P11P22_CHAIN)
        # test save and load network
        self.ntw_chain_out.savenetwork("test_network.pkl")
        self.ntw_chain_in = self.spaghetti.Network.loadnetwork("test_network.pkl")
        self.assertEqual(self.ntw_chain_in.arcs, self.ntw_chain_in.edges)

    def test_network_from_vertical_libpysal_chains(self):
        vert_up = cg.Chain(
            [self.p0505, self.p052]
        )  ################################################################### test vertical stuff...
        self.ntw_up_chain = self.spaghetti.Network(in_data=vert_up)
        self.assertEqual(len(self.ntw_up_chain.arcs), len(vert_up.segments))

        vert_down = cg.Chain(
            [self.p052, self.p0505]
        )  ###################################################################### same
        self.ntw_down_chain = self.spaghetti.Network(in_data=vert_down)
        self.assertEqual(len(self.ntw_down_chain.arcs), len(vert_down.segments))

    def test_network_failures(self):
        # try instantiating network with single point
        with self.assertRaises(TypeError):
            self.spaghetti.Network(in_data=self.p11)
        # try instantiating network with list of single point
        with self.assertRaises(TypeError):
            self.spaghetti.Network(in_data=[self.p11])

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
        self.ntw_shp.snapobservations(examples.get_path("schools.shp"), "schools")
        _, tree = self.ntw_shp.allneighbordistances("schools", gen_tree=True)
        observed_paths = self.ntw_shp.shortest_paths(tree, "schools")
        observed_vertices = len(observed_paths[0][1].vertices)
        self.assertEqual(observed_vertices, known_vertices)

        # asymmetric point pattern
        bounds, h, v = (0, 0, 3, 3), 2, 2
        lattice = self.spaghetti.regular_lattice(bounds, h, nv=v, exterior=False)
        ntw = self.spaghetti.Network(in_data=lattice)

        self.p0505, self.p052

        points1 = [
            self.p0505,
            cg.Point((2.5, 2.5)),
        ]  ##############################################
        points2 = [self.p0525, cg.Point((2.5, 0.5)), cg.Point((0.75, 0.6))]
        ntw.snapobservations(points1, "points1")
        ntw.snapobservations(points2, "points2")
        _, tree = ntw.allneighbordistances("points1", "points2", gen_tree=True)
        observed_paths = ntw.shortest_paths(tree, "points1", pp_dest="points2")

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
            paths = ntw.shortest_paths([], "synth_obs")

    @unittest.skipIf(GEOPANDAS_EXTINCT, "Missing Geopandas")
    def test_element_as_gdf(self):
        vertices, arcs = self.spaghetti.element_as_gdf(
            self.ntw_shp, vertices=True, arcs=True
        )

        known_vertex_wkt = "POINT (728368.04762 877125.89535)"
        obs_vertex = vertices.loc[(vertices["id"] == 0), "geometry"].squeeze()
        obs_vertex_wkt = obs_vertex.wkt
        self.assertEqual(obs_vertex_wkt, known_vertex_wkt)

        known_arc_wkt = (
            "LINESTRING (728368.04762 877125.89535, " + "728368.13931 877023.27186)"
        )
        obs_arc = arcs.loc[(arcs["id"] == (0, 1)), "geometry"].squeeze()
        obs_arc_wkt = obs_arc.wkt
        self.assertEqual(obs_arc_wkt, known_arc_wkt)

        arcs = self.spaghetti.element_as_gdf(self.ntw_shp, arcs=True)
        known_arc_wkt = (
            "LINESTRING (728368.04762 877125.89535, " + "728368.13931 877023.27186)"
        )
        obs_arc = arcs.loc[(arcs["id"] == (0, 1)), "geometry"].squeeze()
        obs_arc_wkt = obs_arc.wkt
        self.assertEqual(obs_arc_wkt, known_arc_wkt)

        # symmetric routes
        known_length, bounds, h, v = 2.6, (0, 0, 3, 3), 2, 2
        lattice = self.spaghetti.regular_lattice(bounds, h, nv=v, exterior=False)
        ntw = self.spaghetti.Network(in_data=lattice)
        synth_obs = [cg.Point([0.2, 1.3]), cg.Point([0.2, 1.7]), cg.Point([2.8, 1.5])]
        ntw.snapobservations(synth_obs, "synth_obs")
        _, tree = ntw.allneighbordistances("synth_obs", gen_tree=True)
        paths = ntw.shortest_paths(tree, "synth_obs")
        paths_gdf = self.spaghetti.element_as_gdf(ntw, routes=paths)
        observed_length = paths_gdf.loc[0, "geometry"].length
        self.assertEqual(observed_length, known_length)

        # asymmetric routes
        known_origins, bounds, h, v = 2, (0, 0, 3, 3), 2, 2
        lattice = self.spaghetti.regular_lattice(bounds, h, nv=v, exterior=False)
        ntw = self.spaghetti.Network(in_data=lattice)
        points1 = [
            cg.Point((0.5, 0.5)),
            cg.Point((2.5, 2.5)),
        ]  ######################################
        points2 = [
            cg.Point((0.5, 2.5)),
            cg.Point((2.5, 0.5)),
        ]  ####################################
        ntw.snapobservations(points1, "points1")
        ntw.snapobservations(points2, "points2")
        _, tree = ntw.allneighbordistances("points1", "points2", gen_tree=True)
        paths = ntw.shortest_paths(tree, "points1", pp_dest="points2")
        paths_gdf = self.spaghetti.element_as_gdf(ntw, routes=paths)
        observed_origins = paths_gdf["O"].nunique()
        self.assertEqual(observed_origins, known_origins)

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

    def test_extract_component(self):

        ntw = copy.deepcopy(self.ntw_from_lattice_ring)

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

    def test_regular_lattice(self):
        # 4x4 regular lattice with the exterior
        known = [
            cg.Point((0.0, 0.0)),
            cg.Point((1.0, 0.0)),
        ]  ####################################
        bounds = (0, 0, 3, 3)
        lattice = self.spaghetti.regular_lattice(bounds, 2, nv=2, exterior=True)
        observed = lattice[0].vertices
        self.assertEqual(observed, known)

        # 5x5 regular lattice without the exterior
        known = [
            cg.Point((3.0, 3.0)),
            cg.Point((3.0, 4.0)),
        ]  #####################################
        bounds = (0, 0, 4, 4)
        lattice = self.spaghetti.regular_lattice(bounds, 3, exterior=False)
        observed = lattice[-1].vertices
        self.assertEqual(observed, known)

        # 7x9 regular lattice from shapefile bounds
        known_vertices = [
            (723414.3683108028, 875929.0396895551),
            (724286.1381211297, 875929.0396895551),
        ]
        shp = io.open(
            STREETS
        )  ####################################################################
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


class TestNetworkPointPattern(unittest.TestCase):
    def setUp(self):
        # path_to_shp = examples.get_path("streets.shp")###################################
        self.ntw = self.spaghetti.Network(in_data=STREETS)
        self.pp1_str = "schools"
        self.pp2_str = "crimes"
        iterator = [(self.pp1_str, "pp1"), (self.pp2_str, "pp2")]
        for (obs, idx) in iterator:
            path_to_shp = examples.get_path("%s.shp" % obs)
            self.ntw.snapobservations(path_to_shp, obs, attribute=True)
            setattr(self, idx, self.ntw.pointpatterns[obs])
        self.known_pp1_npoints = 8

    def tearDown(self):
        pass

    def test_pp_from_libpysal_points(self):
        # known
        crimes = self.ntw.pointpatterns["crimes"]
        known_snapped = set(crimes.snapped_coordinates.values())
        # points from pysal geometries
        points = [cg.Point(crimes.points[i]["coordinates"]) for i in crimes.points]
        for dtype in (list, tuple):
            point_data = dtype(points)
            self.ntw.snapobservations(point_data, "cg_crimes")
            observed = self.ntw.pointpatterns["cg_crimes"]
            observed_snapped = set(observed.snapped_coordinates.values())
            self.assertEqual(observed_snapped, known_snapped)

    def test_pp_from_single_libpysal_point(self):
        # network instantiated from a single libpysal.cg.Chain
        known_dist = 1.4142135623730951
        self.ntw_from_chain = self.spaghetti.Network(in_data=P11P22_CHAIN)
        self.ntw_from_chain.snapobservations(cg.Point((0, 0)), "synth_obs")
        snap_dist = self.ntw_from_chain.pointpatterns["synth_obs"].dist_snapped[0]
        self.assertAlmostEqual(snap_dist, known_dist, places=10)

        # network instantiated from a single vertical (up) libpysal.cg.Chain
        chain = cg.Chain([cg.Point((1, 1)), cg.Point((1, 2))])
        known_dist = 1.0
        self.ntw_from_chain = self.spaghetti.Network(in_data=chain)
        self.ntw_from_chain.snapobservations(cg.Point((0, 1.5)), "synth_obs")
        snap_dist = self.ntw_from_chain.pointpatterns["synth_obs"].dist_snapped[0]
        self.assertEqual(snap_dist, known_dist)

        # network instantiated from a single vertical (down) libpysal.cg.Chain
        chain = cg.Chain([cg.Point((5, 5)), cg.Point((5, 4))])
        known_dist = 1.5
        self.ntw_from_chain = self.spaghetti.Network(in_data=chain)
        self.ntw_from_chain.snapobservations(cg.Point((6.5, 4.5)), "synth_obs")
        snap_dist = self.ntw_from_chain.pointpatterns["synth_obs"].dist_snapped[0]
        self.assertEqual(snap_dist, known_dist)

    def test_pp_failures(self):
        # network instantiated from a single libpysal.cg.Chain
        self.ntw_from_chain = self.spaghetti.Network(in_data=P11P22_CHAIN)
        # try snapping chain
        with self.assertRaises(TypeError):
            self.ntw_from_chain.snapobservations(chain, "chain")
        # try snapping list of chain
        with self.assertRaises(TypeError):
            self.ntw_from_chain.snapobservations([chain], "chain")

    @unittest.skipIf(GEOPANDAS_EXTINCT, "Missing Geopandas")
    def test_pp_from_geopandas(self):
        self.gdf_pp1_str = "schools"
        self.gdf_pp2_str = "crimes"
        iterator = [(self.gdf_pp1_str, "gdf_pp1"), (self.gdf_pp2_str, "gdf_pp2")]
        for (obs, idx) in iterator:
            path_to_shp = examples.get_path("%s.shp" % obs)
            in_data = geopandas.read_file(path_to_shp)
            self.ntw.snapobservations(in_data, obs, attribute=True)
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

    def test_simulate_normal_observations(self):
        sim = self.ntw.simulate_observations(self.known_pp1_npoints)
        self.assertEqual(self.known_pp1_npoints, sim.npoints)

    def test_simulate_poisson_observations(self):
        sim = self.ntw.simulate_observations(
            self.known_pp1_npoints, distribution="poisson"
        )
        self.assertEqual(self.known_pp1_npoints, sim.npoints)

    def test_all_neighbor_distances(self):
        matrix1, tree = self.ntw.allneighbordistances(self.pp1_str, gen_tree=True)
        known_mtx_val = 17682.436988
        known_tree_val = (173, 64)

        self.assertAlmostEqual(numpy.nansum(matrix1[0]), known_mtx_val, places=4)
        self.assertEqual(tree[(6, 7)], known_tree_val)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        matrix2 = self.ntw.allneighbordistances(self.pp1_str, fill_diagonal=0.0)
        observed = matrix2.diagonal()
        known = numpy.zeros(matrix2.shape[0])
        self.assertEqual(observed.all(), known.all())
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        matrix3 = self.ntw.allneighbordistances(self.pp1_str, snap_dist=True)
        known_mtx_val = 3218.2597894
        observed_mtx_val = matrix3
        self.assertAlmostEqual(observed_mtx_val[0, 1], known_mtx_val, places=4)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        matrix4 = self.ntw.allneighbordistances(self.pp1_str, fill_diagonal=0.0)
        observed = matrix4.diagonal()
        known = numpy.zeros(matrix4.shape[0])
        self.assertEqual(observed.all(), known.all())
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        matrix5, tree = self.ntw.allneighbordistances(self.pp2_str, gen_tree=True)
        known_mtx_val = 1484112.694526529
        known_tree_val = (-0.1, -0.1)

        self.assertAlmostEqual(numpy.nansum(matrix5[0]), known_mtx_val, places=4)
        self.assertEqual(tree[(18, 19)], known_tree_val)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

    def test_all_neighbor_distances_multiproccessing(self):
        matrix1, tree = self.ntw.allneighbordistances(
            self.pp1_str, fill_diagonal=0.0, n_processes="all", gen_tree=True
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

        matrix2 = self.ntw.allneighbordistances(self.pp1_str, n_processes=2)
        known_mtx2_val = 17682.436988
        self.assertAlmostEqual(numpy.nansum(matrix2[0]), known_mtx2_val, places=4)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        matrix3, tree = self.ntw.allneighbordistances(
            self.pp1_str, fill_diagonal=0.0, n_processes=2, gen_tree=True
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
        nnd1 = self.ntw.nearestneighbordistances(self.pp1_str)
        nnd2 = self.ntw.nearestneighbordistances(self.pp1_str, destpattern=self.pp1_str)
        nndv1 = numpy.array(list(nnd1.values()))[:, 1].astype(float)
        nndv2 = numpy.array(list(nnd2.values()))[:, 1].astype(float)
        numpy.testing.assert_array_almost_equal_nulp(nndv1, nndv2)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        # nearest neighbor keeping zero test
        known_zero = ([19], 0.0)[0]
        nn_c = self.ntw.nearestneighbordistances(self.pp2_str, keep_zero_dist=True)
        self.assertEqual(nn_c[18][0], known_zero)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        # nearest neighbor omitting zero test
        known_nonzero = ([11], 165.33982412719126)[1]
        nn_c = self.ntw.nearestneighbordistances(self.pp2_str, keep_zero_dist=False)
        self.assertAlmostEqual(nn_c[18][1], known_nonzero, places=4)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

        # nearest neighbor with snap distance
        known_neigh = ([3], 402.5219673922477)[1]
        nn_c = self.ntw.nearestneighbordistances(
            self.pp2_str, keep_zero_dist=True, snap_dist=True
        )
        self.assertAlmostEqual(nn_c[0][1], known_neigh, places=4)
        del self.ntw.distance_matrix
        del self.ntw.network_trees

    @unittest.skipIf(GEOPANDAS_EXTINCT, "Missing Geopandas")
    def test_element_as_gdf(self):
        obs = self.spaghetti.element_as_gdf(self.ntw, pp_name=self.pp1_str)
        snap_obs = self.spaghetti.element_as_gdf(
            self.ntw, pp_name=self.pp1_str, snapped=True
        )

        known_dist = 205.65961300587043
        observed_point = obs.loc[(obs["id"] == 0), "geometry"].squeeze()
        snap_point = snap_obs.loc[(snap_obs["id"] == 0), "geometry"].squeeze()
        observed_dist = observed_point.distance(snap_point)
        self.assertAlmostEqual(observed_dist, known_dist, places=8)

        with self.assertRaises(KeyError):
            self.spaghetti.element_as_gdf(self.ntw, pp_name="i_should_not_exist")


class TestNetworkAnalysis(unittest.TestCase):
    def setUp(self):
        # path_to_shp = examples.get_path("streets.shp")#####################################
        self.ntw = self.spaghetti.Network(in_data=STREETS)
        self.pt_str = "schools"
        path_to_shp = examples.get_path("%s.shp" % self.pt_str)
        self.ntw.snapobservations(path_to_shp, self.pt_str, attribute=True)
        npts = self.ntw.pointpatterns[self.pt_str].npoints
        self.ntw.simulate_observations(npts)
        self.test_permutations = 3
        self.test_steps = 5

    def tearDown(self):
        pass

    def test_network_f(self):
        obtained = self.ntw.NetworkF(
            self.ntw.pointpatterns[self.pt_str],
            permutations=self.test_permutations,
            nsteps=self.test_steps,
        )
        self.assertEqual(obtained.lowerenvelope.shape[0], self.test_steps)

    def test_network_g(self):
        obtained = self.ntw.NetworkG(
            self.ntw.pointpatterns[self.pt_str],
            permutations=self.test_permutations,
            nsteps=self.test_steps,
        )
        self.assertEqual(obtained.lowerenvelope.shape[0], self.test_steps)

    def test_network_k(self):
        obtained = self.ntw.NetworkK(
            self.ntw.pointpatterns[self.pt_str],
            permutations=self.test_permutations,
            nsteps=self.test_steps,
        )
        self.assertEqual(obtained.lowerenvelope.shape[0], self.test_steps)


class TestNetworkUtils(unittest.TestCase):
    def setUp(self):
        # path_to_shp = examples.get_path("streets.shp")#####################################
        self.ntw = self.spaghetti.Network(in_data=STREETS)

    def tearDown(self):
        pass

    def test_compute_length(self):
        self.point1, self.point2 = (0, 0), (1, 1)
        self.length = self.util.compute_length(self.point1, self.point2)
        self.assertAlmostEqual(self.length, 1.4142135623730951, places=4)

    def test_get_neighbor_distances(self):
        self.known_neighs = {1: 102.62353453439829, 2: 660.000001049743}
        self.neighs = self.util.get_neighbor_distances(
            self.ntw, 0, self.ntw.arc_lengths
        )
        self.assertAlmostEqual(self.neighs[1], 102.62353453439829, places=4)
        self.assertAlmostEqual(self.neighs[2], 660.000001049743, places=4)

    def test_generate_tree(self):
        self.known_path = [23, 22, 20, 19, 170, 2, 0]
        self.distance, self.pred = self.util.dijkstra(self.ntw, 0)
        self.tree = self.util.generatetree(self.pred)
        self.assertEqual(self.tree[3], self.known_path)

    def test_dijkstra(self):
        self.distance, self.pred = self.util.dijkstra(self.ntw, 0)
        self.assertAlmostEqual(self.distance[196], 5505.668247, places=4)
        self.assertEqual(self.pred[196], 133)

    def test_dijkstra_mp(self):
        self.distance, self.pred = self.util.dijkstra_mp((self.ntw, 0))
        self.assertAlmostEqual(self.distance[196], 5505.668247, places=4)
        self.assertEqual(self.pred[196], 133)

    def test_squared_distance_point_link(self):
        self.point, self.link = (1, 1), ((0, 0), (2, 0))
        self.sqrd_nearp = self.util.squared_distance_point_link(self.point, self.link)
        self.assertEqual(self.sqrd_nearp[0], 1.0)
        self.assertEqual(self.sqrd_nearp[1].all(), numpy.array([1.0, 0.0]).all())

    def test_snap_points_to_links(self):
        self.points = {0: cg.shapes.Point((1, 1))}
        self.links = [
            cg.shapes.Chain([cg.shapes.Point((0, 0)), cg.shapes.Point((2, 0))])
        ]
        self.snapped = self.util.snap_points_to_links(self.points, self.links)
        self.known_coords = [xy._Point__loc for xy in self.snapped[0][0]]
        self.assertEqual(self.known_coords, [(0.0, 0.0), (2.0, 0.0)])
        self.assertEqual(self.snapped[0][1].all(), numpy.array([1.0, 0.0]).all())
