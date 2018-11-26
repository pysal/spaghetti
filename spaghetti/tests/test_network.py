import unittest
import numpy as np
from libpysal import cg, examples

# dev import structure
from .. import network as spgh
from .. import util

try:
    import geopandas
    GEOPANDAS_EXTINCT = False
except ImportError:
    GEOPANDAS_EXTINCT = True


@unittest.skipIf(GEOPANDAS_EXTINCT, 'Missing Geopandas')
class TestNetwork(unittest.TestCase):
    
    def setUp(self):
        path_to_shp = examples.get_path('streets.shp')
        
        # network instantiated from shapefile
        self.ntw_from_shp = spgh.Network(in_data=path_to_shp)
        
        # network instantiated from geodataframe
        gdf = geopandas.read_file(path_to_shp)
        self.ntw_from_gdf = spgh.Network(in_data=gdf)
    
    def tearDown(self):
        pass
    
    def test_network_data_read(self):
        n_known_edges,  n_known_nodes= 303, 230
        
        # shp test against known
        self.assertEqual(len(self.ntw_from_shp.edges), n_known_edges)
        self.assertEqual(len(self.ntw_from_shp.nodes), n_known_nodes)
        # gdf test against known
        self.assertEqual(len(self.ntw_from_gdf.edges), n_known_edges)
        self.assertEqual(len(self.ntw_from_gdf.nodes), n_known_nodes)
        # shp against gdf
        self.assertEqual(len(self.ntw_from_shp.edges),
                         len(self.ntw_from_gdf.edges))
        self.assertEqual(len(self.ntw_from_shp.nodes),
                         len(self.ntw_from_gdf.nodes))
    
    def test_extract_network(self):
        self.assertEqual(len(self.ntw_from_shp.edges), 303)
        self.assertEqual(len(self.ntw_from_shp.nodes), 230)
        edgelengths = self.ntw_from_shp.edge_lengths.values()
        self.assertAlmostEqual(sum(edgelengths), 104414.0920159, places=5)
        self.assertIn(0, self.ntw_from_shp.adjacencylist[1])
        self.assertIn(0, self.ntw_from_shp.adjacencylist[2])
        self.assertNotIn(0, self.ntw_from_shp.adjacencylist[3])
    
    def test_contiguity_weights(self):
        w = self.ntw_from_shp.contiguityweights(graph=False)
        self.assertEqual(w.n, 303)
        self.assertEqual(w.histogram,
                         [(2, 35), (3, 89), (4, 105), (5, 61), (6, 13)])
    
    def test_contiguity_weights_graph(self):
        w = self.ntw_from_shp.contiguityweights(graph=True)
        self.assertEqual(w.n, 179)
        self.assertEqual(w.histogram,
                         [(2, 2), (3, 2), (4, 45), (5, 82), (6, 48)])
    
    def test_distance_band_weights(self):
        # I do not trust this result, test should be manually checked.
        w = self.ntw_from_shp.distancebandweights(threshold=500)
        self.assertEqual(w.n, 230)
        self.assertEqual(w.histogram,
                         [(1, 22), (2, 58), (3, 63), (4, 40),
                          (5, 36), (6, 3), (7, 5), (8, 3)])
    
    def test_edge_segmentation(self):
        n200 = self.ntw_from_shp.segment_edges(200.0)
        self.assertEqual(len(n200.edges), 688)
        n200 = None
    
    def test_enum_links_node(self):
        coincident = self.ntw_from_shp.enum_links_node(24)
        self.assertIn((24, 48), coincident)
    
    def test_element_as_gdf(self):
        nodes, edges = spgh.element_as_gdf(self.ntw_from_shp,
                                              nodes=True,
                                              edges=True)
        
        known_node_wkt = 'POINT (728368.04762 877125.89535)'
        obs_node = nodes.loc[(nodes['id'] == 0), 'geometry'].squeeze()
        obs_node_wkt = obs_node.wkt
        self.assertEqual(obs_node_wkt, known_node_wkt)
        
        known_edge_wkt = 'LINESTRING (728368.04762 877125.89535, '\
                         + '728368.13931 877023.27186)'
        obs_edge = edges.loc[(edges['id'] == (0,1)), 'geometry'].squeeze()
        obs_edge_wkt = obs_edge.wkt
        self.assertEqual(obs_edge_wkt, known_edge_wkt)
        
        edges = spgh.element_as_gdf(self.ntw_from_shp, edges=True)
        known_edge_wkt = 'LINESTRING (728368.04762 877125.89535, '\
                         + '728368.13931 877023.27186)'
        obs_edge = edges.loc[(edges['id'] == (0,1)), 'geometry'].squeeze()
        obs_edge_wkt = obs_edge.wkt
        self.assertEqual(obs_edge_wkt, known_edge_wkt)
        
    def test_round_sig(self):
        # round to 2 significant digits test
        x_round2, y_round2 = 1200, 1900
        self.ntw_from_shp.node_sig = 2
        obs_xy_round2 = self.ntw_from_shp._round_sig((1215, 1865))
        self.assertEqual(obs_xy_round2, (x_round2, y_round2))
        
        # round to no significant digits test
        x_roundNone, y_roundNone = 1215, 1865
        self.ntw_from_shp.node_sig = None
        obs_xy_roundNone = self.ntw_from_shp._round_sig((1215, 1865))
        self.assertEqual(obs_xy_roundNone, (x_roundNone, y_roundNone))


@unittest.skipIf(GEOPANDAS_EXTINCT, 'Missing Geopandas')
class TestNetworkPointPattern(unittest.TestCase):
    
    def setUp(self):
        path_to_shp = examples.get_path('streets.shp')
        self.ntw = spgh.Network(in_data=path_to_shp)
        self.pp1_str = 'schools'
        self.pp2_str = 'crimes'
        for (obs, idx) in [(self.pp1_str, 'pp1'), (self.pp2_str, 'pp2')]:
            path_to_shp = examples.get_path('%s.shp' % obs)
            in_data = geopandas.read_file(path_to_shp)
            self.ntw.snapobservations(in_data, obs, attribute=True)
            setattr(self, idx, self.ntw.pointpatterns[obs])
        self.known_pp1_npoints = 8
    
    def tearDown(self):
        pass
    
    def test_add_point_pattern(self):
        self.assertEqual(self.pp1.npoints, self.known_pp1_npoints)
        self.assertIn('properties', self.pp1.points[0])
        self.assertIn([1], self.pp1.points[0]['properties'])
    
    def test_count_per_edge(self):
        counts = self.ntw.count_per_edge(self.pp1.obs_to_edge, graph=False)
        meancounts = sum(counts.values()) / float(len(counts.keys()))
        self.assertAlmostEqual(meancounts, 1.0, places=5)
    
    def test_count_per_graph_edge(self):
        counts = self.ntw.count_per_edge(self.pp1.obs_to_edge, graph=True)
        meancounts = sum(counts.values()) / float(len(counts.keys()))
        self.assertAlmostEqual(meancounts, 1.0, places=5)
    
    def test_simulate_normal_observations(self):
        sim = self.ntw.simulate_observations(self.known_pp1_npoints)
        self.assertEqual(self.known_pp1_npoints, sim.npoints)
    
    def test_simulate_poisson_observations(self):
        sim = self.ntw.simulate_observations(self.known_pp1_npoints,
                                             distribution='poisson')
        self.assertEqual(self.known_pp1_npoints, sim.npoints)
    
    def test_all_neighbor_distances_mp(self):
        matrix1, tree = self.ntw.allneighbordistances(self.pp1_str,
                                                      fill_diagonal=0.,
                                                      n_processes='all',
                                                      gen_tree=True)
        known_mtx1_val = 17682.436988
        known_tree_val = (173, 64)
        
        observed = matrix1.diagonal()
        known = np.zeros(matrix1.shape[0])
        self.assertEqual(observed.all(), known.all())
        self.assertAlmostEqual(np.nansum(matrix1[0]), known_mtx1_val, places=4)
        self.assertEqual(tree[(6, 7)], known_tree_val)
        
        matrix2 = self.ntw.allneighbordistances(self.pp1_str, n_processes=2)
        known_mtx2_val = 17682.436988
        self.assertAlmostEqual(np.nansum(matrix2[0]), known_mtx2_val, places=4)
        
        matrix3, tree = self.ntw.allneighbordistances(self.pp1_str,
                                                      fill_diagonal=0.,
                                                      n_processes=2,
                                                      gen_tree=True)
        known_mtx3_val = 17682.436988
        known_tree_val = (173, 64)
        
        self.assertAlmostEqual(np.nansum(matrix3[0]), known_mtx3_val, places=4)
        self.assertEqual(tree[(6, 7)], known_tree_val)
    
    def test_all_neighbor_distances(self):
        matrix1, tree = self.ntw.allneighbordistances(self.pp1_str,
                                                      gen_tree=True)
        known_mtx_val = 17682.436988
        known_tree_val = (173, 64)
        
        self.assertAlmostEqual(np.nansum(matrix1[0]), known_mtx_val, places=4)
        self.assertEqual(tree[(6, 7)], known_tree_val)
        '''
        for k, (distances, predlist) in self.ntw.alldistances.items():
            self.assertEqual(distances[k], 0)
            for p, plists in predlist.items():
                self.assertEqual(plists[-1], k)
                self.assertEqual(self.ntw.node_list, list(predlist.keys()))
        '''
        matrix2 = self.ntw.allneighbordistances(self.pp1_str, fill_diagonal=0.)
        observed = matrix2.diagonal()
        known = np.zeros(matrix2.shape[0])
        self.assertEqual(observed.all(), known.all())
        
        matrix3 = self.ntw.allneighbordistances(self.pp1_str, snap_dist=True)
        known_mtx_val = 3218.2597894
        observed_mtx_val = matrix3
        self.assertAlmostEqual(observed_mtx_val[0, 1], known_mtx_val, places=4)
        
        matrix4 = self.ntw.allneighbordistances(self.pp1_str,
                                                fill_diagonal=0.)
        observed = matrix4.diagonal()
        known = np.zeros(matrix4.shape[0])
        self.assertEqual(observed.all(), known.all())
    
    def test_nearest_neighbor_distances(self):
        # general test
        with self.assertRaises(KeyError):
            self.ntw.nearestneighbordistances('i_should_not_exist')
        nnd1 = self.ntw.nearestneighbordistances(self.pp1_str)
        nnd2 = self.ntw.nearestneighbordistances(self.pp1_str,
                                                 destpattern=self.pp1_str)
        nndv1 = np.array(list(nnd1.values()))[:,1].astype(float)
        nndv2 = np.array(list(nnd2.values()))[:,1].astype(float)
        np.testing.assert_array_almost_equal_nulp(nndv1, nndv2)
        
        # nearest neighbor keeping zero test
        known_zero = ([19], 0.0)[0]
        nn_c = self.ntw.nearestneighbordistances(self.pp2_str,
                                                 keep_zero_dist=True)
        self.assertEqual(nn_c[18][0], known_zero)
        
        # nearest neighbor omitting zero test
        known_nonzero = ([11], 165.33982412719126)[1]
        nn_c = self.ntw.nearestneighbordistances(self.pp2_str,
                                                 keep_zero_dist=False)
        self.assertAlmostEqual(nn_c[18][1], known_nonzero, places=4)
        
        # nearest neighbor with snap distance
        known_neigh = ([3], 402.5219673922477)[1]
        nn_c = self.ntw.nearestneighbordistances(self.pp2_str,
                                                 keep_zero_dist=True,
                                                 snap_dist=True)
        self.assertAlmostEqual(nn_c[0][1], known_neigh, places=4)
    
    def test_element_as_gdf(self):
        obs = spgh.element_as_gdf(self.ntw, pp_name=self.pp1_str)
        snap_obs = spgh.element_as_gdf(self.ntw, pp_name=self.pp1_str,
                                       snapped=True)
        
        known_dist = 205.65961300587043
        observed_point = obs.loc[(obs['id']==0), 'geometry'].squeeze()
        snap_point = snap_obs.loc[(snap_obs['id']==0), 'geometry'].squeeze()
        observed_dist = observed_point.distance(snap_point)
        self.assertAlmostEqual(observed_dist, known_dist, places=8)
        
        with self.assertRaises(KeyError):
            spgh.element_as_gdf(self.ntw, pp_name='i_should_not_exist')


@unittest.skipIf(GEOPANDAS_EXTINCT, 'Missing Geopandas')
class TestNetworkAnalysis(unittest.TestCase):
    
    def setUp(self):
        path_to_shp = examples.get_path('streets.shp')
        gdf = geopandas.read_file(path_to_shp)
        self.ntw = spgh.Network(in_data=gdf)
        self.pt_str = 'schools'
        path_to_shp = examples.get_path('%s.shp' % self.pt_str )
        in_data = geopandas.read_file(path_to_shp)
        self.ntw.snapobservations(in_data, self.pt_str , attribute=True)
        npts = self.ntw.pointpatterns[self.pt_str].npoints
        self.ntw.simulate_observations(npts)
        self.test_permutations = 3
        self.test_steps = 5
    
    def tearDown(self):
        pass
    
    def test_network_f(self):
        obtained = self.ntw.NetworkF(self.ntw.pointpatterns[self.pt_str],
                                     permutations= self.test_permutations,
                                     nsteps=self.test_steps)
        self.assertEqual(obtained.lowerenvelope.shape[0], self.test_steps)
    
    def test_network_g(self):
        obtained = self.ntw.NetworkG(self.ntw.pointpatterns[self.pt_str],
                                     permutations= self.test_permutations,
                                     nsteps=self.test_steps)
        self.assertEqual(obtained.lowerenvelope.shape[0], self.test_steps)
    
    def test_network_k(self):
        obtained = self.ntw.NetworkK(self.ntw.pointpatterns[self.pt_str],
                                     permutations= self.test_permutations,
                                     nsteps=self.test_steps)
        self.assertEqual(obtained.lowerenvelope.shape[0], self.test_steps)

'''
@unittest.skipIf(GEOPANDAS_EXTINCT, 'Missing Geopandas')
class TestNetworkUtils(unittest.TestCase):
    
    def setUp(self):
        path_to_shp = examples.get_path('streets.shp')
        gdf = geopandas.read_file(path_to_shp)
        self.ntw = spgh.Network(in_data=gdf)
    
    def tearDown(self):
        pass
    
    def test_compute_length(self):
        self.point1, self.point2 = (0,0), (1,1)
        self.length = util.compute_length( self.point1, self.point2)
        self.assertAlmostEqual(self.length, 1.4142135623730951, places=4)
    
    def test_get_neighbor_distances(self):
        self.known_neighs = {1: 102.62353453439829, 2: 660.000001049743}
        self.neighs = util.get_neighbor_distances(self.ntw, 0,
                                                  self.ntw.edge_lengths)
        self.assertAlmostEqual(self.neighs[1], 102.62353453439829, places=4)
        self.assertAlmostEqual(self.neighs[2], 660.000001049743, places=4)
    
    def test_generate_tree(self):
        self.known_path = [23, 22, 20, 19, 170, 2, 0]
        self.distance, self.pred = util.dijkstra(self.ntw,
                                                 self.ntw.edge_lengths, 0)
        self.tree = util.generatetree(self.pred)
        self.assertEqual(self.tree[3], self.known_path)
    
    def test_dijkstra(self):
        self.distance, self.pred = util.dijkstra(self.ntw,
                                                 self.ntw.edge_lengths, 0)
        self.assertAlmostEqual(self.distance[196], 5505.668247, places=4)
        self.assertEqual(self.pred[196], 133)
    
    def test_dijkstra_mp(self):
        self.distance, self.pred = util.dijkstra_mp((self.ntw,
                                                     self.ntw.edge_lengths, 0))
        self.assertAlmostEqual(self.distance[196], 5505.668247, places=4)
        self.assertEqual(self.pred[196], 133)
    
    def test_square_distance_point_segment(self):
        self.point, self.segment = (1,1), ((0,0), (2,0))
        self.sqrd_nearp = util.squared_distance_point_segment(self.point,
                                                              self.segment)
        self.assertEqual(self.sqrd_nearp[0], 1.0)
        self.assertEqual(self.sqrd_nearp[1].all(), np.array([1., 0.]).all())
    
    def test_snap_points_on_segments(self):
        self.points = {0: cg.shapes.Point((1,1))}
        self.segments = [cg.shapes.Chain([cg.shapes.Point((0,0)),
                                          cg.shapes.Point((2,0))])]
        self.snapped = util.snap_points_on_segments(self.points, self.segments)
        self.known_coords = [xy._Point__loc for xy in self.snapped[0][0]]
        self.assertEqual(self.known_coords, [(0.0, 0.0), (2.0, 0.0)])
        self.assertEqual(self.snapped[0][1].all(), np.array([1., 0.]).all())
'''

if __name__ == '__main__':
    unittest.main()
