""" Testing for the spaghetti api import structure.
"""

import unittest

# from .network_unittest_classes import TestNetwork
# from .network_unittest_classes import TestNetworkPointPattern
from .network_unittest_classes import TestNetworkAnalysis

# api import structure
import spaghetti

# run tests on spaghetti.network.Network
# TestNetwork.spaghetti = spaghetti
# TestNetwork()

# run tests on spaghetti.network.PointPattern
# TestNetworkPointPattern.spaghetti = spaghetti
# TestNetworkPointPattern()

# run tests on spaghetti.analysis
TestNetworkAnalysis.spaghetti = spaghetti
TestNetworkAnalysis()

if __name__ == "__main__":
    unittest.main()
