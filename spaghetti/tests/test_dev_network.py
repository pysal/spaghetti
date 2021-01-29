""" Testing for the spaghetti development import structure.
"""

import unittest

from .network_unittest_classes import TestNetwork
from .network_unittest_classes import TestNetworkPointPattern
from .network_unittest_classes import TestNetworkAnalysis

from .network_unittest_classes import TestNetworkUtils

# dev import structure
from .. import network as spaghetti
from .. import util

# run tests on spaghetti.network.Network
TestNetwork.spaghetti = spaghetti
TestNetwork()

# run tests on spaghetti.network.PointPattern
TestNetworkPointPattern.spaghetti = spaghetti
TestNetworkPointPattern()

# run tests on spaghetti.analysis
TestNetworkAnalysis.spaghetti = spaghetti
TestNetworkAnalysis()

# run tests on spaghetti.util
TestNetworkUtils.spaghetti = spaghetti
TestNetworkUtils.util = util
TestNetworkUtils()

if __name__ == "__main__":
    unittest.main()
