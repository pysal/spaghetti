""" Testing for the spaghetti api import structure.
"""

from .network_test_classes import TestNetwork
from .network_test_classes import TestNetworkPointPattern
from .network_test_classes import TestNetworkAnalysis

# api import structure
import spaghetti

# run tests on spaghetti.network.Network
TestNetwork.spaghetti = spaghetti
TestNetwork()

# run tests on spaghetti.network.PointPattern
TestNetworkPointPattern.spaghetti = spaghetti
TestNetworkPointPattern()

# run tests on spaghetti.analysis
TestNetworkAnalysis.spaghetti = spaghetti
TestNetworkAnalysis()
