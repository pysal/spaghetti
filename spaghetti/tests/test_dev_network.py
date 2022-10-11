""" Testing for the spaghetti development import structure.
"""

# dev import structure
from .. import network as spaghetti
from .. import util
from .network_test_classes import (
    TestNetwork,
    TestNetworkAnalysis,
    TestNetworkPointPattern,
    TestNetworkUtils,
)

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
