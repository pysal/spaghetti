""" Testing for the spaghetti api import structure.
"""

# api import structure
import spaghetti

from .network_test_classes import (
    TestNetwork,
    TestNetworkAnalysis,
    TestNetworkPointPattern,
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
