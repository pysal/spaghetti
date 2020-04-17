import numpy


class FuncBase(object):
    """Base object for performing network analysis on a
    ``spaghetti.Network`` object.
    
    Parameters
    ----------
    
    ntw : spaghetti.Network
        A spaghetti network object.
    
    pointpattern : spaghetti.network.PointPattern
        A spaghetti point pattern object.
    
    nsteps : int
        The number of steps at which the count of the nearest
        neighbors is computed.
        
    permutations : int
        The number of permutations to perform. Default 99.
    
    threshold : float
        The level at which significance is computed.
        (0.5 would be 97.5% and 2.5%).
    
    distribution : str
        The distribution from which random points are sampled
        Either ``"uniform"`` or ``"poisson"``.
    
    upperbound : float
        The upper bound at which the `K`-function is computed.
        Defaults to the maximum observed nearest neighbor distance.
    
    Attributes
    ----------
    
    sim : numpy.ndarray
        A simulated distance matrix.
    
    npts : int
        The number of points (``pointpattern.npoints``).
    
    xaxis : numpy.ndarray
        The observed x-axis of values.
    
    observed : numpy.ndarray
        The observed y-axis of values.
    
    """

    def __init__(
        self,
        ntw,
        pointpattern,
        nsteps=10,
        permutations=99,
        threshold=0.5,
        distribution="poisson",
        upperbound=None,
    ):

        # set initial class attributes
        self.ntw = ntw
        self.pointpattern = pointpattern
        self.nsteps = nsteps
        self.permutations = permutations
        self.threshold = threshold

        # set and validate the distribution
        self.distribution = distribution
        self.validatedistribution()

        # create an empty array to store the simulated points
        self.sim = numpy.empty((permutations, nsteps))
        self.npts = self.pointpattern.npoints

        # set the upper bounds
        self.upperbound = upperbound

        # compute the statistic K
        self.computeobserved()
        self.computepermutations()

        # compute the envelope vectors
        self.computeenvelope()

    def validatedistribution(self):
        """enusure statistical distribution is supported
        """

        valid_distributions = ["uniform", "poisson"]
        assert self.distribution in valid_distributions, (
            "Distribution not in %s" % valid_distributions
        )

    def computeenvelope(self):
        """compute upper and lower bounds of envelope
        """

        upper = 1.0 - self.threshold / 2.0
        lower = self.threshold / 2.0

        self.upperenvelope = numpy.nanmax(self.sim, axis=0) * upper
        self.lowerenvelope = numpy.nanmin(self.sim, axis=0) * lower

    def setbounds(self, nearest):
        """set upper and lower bounds
        """
        if self.upperbound is None:
            self.upperbound = numpy.nanmax(nearest)


class KFunc(FuncBase):
    """Compute a network constrained `K` statistic. This requires the
    capability to compute a distance matrix between two point patterns.
    In this case one will be observed and one will be simulated.
    
    
    Compute a network constrained `K` statistic. This requires the
    capability to compute a distance matrix between two point patterns.
    In this case one will be observed and one will be simulated....................................
    
    
    
    Attributes
    ----------
    
    lam : float
        The ``lambda`` value.
    
    Notes
    -----
    
    Based on :cite:`Okabe2001`.
    
    """

    def computeobserved(self):
        """Compute the observed nearest.
        """

        # pairwise distances
        distances = self.ntw.allneighbordistances(self.pointpattern)
        self.setbounds(distances)

        # set the intensity (lambda)
        self.lam = self.npts / sum(self.ntw.arc_lengths.values())

        # compute a K-Function
        observedx, observedy = kfunction(
            self.npts, distances, self.upperbound, self.lam, nsteps=self.nsteps
        )

        # set observed values
        self.observed = observedy
        self.xaxis = observedx

    def computepermutations(self):
        """Compute permutations of the points.
        """

        # for each round of permutations
        for p in range(self.permutations):

            # simulate a point pattern
            sim = self.ntw.simulate_observations(
                self.npts, distribution=self.distribution
            )

            # distances
            distances = self.ntw.allneighbordistances(sim)

            # compute a K-Function
            simx, simy = kfunction(
                self.npts, distances, self.upperbound, self.lam, nsteps=self.nsteps
            )

            # label the permutation
            self.sim[p] = simy


def kfunction(n_obs, dists, upperbound, intensity, nsteps=10):
    """Compute a `K`-function.

    Parameters
    ----------
    
    n_obs : int
        The number of observations. See ``self.npts``.
    
    dists : numpy.ndarray
        A matrix of pairwise distances.
    
    upperbound : int or float
        The end value of the sequence.
    
    intensity : float
        lambda value
    
    nsteps : int
        The number of distance bands. Default is 10. Must be
        non-negative.
    
    Returns
    -------
    
    x : numpy.ndarray
        The x-axis of values.
    
    y : numpy.ndarray
        The y-axis of values.
    
    """

    # create interval for x-axis
    x = numpy.linspace(0, upperbound, nsteps)

    # create empty y-axis vector
    y = numpy.empty(x.shape[0])

    # iterate over x-axis interval
    for i, r in enumerate(x):

        # slice out and count neighbors within radius
        with numpy.errstate(invalid="ignore"):
            y[i] = dists[dists <= r].shape[0]

    # compute k for y-axis vector
    y /= n_obs * intensity

    return x, y
