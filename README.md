<p align="center">
<img src="figs/pysal_logo.png" width="200" height="200" />
</p>

[pysal/spaghetti](http://pysal.org/spaghetti/)
=================================

*SPA*tial *G*rap*H*s: n*ET*works, *T*opology, & *I*nference
============================================

Spaghetti is an open-source Python library for the analysis of network-based spatial data. Originating from the `network` module in [PySAL (Python Spatial Analysis Library)](http://pysal.org), it is under active development for the inclusion of newly proposed methods for building graph-theoretic networks and the analysis of network events. This package is part of a [refactoring of PySAL](https://github.com/pysal/pysal/wiki/PEP-13:-Refactor-PySAL-Using-Submodules).

*An example of shortest path plotting, routes originating from observation 0:*
<p align="center">
<img src="figs/shortest_path_plot.png" width="360" height="360" />
</p>


|[![PyPI version](https://badge.fury.io/py/spaghetti.svg)](https://badge.fury.io/py/spaghetti)| [![Conda Version](https://img.shields.io/conda/vn/conda-forge/spaghetti.svg)](https://anaconda.org/conda-forge/spaghetti) | ![tag](https://img.shields.io/github/v/release/pysal/spaghetti?include_prereleases&sort=semver) | [![GitHub issues open](https://img.shields.io/github/issues/pysal/spaghetti.svg?maxAge=3600)](https://github.com/pysal/spaghetti/issues) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pysal/spaghetti/master)
|:---:|:---:|:---:|:---:|:---:|
|[![Downloads](https://pepy.tech/badge/spaghetti)](https://pepy.tech/project/spaghetti) | [![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/spaghetti.svg)](https://anaconda.org/conda-forge/spaghetti) | [![Documentation](https://img.shields.io/static/v1.svg?label=docs&message=current&color=9cf)](http://pysal.org/spaghetti/) | [![GitHub issues closed](https://img.shields.io/github/issues-closed/pysal/spaghetti.svg?maxAge=3600)](https://github.com/pysal/spaghetti/issues) | [![Gitter](https://badges.gitter.im/pysal/Spaghetti.svg)](https://gitter.im/pysal/Spaghetti?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
| ![Pypi python versions](https://img.shields.io/pypi/pyversions/spaghetti.svg) | [![Conda Recipe](https://img.shields.io/badge/recipe-spaghetti-red.svg)](https://github.com/conda-forge/spaghetti-feedstock) | [![codecov](https://codecov.io/gh/pysal/spaghetti/branch/master/graph/badge.svg)](https://codecov.io/gh/pysal/spaghetti) | ![Github pull requests open](https://img.shields.io/github/issues-pr/pysal/spaghetti.svg) | [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
| [![image](https://travis-ci.org/pysal/spaghetti.svg)](https://travis-ci.org/pysal/spaghetti) | [![Build status](https://ci.appveyor.com/api/projects/status/eymi8wxdcmod95ge?svg=true)](https://ci.appveyor.com/project/pysal/spaghetti) | [![DOI](https://zenodo.org/badge/88305306.svg)](https://zenodo.org/badge/latestdoi/88305306) | ![Github pull requests closed](https://img.shields.io/github/issues-pr-closed/pysal/spaghetti.svg) | [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)


Examples
-----------
The following are a selection of some examples that can be launched individually as interactive binders from the links on their respective pages. Additional examples can be found in the [Tutorials](https://pysal.org/spaghetti/tutorials.html) section of the documentation. See the [`pysal/notebooks`](http://pysal.org/notebooks) project for a [`jupyter-book`](https://github.com/choldgraf/jupyter-book) version of this repository.
* [Network Representation](https://pysal.org/spaghetti/notebooks/quickstart.html)
* [Spatial Network Analysis](https://pysal.org/spaghetti/notebooks/network-analysis.html)
* [Optimal Facility Location](https://pysal.org/spaghetti/notebooks/facility-location.html)


Installation
------------

As of version 1.4.2, `spaghetti` officially supports Python [3.6](https://docs.python.org/3.6/), [3.7](https://docs.python.org/3.7/), and [3.8](https://docs.python.org/3.8/). Please make sure that you are operating in a Python >= 3.6 environment.

**Installing with `conda` via [`conda-forge`](https://github.com/conda-forge/spaghetti-feedstock) (highly recommended)**

To install `spaghetti` and all its dependencies, we recommend using the [`conda`](https://docs.conda.io/en/latest/)
manager, specifically with the [`conda-forge`](https://conda-forge.org) channel. This can be obtained by installing the [`Anaconda Distribution`](https://docs.continuum.io/anaconda/) (a free Python distribution for data science), or through [`miniconda`](https://docs.conda.io/en/latest/miniconda.html) (minimal distribution only containing Python and the conda package manager). 

Using `conda`, `spaghetti` can be installed as follows:
```
$ conda config --set channel_priority strict
$ conda install --channel conda-forge spaghetti
```

**Installing with [`PyPI`](https://pypi.org/project/spaghetti/)**
```
$ pip install spaghetti
```
*or* download the source distribution (`.tar.gz`) and decompress it to your selected destination. Open a command shell and navigate to the decompressed folder.
```
$ pip install .
```

***Warning***

When installing via `pip`, you have to ensure that the required dependencies for `spaghetti` are installed on your operating system. Details on how to install these packages are linked below. Using `conda` (above) avoids having to install the dependencies separately.

Install the most current development version of `spaghetti` by running:

```
$ pip install git+https://github.com/pysal/spaghetti
```


Requirements
----------------
- [`esda`](https://esda.readthedocs.io/en/latest/)
- [`libspatialindex`](https://libspatialindex.org/index.html)
- [`numpy`](https://numpy.org/devdocs/)
- [`rtree`](http://toblerity.org/rtree/install.html)
- [`scipy`](http://scipy.github.io/devdocs/)

Soft Dependencies
----------------------
- [`geopandas`](http://geopandas.org/install.html)
- [`shapely`](https://shapely.readthedocs.io/en/latest/)

Contribute
------------

PySAL-spaghetti is under active development and contributors are welcome.

If you have any suggestion, feature request, or bug report, please open a new [issue](https://github.com/pysal/spaghetti/issues) on GitHub. To submit patches, please review [PySAL: Getting Started](http://pysal.org/getting_started#for-developers), the PySAL [development guidelines](https://github.com/pysal/pysal/wiki), the `spaghetti` [contributing guidelines](https://github.com/pysal/spaghetti/blob/master/.github/CONTRIBUTING.md) before  opening a [pull request](https://github.com/pysal/spaghetti/pulls). Once your changes get merged, you’ll automatically be added to the [Contributors List](https://github.com/pysal/spaghetti/graphs/contributors).

Support
---------

If you are having issues, please [create an issue](https://github.com/pysal/spaghetti/issues) or talk to us in the [gitter room](https://gitter.im/pysal/spaghetti).


Code of Conduct
--------------------

As a PySAL-federated project, `spaghetti` follows the [Code of Conduct](https://github.com/pysal/governance/blob/master/conduct/code_of_conduct.rst) under the [PySAL governance model](https://github.com/pysal/governance).


License
---------

The project is licensed under the [BSD license](https://github.com/pysal/spaghetti/blob/master/LICENSE.txt).

BibTeX Citation
------------------

If you use PySAL-spaghetti in a scientific publication, we would appreciate using the following citation:

```
@misc{Gaboardi2018,
    author    = {Gaboardi, James D. and Laura, Jay and Rey, Sergio and Wolf, Levi John and Folch, David C. and Kang, Wei and Stephens, Philip and Schmidt, Charles},
    month     = {oct},
    year      = {2018},
    title     = {pysal/spaghetti},
    url       = {https://github.com/pysal/spaghetti},
    doi       = {10.5281/zenodo.1343650},
    keywords  = {graph-theory,network-analysis,python,spatial-networks,topology}
}
```
 
