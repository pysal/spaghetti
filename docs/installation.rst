.. Installation

Python 3.10_ is tested for support by `spaghetti`. Please make sure that you are operating in a Python >= 3.10 environment.

Installation
============

Installing with ``conda`` via `spaghetti-feedstock`_ (highly recommended)
-------------------------------------------------------------------------

To install `spaghetti` and all its dependencies, we recommend using the conda_ manager, specifically with the conda-forge_ channel. This can be obtained by installing the `Anaconda Distribution`_ (a free Python distribution for data science), or through miniconda_ (minimal distribution only containing Python and the ``conda`` package manager). 

Using ``conda``, `spaghetti` can be installed as follows::

  $ conda config --set channel_priority strict
  $ conda install --channel conda-forge spaghetti

Also, ``geopandas`` provides `a nice example`_ to create a fresh environment for working with spatial data.

Installing with `Python Package Index`_
---------------------------------------
::

  $ pip install spaghetti

*or* download the source distribution (``.tar.gz``) and decompress it to your selected destination. Open a command shell and navigate to the decompressed folder. ::

  $ pip install .

.. role:: rubric

**Warning**

When installing via ``pip``, you have to ensure that the required dependencies for `spaghetti` are installed on your operating system. Details on how to install these packages are linked here_. Using ``conda`` (above) avoids having to install the dependencies separately.

Development Version
-------------------

Install the most current development version of `spaghetti` by running::

  $ pip install git+https://github.com/pysal/spaghetti

You can  also fork_ the `pysal/spaghetti`_ repo and create a local clone of your fork. By making changes to your local clone and submitting a pull request to `pysal/spaghetti`_, you can contribute to the spaghetti development.

|

.. _3.10: https://docs.python.org/3.10/
.. _spaghetti-feedstock: https://github.com/conda-forge/spaghetti-feedstock
.. _a nice example: https://geopandas.readthedocs.io/en/latest/getting_started/install.html#creating-a-new-environment
.. _conda: https://docs.conda.io/en/latest/
.. _conda-forge: https://conda-forge.org
.. _Anaconda Distribution: https://docs.continuum.io/anaconda/
.. _miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _Python Package Index: https://pypi.org/project/spaghetti/
.. _pysal/spaghetti: https://github.com/pysal/spaghetti
.. _fork: https://help.github.com/articles/fork-a-repo/
.. _here: https://github.com/pysal/spaghetti#requirements
