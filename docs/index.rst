.. documentation master file

.. raw:: html

    <img 
        src="_static/images/spaghetti_nav_logo.svg" 
        class="img-responsive center-block" 
        alt="spaghetti logo" 
        width="470" 
        height="200"
    >

`spaghetti`
===========

**SPA**\ tial **G**\ rap\ **H**\ s: n\ **ET**\ works, **T**\ opology, & **I**\ nference
---------------------------------------------------------------------------------------

`Spaghetti` is an open-source Python library for the analysis of network-based
spatial data. Originating from the `network` module in `PySAL (Python Spatial 
Analysis Library) <http://pysal.org>`_, it is under active development for the 
inclusion of newly proposed methods for building graph-theoretic networks and 
the analysis of network events. An installation guide, API reference, 
and usage tutorials are provided here through the links above.

.. raw:: html

    <div class="container-fluid">
      <div class="row equal-width">
        <div class="col-sm-.5 col-xs-hidden">
        </div>
        <div class="col-md-6 col-xs-15">
            <a 
                href="https://pysal.org/spaghetti/notebooks/network-segmentation.html" 
                class="thumbnail"
            >
                <img 
                    src="_static/images/crime_counts.png" 
                    class="img-responsive center-block"
                >
                <div class="caption text-center">
                <h6>Network Representations</h6>
                </div>
            </a>
        </div>
        <div class="col-sm-6 col-xs-15">
            <a 
                href="https://pysal.org/spaghetti/notebooks/network-spatial-dependence.html" 
                class="thumbnail"
            >
                <img 
                    src="_static/images/network_k.png" 
                    class="img-responsive center-block"
                >
                <div class="caption text-center">
                <h6>Network Spatial Dependence</h6>
                </div>
            </a>
        </div>
        </div>
        <div class="col-sm-.5 col-xs-hidden">
        </div>
      </div>
    </div>

History
-------

`Spaghetti` was 
created and has evolved in line with the Python Spatial Analysis Library ecosystem for 
the specific purpose of utilizing the functionality of spatial weights in 
`libpysal <https://pysal.org/libpysal/>`_ for generating network segment contiguity objects. 
The PySAL project was started in the mid-2000s when installation was difficult to maintain. 
Due to the non-triviality of relying on dependencies to secondary packages, a conscious 
decision was made to limit dependencies and build native PySAL data structures in cases 
where at all possible. Therefore, the original `pysal.network` submodule was developed to 
address the need for integrating support for network data structures with PySAL weights 
data structures, with the target audience being spatial data scientists and anyone 
interested in investigating network-centric phenomena within PySAL. Owing to the 
co-development of network functionality found within `spaghetti` and the evolution of 
the wider PySAL ecosystem, today, the package provides specialized network functionality 
that easily integrates with the rest of PySAL. This allows users of `spaghetti`’s network 
functionality to access spatial analysis functionality that complements network analysis, 
such as spatial statistical tools with `esda` and integration with core components of 
`libpysal`: `libpysal.weights` (mentioned above), 
`libpysal.cg` (computational geometry and data structures), 
`libpysal.io` (input-output), and `libpysal.examples` (built-in example data).

Development
-----------

Development of `spaghetti` is hosted on GitHub_.

Support
-------

All questions, comments, & discussions should happen in a public forum, where possible. 
Please start a `discussion <https://github.com/pysal/spaghetti/discussions>`_ for questions, talk to us in `PySAL's Discord channel <https://discord.gg/BxFTEPFFZn>`_, or open an `issue <https://github.com/pysal/spaghetti/issues>`_ if there appears to be a 
bug. Private messages and emails will not be answered in a substantive manner.

Citing `spaghetti`
------------------

If you use PySAL-spaghetti in a scientific publication, we would appreciate using the following BibTeX citations::

  @article{Gaboardi2021,
    doi         = {10.21105/joss.02826},
    url         = {https://doi.org/10.21105/joss.02826},
    year        = {2021},
    publisher   = {The Open Journal},
    volume      = {6},
    number      = {62},
    pages       = {2826},
    author      = {James D. Gaboardi and Sergio Rey and Stefanie Lumnitz},
    title       = {spaghetti: spatial network analysis in PySAL},
    journal     = {Journal of Open Source Software}
  }
  
  @misc{Gaboardi2018,
    author      = {Gaboardi, James D. and Laura, Jay and Rey, Sergio and
                   Wolf, Levi John and Folch, David C. and Kang, Wei and 
                   Stephens, Philip and Schmidt, Charles},
    month       = {oct},
    year        = {2018},
    title       = {pysal/spaghetti},
    url         = {https://github.com/pysal/spaghetti},
    doi         = {10.5281/zenodo.1343650},
    keywords    = {graph-theory,network-analysis,python,spatial-networks,topology}
  }

Citing Work
-----------

* **Lovelace, R**. `Open source tools for geographic analysis in transport planning`. Journal of Geographical Systems (2021): 1-32. https://doi.org/10.1007/s10109-020-00342-2.
* **Rey, Sergio J., et al**. `The PySAL Ecosystem: Philosophy and Implementation`. Geographical Analysis (2022): 467-487. https://doi.org/10.1111/gean.12276.
* **Barboza-Salerno, Gia E., and Jacquelyn CA Meshelemiah**. `Gun Violence on Walkable Routes to and from School: Recommendations for Policy and Practice`. Journal of Urban Health (2023): 1-16. https://doi.org/10.1007/s11524-023-00802-2
* **Barboza, Gia, and Jacquelyn Meshelemiah**. `Danger, Students Beware, School Ahead! Gun Violence Exposure Near Schools in Compton, California`. (2023). https://doi.org/10.21203/rs.3.rs-2976516/v1

Funding
-------

This project is/was partially funded through:

.. figure:: _static/images/ardc_logo.png
    :target: https://atlantardc.wordpress.com
    :width: 150
    :align: left

    The Atlanta Research Data Center: `A Polygon-Based Approach to Spatial Network Allocation <https://atlantardc.files.wordpress.com/2018/05/ardc-newsletter_2018_2.pdf>`_

.. figure:: _static/images/nsf_logo.png
    :target: https://www.nsf.gov/index.jsp
    :width: 100
    :align: left

    National Science Foundation Award #1825768: `National Historical Geographic Information System <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1825768&HistoricalAwards=false>`_

.. raw:: html

    <img 
        src="_static/images/pysal_logo.svg" 
        class="img-responsive center-block" 
        alt="PySAL Logo" 
        width="400" 
        height="400"
    >

.. toctree::
   :hidden:
   :maxdepth: 3
   :caption: Contents:

   Installation <installation>
   Tutorials <tutorials>
   API <api>
   References <references>

.. _PySAL: https://github.com/pysal/pysal
.. _GitHub: https://github.com/pysal/spaghetti
