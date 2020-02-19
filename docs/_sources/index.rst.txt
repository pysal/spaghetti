.. documentation master file

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
        <div class="col-md-4 col-xs-15">
            <a href="https://pysal.org/spaghetti/notebooks/quickstart.html" class="thumbnail">
                <img src="_static/images/net_rep.png" class="img-responsive center-block">
                <div class="caption text-center">
                <h6>Network Representation</h6>
                </div>
            </a>
        </div>
        <div class="col-sm-4 col-xs-15">
            <a href="https://pysal.org/spaghetti/notebooks/network-analysis.html" class="thumbnail">
                <img src="_static/images/network_k.png" class="img-responsive center-block">
                <div class="caption text-center">
                <h6>Spatial Network Analysis</h6>
                </div>
            </a>
        </div>
        <div class="col-sm-4 col-xs-15">
            <a href="https://pysal.org/spaghetti/notebooks/facility-location.html" class="thumbnail">
                <img src="_static/images/facility_location.png"
                class="img-responsive center-block">
                <div class="caption text-center">
                <h6>Optimal Facility Location</h6>
                </div>
            </a>
        </div>
        <div class="col-sm-.5 col-xs-hidden">
        </div>
      </div>
    </div>


Development
-----------

Development of `spaghetti` is hosted on GitHub_.


Citing `spaghetti`
------------------

If you use PySAL-spaghetti in a scientific publication, we would appreciate using the following citation:

  Bibtex entry::

      @misc{Gaboardi2018,
        author   = {Gaboardi, James D. and Laura, Jay and Rey, Sergio and
                    Wolf, Levi John and Folch, David C. and Kang, Wei and 
                    Stephens, Philip and Schmidt, Charles},
        month    = {oct},
        year     = {2018},
        title    = {pysal/spaghetti},
        url      = {https://github.com/pysal/spaghetti},
        doi      = {10.5281/zenodo.1343650},
        keywords = {graph-theory,network-analysis,python,spatial-networks,topology}
      }


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
