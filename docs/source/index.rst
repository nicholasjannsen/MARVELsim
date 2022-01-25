.. MARVELsim documentation master file, created by
   sphinx-quickstart on Sat Jan 22 20:31:34 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MARVELsim's documentation!
=====================================

.. image:: marvel_setup.png
   :align: center
   :width: 700

This documentation is an easy guide how to sucessfully generate realistic MARVEL data products. The simulations are build from two popular main modules named `PyEchelle <https://stuermer.gitlab.io/pyechelle/index.html>`_ and `Pyxel <https://esa.gitlab.io/pyxel/>`_. `PyEchelle <https://stuermer.gitlab.io/pyechelle/index.html>`_ is a Python based fast generic spectrum simulator that is used to generate each object or calibration spectrum, whereas `Pyxel <https://esa.gitlab.io/pyxel/>`_ is a general imaging detector simulation framework in Python that is used to efficiently add-on all necessary CCDs effects.
	 
More information
----------------

Full list of PyEchelle resources:
  - `PyEchelle's documentation <https://stuermer.gitlab.io/pyechelle/index.html>`_
  - `PyEchelle's GitLab repository <https://gitlab.com/Stuermer/pyechelle>`_
  - `PyEchelle's precusor's publication (Echelle++) <https://iopscience.iop.org/article/10.1088/1538-3873/aaec2e/pdf>`_

Full list of Pyxel resources:
  - `Pyxel's Documentation <https://esa.gitlab.io/pyxel/>`_
  - `Pyxel's Simulator GitLab repository <https://gitlab.com/esa/pyxel>`_
  - `Pyxel's Data-example GitLab repository <https://gitlab.com/esa/pyxel-data>`_
  - `Pyxel's Gitter helpdesk community <https://gitter.im/pyxel-framework/community>`_

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   overview
   installation
   tutorial
   troubleshooting
   acknowledgements
