# Change-log

All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## Version 0.1.1 / 2022-09-28

### Added

* Updated to PyEchelle v. 0.3.1
* Updated to Pyxel v. 1.3.0
* It is now possible to make multi-target simulations (PyEchelle)
* Dark frames has now been added to the calibration image products
* Fast and slow readout mode of the CCD has been added as an option
* Added a LICENSE and this CHANGELOG.md
* New example: A Pyxel validation Jupyter notebook
* New example: Scripts to simulate using High Performance Computing (HPC)

### Changed

* PyEchelle's new Python interface is used instead the previous bash interface
* New parse arguments to select the fast and slow CCD readout mode
* Improved fits header saved to each file (software versions, exposure time, readmode, etc.)

### Bugfix

* Froze PyEchelle v. 0.3.1 and numba v. 0.55.1 in ``pyproject.toml`` due to cluster specfic NVIDIA drives
* Bugfix of package ``zipfile` to correctly deflating files and not only store them

## Version 0.1.0 / 2021-06-01

### Added

* Initial release of MARVELsim
* Usage of PyEchelle v. 0.2.0
* Usage of Pyxel v. 0.11.7
* Documentation for MARVELsim
