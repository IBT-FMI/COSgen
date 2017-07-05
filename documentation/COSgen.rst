======
COSgen
======

Overview
========
Contrast optimized stimulus generator (COSgen) is a highly parameterized genetic algorithm implementation to produce optimized stimulus sequences for clinical and preclinical fMRI.
Our package is fully compatible with a pure Micro/Python stimulus delivery solution (COSplay), and provides a convenient and well documented API. Because of its modular structure, the implementation is highly adaptable for specific use cases.


Features
========

- High adaptability for specific use cases
- Full compatability with COSgen
- Custom model specification (design matrix construction, covariance matrix computation)
- Custom fitness measure specification

Dependencies
============

- Python 2.7 or Python 3.5 and newer
- Numpy
- Scipy
- argh
- nibabel (only for `make test data` and `estimate autocorr`)
- matplotlib (only for `cosgen.models.plot_design_matrix` and the example :ref:`block design example`)

Installation
============
The module can be installed using setuptools. 
Run ``python setup.py install`` inside the COSgen folder you downloaded. 
After installation you can use ``COSgen``, ``make_test_data`` and ``estimate_autocorr`` directly as commands.


References
==========
* COSplay: https://readthedocs.org/projects/cosplay/
* argh: https://pypi.python.org/pypi/argh
* Python: http://www.python.org/
* Setuptools: https://setuptools.readthedocs.io/en/latest/

