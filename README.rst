COSgen
=======

.. image:: https://readthedocs.org/projects/cosgen/badge/?version=latest
  :target: http://cosgen.readthedocs.io/en/latest/?badge=latest
  :alt: Documentation Status
.. image:: https://travis-ci.org/IBT-FMI/COSgen.svg?branch=master
  :target: https://travis-ci.org/IBT-FMI/COSgen

COSgen (Contrast optimized stimulus generator) is a Python package for the optimization of stimulus sequences in stimulus evoked fMRI experiments.
It is fully compatible with the stimulus sequence delivery solution COSplay_.

Features
--------

- Highly adaptable for specific use cases
- Allows full model (design matrix construction, convarinace matrix compution) specification
- Allows custom fitness measure specification
- Extensive API Documentation_
- Fully compatible with a pure Micro/Python stimulus delivery solution (COSplay_)

Installation
------------

COSgen can be installed using Python's `setuptools`.

.. code-block:: console

   git clone https://github.com/IBT-FMI/COSgen.git
   cd COSgen
   python setup.py install --user

*Note:* `setuptools` will not manage any dependencies.
You will have to install Dependencies_ manually, e.g. using ``pip install argh``.

Dependencies
------------

- Python_ 2.7 or Python 3.5 and newer
- Numpy_ 1.13 or newer
- Scipy_ 0.19 or newer
- Argh_ 0.26 or newer

.. _Python: https://www.python.org/
.. _COSplay: https://github.com/IBT-FMI/COSplay
.. _Documentation: http://cosgen.readthedocs.io/en/latest/
.. _Numpy: http://www.numpy.org/
.. _Scipy: https://www.scipy.org/
.. _Argh: https://pypi.python.org/pypi/argh
