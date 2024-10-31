# CryoXKit - Python API

### Requirements

Minimum requirements:
* Python (>=3.5)
* SWIG

### Installation (from source)

We recommend to install the Anaconda distribution (https://www.continuum.io/downloads) for a clean python environment with all prerequisites already installed. To install:
```bash
$ conda create -n cryoXkit python=3
$ conda activate cryoXkit
$ conda install -c conda-forge swig
```

To install the Python bindings:
```bash
$ python setup.py build install
```

