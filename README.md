# ![polystar]
A C++ library for polygon and polyhedron operations.
Wrapped for use in python using [pybind11](https://github.com/pybind/pybind11).

[polystar]: https://raw.githubusercontent.com/g5t/polystar/master/polystar.svg

# Dependencies
## TetGen
A modified version of [TetGen](http://tetgen.org) is used to create
refined tetrahedral meshes in bounding polyhedra.

The modified version is included as part of this repository.

# Installation
Install via `pip`, e.g.,
```cmd
python -m pip install polystar
```

This repository can be installed locally via
```cmd
python -m pip install .
# or 
python setup.py install
```

Alternatively, the python module, C++ library, and [catch2](https://github.com/catchorg/Catch2) based tests can be built directly using `cmake`.
