## Issues

GitHub Issues are a great way to report problems or request features.
This project adheres to [Semantic Versioning](https://semver.org).
The goal of PyMap3D is to provide a Matlab Mapping Toolbox - like API for 3D coordinate conversion,
allowing for arbitrarily shaped input/output arrays where feasible.

## Pull Requests

It's always nice to get pull requests.
PyMap3D is intended for non-interactive use from massively parallel (HPC) down to embedded systems.

Most PyMap3D functions require Numpy for arrays.
Some of the most essential functions fall back to Python stdlib `math` for scalar cases when Numpy isn't available.
