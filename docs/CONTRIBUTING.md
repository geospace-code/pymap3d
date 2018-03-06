## Issues 

GitHub Issues are a great way to report problems or request features.
The goal of PyMap3D is to provide a Matlab Mapping Toolbox - like API for 3D coordinate conversion,
allowing for arbitrarily shaped input/output arrays where feasible.

## Pull Requests
It's always nice to get pull requests. Keeping in mind the above goals 
and that PyMap3D is intended for non-interactive use on massively parallel (HPC) and embedded systems.


### Python
Should use Numpy for arrays, falling back to `math` for scalar case when Numpy isn't installed.

### Matlab/Octave
Functions must work in latest GNU Octave release and Matlab.

### Fortran
Generally procedures should be Pure Elemental or at least Impure Elemental to allow
massively parallel HPC usage.
