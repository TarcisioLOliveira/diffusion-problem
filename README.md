# Diffusion problem analysis

Some tests done for a paper, as an analogy for linear elasticity in density-based
topology optimization. See also [min-diag](https://github.com/TarcisioLOliveira/min-diag).

## Dependencies
- [BLAS/LAPACK](https://www.netlib.org/)
- [Eigen3](https://eigen.tuxfamily.org/)
- [SFML](https://www.sfml-dev.org/)

## Implemented tests
- `test1`: Minimal material analogy using only Dirichlet boundary conditions.
- `test2`: Minimal material analogy using Dirichlet and Neumann boundary conditions.
- `test3`: Minimal diagonal approach using only Dirichlet boundary conditions.
- `test4`: Minimal diagonal approach using Dirichlet and Neumann boundary conditions.
