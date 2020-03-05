# Stochastic SVD

For further reading on this method, see [*Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions*].(https://arxiv.org/pdf/0909.4061.pdf)

The stochastic SVD is applied to square matrices of arbitrary dimension. To evaluate its efficacy, one may specify the elements of a diagonal matrix and apply arbitrary rotations to smear information away from the diagonal. A successful factorization ought to return the original eigenvalues. This work differs somewhat from that of Halko et al. in that this program includes an algorithm that attempts to determine the minimum number of components that can accurately represent the original matrix.

## Altering the C++ file
Currently, script creates a 7-dimensional matrix whose upper right 4x4 subspace contains relevant physics. The 0th and 3rd eigenvalues are -50 and +50 while the middle two are +1. This can be changed by modifying `N`, `n`, and `E` within `SVD_ND_Test.cpp`. These variables correspond to the full matrix dimension, the diagonal elements to set as 1, and the nonzero, nonunity eigenvalues respectively. Modifying any other portion of the code without proficiency in c++ is not recommended.

## Compiling
The clang compiler does not play nicely with this script. The GNU compiler is recommended (the executable in this repo was produced with g++-mp-6). Naturally, one must have the GNU Scientific Library installed and pass `-I /path/to/gsl/*.h` as well as the `-L /path/to/libgsl` and `-lgsl` linker flags.