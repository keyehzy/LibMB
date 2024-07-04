![Build Status](https://github.com/keyehzy/LibMB/actions/workflows/cmake.yml/badge.svg)
# LibMB

This project provides a C++ library for second-quantization calculations, enabling the representation and manipulation of quantum many-body systems using creation and annihilation operators. 

## Features

- **Symbolic representation of second-quantized operators:** Define Hamiltonians and other operators using a flexible and intuitive syntax.
- **Basis generation:**  Construct custom basis sets tailored to your problem, including restrictions on particle number, spin, and other quantum numbers.
- **Matrix representation:**  Efficiently compute matrix elements of operators in the chosen basis, facilitating numerical diagonalization.
- **Integration with external libraries:** Seamlessly interface with linear algebra libraries like Armadillo for powerful numerical computations.

## Usage

```bash
# Build the project using cmake
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build

# Run the tests
./build/tests/libwick-test
```

## Examples

In the following example, we construct a Hubbard chain model and compute the ground state using the Armadillo library.

```cpp
// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include <armadillo>  //  for eigensolver
#include <iostream>

#include "Assert.h"
#include "BasisFilter.h"
#include "FermionicBasis.h"
#include "Model.h"

using enum Operator::Type;        // for Creation, Annihilation
using enum Operator::Statistics;  // for Fermion
using enum Operator::Spin;        // Up, Down

class HubbardChain : public Model {
 public:
  HubbardChain(double mu, double t, double u, size_t n)
      : m_mu(mu), m_t(t), m_u(u), m_size(n) {}

  ~HubbardChain() override {}

 private:
  Expression hamiltonian() const override {
    Expression result;

    // Chemical potential and Hopping
    for (Operator::Spin spin : {Up, Down}) {
      for (std::size_t i = 0; i < m_size; i++) {
        result += density<Fermion>(-m_mu, spin, i);
        result += hopping<Fermion>(-m_t, spin, i, (i + 1) % m_size);
      }
    }

    // Hubbard U
    for (size_t i = 0; i < m_size; i++) {
      result += density_density<Fermion>(m_u, Up, i, Down, i);
    }

    return result;
  }

  double m_mu;
  double m_t;
  double m_u;
  size_t m_size;
};

int main() {
  const std::size_t size = 2;
  const std::size_t particles = 2;
  const double mu = 0.0;
  const double t = 1.0;
  const double u = 2.0;

  // Create the model
  HubbardChain model(mu, t, u, size);

  // Construct a basis
  FermionicBasis basis(size, particles);

  // Compute matrix elements
  arma::SpMat<arma::cx_double> m(basis.size(), basis.size());
  model.compute_matrix_elements(basis, m);
  LIBMB_ASSERT(m.is_hermitian());

  // Compute ground state using, e.g. Armadillo library
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  const std::size_t eigval_count = 1;
  arma::eigs_gen(eigval, eigvec, m, eigval_count, "sr");

  double gs_energy = eigval(0).real();
  const arma::cx_vec& ground_state = eigvec.col(0);

  // Perform some further analysis here...
}
```

1. The code defines a `HubbardChain` class, inheriting from a `Model` class, to
represent the Hubbard Hamiltonian. The hopping and interaction terms are
constructed using second-quantized operators.

2. A `Basis` object is created, specifying the Hilbert space for the calculation.
This example restricts the basis to states with a fixed particle number and
total spin projection.

3. The Hamiltonian's matrix representation in the chosen basis is computed. The
Armadillo library's `eigs_sym` function efficiently finds the lowest eigenvalues
and eigenvectors.

This example highlights LibMB's core functionalities, demonstrating its
flexibility in tackling quantum many-body problems. You can easily adapt this
code to explore different models, parameter regimes, or analysis techniques by
modifying the model definition, basis construction, or post-processing steps.

## License

This project is licensed under the BSD 2-Clause License - see the [LICENSE](LICENSE) file for details.
