// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include <benchmark/benchmark.h>

#include "FermionicBasis.h"
#include "Models/HubbardChain.h"
#include "SparseMatrix.h"

static auto basis_range = benchmark::CreateDenseRange(8, 12, 2);

static void BM_CreateHubbardChainMatrixElements(benchmark::State& state) {
  for (auto _ : state) {
    state.PauseTiming();
    const std::size_t size = state.range(0);
    const std::size_t particles = state.range(1);
    const double mu = 0.0;
    const double t = 1.0;
    const double u = 2.0;

    HubbardChain model(mu, t, u, size);
    FermionicBasis basis(size, particles);
    SparseMatrix<std::complex<double>> m;
    state.ResumeTiming();
    model.compute_matrix_elements(basis, m);
  }
}

BENCHMARK(BM_CreateHubbardChainMatrixElements)
    ->ArgsProduct({basis_range, basis_range});
