/**
 * @file highfm_ode_solver.cpp
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief A class to solve ODEs using the trapezoidal method through HighFM and PARDISO
 * @version 0.1
 * @date 2025-07-08
 *
 * @copyright Copyright (c) 2025
 *
 */

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE long

#include "highfm_ode_solver.h"

#include <fmt/core.h>

#include <Eigen/Sparse>
#include <highfm/linsys/amgcl.hpp>

HigFMODESolver::HigFMODESolver(const Eigen::SparseMatrix<float>& A,
                               const std::vector<std::unique_ptr<TimeFunction>>& b,
                               const std::vector<float>& x_0, float t, size_t row, size_t N)
    : _A(A), _b(b), _x_0(x_0), _t(t), _row(row), _N(N) {}

std::vector<float> HigFMODESolver::solve_sequence() const {
    size_t dim = _x_0.size();
    float dt = _t / static_cast<float>(_N);

    /* Copy x_0 into mutable x */
    HighFM::Vector<float> x(dim);
    std::memcpy(x.data(), _x_0.data(), dim * sizeof(float));

    // Initialize result vector with
    std::vector<float> result;
    result.reserve(_N + 1);
    result.push_back(x(_row));

    Eigen::SparseMatrix<float, Eigen::RowMajor, HighFM::index_t> A_csr = _A;

    HighFM::CSRMap<float> A(HighFM::CSRMatrixData<float>{.nrows = static_cast<HighFM::index_t>(dim),
                                                         .ncols = static_cast<HighFM::index_t>(dim),
                                                         .nnz = A_csr.nonZeros(),
                                                         .row_offset = A_csr.outerIndexPtr(),
                                                         .nz_columns = A_csr.innerIndexPtr(),
                                                         .nz_values = A_csr.valuePtr()});

    HighFM::CSRMatrix<float> LHS(dim, dim), RHS_mat(dim, dim);

    LHS = HighFM::Eye<float>() - 0.5f * dt * A;
    RHS_mat = HighFM::Eye<float>() + 0.5f * dt * A;

    HighFM::Pardiso<float> solver;
    solver.factorize(LHS);

    float t = 0.0f;
    HighFM::Vector<float> b_n(dim), b_np1(dim), rhs(dim);
    for (size_t j = 0; j < dim; ++j) b_n(j) = _b[j]->operator()(t);

    for (size_t i = 0; i < _N; ++i) {
        t += dt;

        // Evaluate b(t_{n+1})
        for (size_t j = 0; j < dim; ++j) b_np1(j) = _b[j]->operator()(t);

        // Compute RHS
        rhs = RHS_mat * x;
        rhs = rhs + 0.5f * dt * (b_n + b_np1);

        // Solve system
        x = solver.solve(rhs);

        result.push_back(x(_row));
        b_n.swap(b_np1);  // Avoid copy, reuse buffer
    }

    return result;
}

std::vector<float> HigFMODESolver::solve_sequence_cg() const {
    size_t dim = _x_0.size();
    float dt = _t / static_cast<float>(_N);

    /* Copy x_0 into mutable x */
    HighFM::Vector<float> x(dim);
    std::memcpy(x.data(), _x_0.data(), dim * sizeof(float));

    // Initialize result vector with
    std::vector<float> result;
    result.reserve(_N + 1);
    result.push_back(x(_row));

    Eigen::SparseMatrix<float, Eigen::RowMajor, HighFM::index_t> A_csr = _A;

    HighFM::CSRMap<float> A(HighFM::CSRMatrixData<float>{.nrows = static_cast<HighFM::index_t>(dim),
                                                         .ncols = static_cast<HighFM::index_t>(dim),
                                                         .nnz = A_csr.nonZeros(),
                                                         .row_offset = A_csr.outerIndexPtr(),
                                                         .nz_columns = A_csr.innerIndexPtr(),
                                                         .nz_values = A_csr.valuePtr()});

    HighFM::CSRMatrix<float> LHS(dim, dim), RHS_mat(dim, dim);

    LHS = HighFM::Eye<float>() - 0.5f * dt * A;
    RHS_mat = HighFM::Eye<float>() + 0.5f * dt * A;

    HighFM::AMGSA<float, HighFM::SPAI0> P(LHS);
    HighFM::ConjugateGradient<float> cg(LHS.rows());

    float t = 0.0f;
    HighFM::Vector<float> b_n(dim), b_np1(dim), rhs(dim);
    for (size_t j = 0; j < dim; ++j) b_n(j) = _b[j]->operator()(t);

    for (size_t i = 0; i < _N; ++i) {
        t += dt;

        // Evaluate b(t_{n+1})
        for (size_t j = 0; j < dim; ++j) b_np1(j) = _b[j]->operator()(t);

        // Compute RHS
        rhs = RHS_mat * x;
        rhs = rhs + 0.5f * dt * (b_n + b_np1);

        // Solve system using conjugate gradient
        [[maybe_unused]] auto [steps, error] = HighFM::solve_linsys(P, LHS, x, rhs, cg);

        result.push_back(x(_row));
        b_n.swap(b_np1);  // Avoid copy, reuse buffer
    }

    return result;
}
