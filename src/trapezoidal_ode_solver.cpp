#include "trapezoidal_ode_solver.h"

TrapezoidalODESolver::TrapezoidalODESolver(const Eigen::SparseMatrix<float>& A,
                                           const std::vector<std::unique_ptr<TimeFunction>>& b,
                                           const std::vector<float>& x_0, float t, size_t row, size_t N)
    : _A(A), _b(b), _x_0(x_0), _t(t), _row(row), _N(N) {}

std::vector<float> TrapezoidalODESolver::solve_sequence() {
    size_t dim = _x_0.size();
    float dt = _t / static_cast<float>(_N);

    Eigen::VectorXf x(Eigen::Map<const Eigen::VectorXf>(_x_0.data(), dim));

    // Initialize result vector with
    std::vector<float> result;
    result.reserve(_N + 1);
    result.push_back(x(_row));

    // Build identity matrix
    Eigen::SparseMatrix<float> I(dim, dim);
    I.setIdentity();

    // Precompute matrices
    Eigen::SparseMatrix<float> LHS = I - 0.5f * dt * _A;
    Eigen::SparseMatrix<float> RHS_mat = I + 0.5f * dt * _A;

    // Solver setup
    Eigen::SparseLU<Eigen::SparseMatrix<float>> solver;
    solver.compute(LHS);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Decomposition failed in Trapezoidal solver.");
    }

    float t = 0.0f;
    Eigen::VectorXf b_n(dim), b_np1(dim);
    for (size_t j = 0; j < dim; ++j) b_n(j) = _b[j]->operator()(t);

    for (size_t i = 0; i < _N; ++i) {
        t += dt;

        // Evaluate b(t_{n+1})
        for (size_t j = 0; j < dim; ++j) b_np1(j) = _b[j]->operator()(t);

        // Compute RHS
        Eigen::VectorXf rhs = RHS_mat * x + 0.5f * dt * (b_n + b_np1);

        // Solve system
        x = solver.solve(rhs);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("System solve failed.");
        }

        result.push_back(x(_row));
        b_n.swap(b_np1);  // Avoid copy, reuse buffer
    }

    return result;
}
