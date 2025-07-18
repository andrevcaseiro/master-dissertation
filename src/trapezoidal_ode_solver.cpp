#include "trapezoidal_ode_solver.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include "super_lu/solver.h"

TrapezoidalODESolver::TrapezoidalODESolver(const Eigen::SparseMatrix<float>& A,
                                           const std::vector<std::unique_ptr<TimeFunction>>& b,
                                           const std::vector<float>& x_0, float t, size_t row,
                                           size_t N)
    : _A(A), _b(b), _x_0(x_0), _t(t), _row(row), _N(N) {}

std::vector<float> TrapezoidalODESolver::solve_sequence(Method method) const {
    switch (method) {
        case Method::LU:
            return solve_sequence_LU();
        case Method::CG:
            return solve_sequence_CG();
        case Method::SUPERLU_MT:
            return solve_sequence_superlu_mt();
        default:
            throw std::runtime_error("Unknown solver method");
    }
}

std::vector<float> TrapezoidalODESolver::solve_sequence_LU() const {
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
    Eigen::SparseMatrix<float, Eigen::RowMajor> RHS_mat = I + 0.5f * dt * _A;

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

std::vector<float> TrapezoidalODESolver::solve_sequence_CG() const {
    size_t dim = _x_0.size();
    float dt = _t / static_cast<float>(_N);

    Eigen::VectorXf x(Eigen::Map<const Eigen::VectorXf>(_x_0.data(), dim));

    // Initialize result vector
    std::vector<float> result;
    result.reserve(_N + 1);
    result.push_back(x(_row));

    // Build identity matrix
    Eigen::SparseMatrix<float> I(dim, dim);
    I.setIdentity();

    // Precompute matrices
    Eigen::SparseMatrix<float> LHS = I - 0.5f * dt * _A;
    Eigen::SparseMatrix<float, Eigen::RowMajor> RHS_mat = I + 0.5f * dt * _A;

    // Solver setup
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> solver;
    solver.setMaxIterations(1000);
    solver.setTolerance(1e-10);
    solver.compute(LHS);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Failed to initialize CG solver");
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
            throw std::runtime_error("Failed to solve system using CG method");
        }

        result.push_back(x(_row));
        b_n.swap(b_np1);  // Avoid copy, reuse buffer
    }

    return result;
}

std::vector<float> TrapezoidalODESolver::solve_sequence_superlu_mt() const {
    size_t dim = _x_0.size();
    float dt = _t / static_cast<float>(_N);

    Eigen::VectorXf x(Eigen::Map<const Eigen::VectorXf>(_x_0.data(), dim));

    // Initialize result vector
    std::vector<float> result;
    result.reserve(_N + 1);
    result.push_back(x(_row));

    // Build identity matrix
    Eigen::SparseMatrix<float> I(dim, dim);
    I.setIdentity();

    // Precompute matrices
    Eigen::SparseMatrix<float> LHS = I - 0.5f * dt * _A;
    Eigen::SparseMatrix<float, Eigen::RowMajor> RHS_mat = I + 0.5f * dt * _A;

    // Solver setup
    SuperLUSolver solver;
    solver.compute(LHS);
    if (solver.info() != SuperLUSolver::Success) {
        throw std::runtime_error("Failed to initialize SLU solver");
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
        if (solver.info() != SuperLUSolver::Success) {
            throw std::runtime_error("Failed to solve system using SLU method");
        }

        result.push_back(x(_row));
        b_n.swap(b_np1);  // Avoid copy, reuse buffer
    }

    return result;
}