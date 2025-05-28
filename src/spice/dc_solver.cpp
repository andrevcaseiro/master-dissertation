#include "dc_solver.h"

#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>

Eigen::VectorXf DCSolver::solve(Method method) const {
    switch (method) {
        case Method::LU:
            return solve_LU();
        case Method::CG:
            return solve_CG();
        default:
            throw std::runtime_error("Unknown solver method");
    }
}

Eigen::VectorXf DCSolver::solve_LU() const {
    // Convert b vector to values at t=0
    Eigen::VectorXf b0(_size);
    for (size_t i = 0; i < _size; i++) {
        b0[i] = _b[i] ? (*_b[i])(0) : 0.0f;
    }

    // Solve system Gx = b(0)
    Eigen::SparseLU<Eigen::SparseMatrix<float>> solver;
    solver.compute(_G);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Failed to decompose G matrix for DC analysis");
    }

    auto x = solver.solve(b0);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Failed to solve DC analysis");
    }

    return x;
}

Eigen::VectorXf DCSolver::solve_CG() const {
    // Convert b vector to values at t=0
    Eigen::VectorXf b0(_size);
    for (size_t i = 0; i < _size; i++) {
        b0[i] = _b[i] ? (*_b[i])(0) : 0.0f;
    }

    // Solve system Gx = b(0) using Conjugate Gradient
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> solver;
    solver.setMaxIterations(1000);
    solver.setTolerance(1e-10);
    solver.compute(_G);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Failed to initialize CG solver for DC analysis");
    }

    auto x = solver.solve(b0);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Failed to solve DC analysis using CG method");
    }

    return x;
}
