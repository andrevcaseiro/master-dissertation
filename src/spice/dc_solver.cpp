#include "dc_solver.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include "../super_lu/solver.h"
#include "../highfm/pardiso_solver.h"

Eigen::VectorXf DCSolver::solve(Method method) const {
    switch (method) {
        case Method::LU:
            return solve_LU();
        case Method::CG:
            return solve_CG();
        case Method::SLU:
            return solver_SLU();
        case Method::PARDISO:
            return solve_PARDISO();
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

Eigen::VectorXf DCSolver::solver_SLU() const {
    // Convert b vector to values at t=0
    Eigen::VectorXf b0(_size);
    for (size_t i = 0; i < _size; i++) {
        b0[i] = _b[i] ? (*_b[i])(0) : 0.0f;
    }

    // Create a copy of the G matrix for SuperLU
    SuperLUSolver::SparseMatrixType G_copy = _G;

    // Use SuperLU to solve the system
    SuperLUSolver solver;
    SuperLUSolver::Status status = solver.compute(G_copy);
    if (status != SuperLUSolver::Success) {
        throw std::runtime_error("Failed to decompose G matrix using SuperLU for DC analysis");
    }

    auto x = solver.solve(b0);
    if (solver.info() != SuperLUSolver::Success) {
        throw std::runtime_error("Failed to solve DC analysis using SuperLU");
    }

    return x;
}

Eigen::VectorXf DCSolver::solve_PARDISO() const {
    // Convert b vector to values at t=0
    Eigen::VectorXf b0(_size);
    for (size_t i = 0; i < _size; i++) {
        b0[i] = _b[i] ? (*_b[i])(0) : 0.0f;
    }

    // Use HighFM Pardiso to solve the system
    PardisoSolver solver;
    PardisoSolver::Status status = solver.compute(_G);
    if (status != PardisoSolver::Success) {
        throw std::runtime_error("Failed to decompose G matrix using Pardiso for DC analysis");
    }

    auto x = solver.solve(b0);
    if (solver.info() != PardisoSolver::Success) {
        throw std::runtime_error("Failed to solve DC analysis using Pardiso");
    }

    return x;
}
