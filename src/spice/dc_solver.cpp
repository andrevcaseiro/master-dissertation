#include "dc_solver.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <iostream>
#include <iomanip>
#include "../super_lu/solver.h"
#include "../highfm/pardiso_solver.h"

Eigen::VectorXf DCSolver::create_b0() const {
    Eigen::VectorXf b0(_size);
    for (size_t i = 0; i < _size; i++) {
        b0[i] = _b[i] ? (*_b[i])(0) : 0.0f;
    }
    return b0;
}

float DCSolver::calculate_residual_error(const Eigen::VectorXf& x) const {
    // Create fresh b0 vector (solvers may have modified the original)
    Eigen::VectorXf b0 = create_b0();
    
    // Calculate Gx
    Eigen::VectorXf Gx = _G * x;
    
    // Calculate residual: r = Gx - b0
    Eigen::VectorXf residual = Gx - b0;

    // Print G*x and b0 side by side for comparison
    std::cout << "\nDC Analysis Verification:" << std::endl;
    std::cout << std::setw(6) << "Index" << std::setw(15) << "x (solution)" << std::setw(15) << "G*x" 
              << std::setw(15) << "b0" << std::setw(15) << "Residual" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    for (Eigen::Index i = 0; i < x.size(); i++) {
        std::cout << std::setw(6) << i 
                  << std::setw(15) << std::scientific << std::setprecision(6) << x[i]
                  << std::setw(15) << std::scientific << std::setprecision(6) << Gx[i]
                  << std::setw(15) << std::scientific << std::setprecision(6) << b0[i]
                  << std::setw(15) << std::scientific << std::setprecision(6) << residual[i]
                  << std::endl;
    }
    
    std::cout << std::endl;
    
    // Calculate normalized error: ||r|| / ||b0||
    float residual_norm = residual.norm();
    float b0_norm = b0.norm();
    
    float normalized_error = (b0_norm > 1e-15f) ? residual_norm / b0_norm : residual_norm;
    
    std::cout << "DC solver residual error: " << std::scientific << std::setprecision(6) 
              << normalized_error << std::endl << std::endl;
    
    return normalized_error;
}

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
    Eigen::VectorXf b0 = create_b0();

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
    Eigen::VectorXf b0 = create_b0();

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
    Eigen::VectorXf b0 = create_b0();

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
    Eigen::VectorXf b0 = create_b0();

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
