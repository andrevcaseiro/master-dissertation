/**
 * @file pardiso_solver.cpp
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief A wrapper class for HighFM Pardiso solver with Eigen interface
 * @version 0.1
 * @date 2025-07-25
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "pardiso_solver.h"
#include <iostream>

PardisoSolver::PardisoSolver() { _info = Success; }

PardisoSolver::~PardisoSolver() = default;

PardisoSolver::Status PardisoSolver::compute(const SparseMatrixType& matrix) {
    try {
        // Reset state
        _is_factorized = false;
        _info = Success;

        // Check matrix validity
        if (matrix.rows() != matrix.cols()) {
            _info = InvalidInput;
            return _info;
        }

        if (matrix.rows() == 0) {
            _info = InvalidInput;
            return _info;
        }

        _dim = matrix.rows();

        // Convert matrix from csc to csr
        _csr_matrix = matrix;

        // Create a HighFM map to the eigen matrix without copying
        _A = HighFM::CSRMap<float>(
            HighFM::CSRMatrixData<float>{.nrows = static_cast<HighFM::index_t>(_dim),
                                         .ncols = static_cast<HighFM::index_t>(_dim),
                                         .nnz = _csr_matrix.nonZeros(),
                                         .row_offset = _csr_matrix.outerIndexPtr(),
                                         .nz_columns = _csr_matrix.innerIndexPtr(),
                                         .nz_values = _csr_matrix.valuePtr()});

        // Perform factorization
        _solver.factorize(_A);
        _is_factorized = true;
        _info = Success;
    } catch (const std::exception& e) {
        std::cerr << "Pardiso factorization error: " << e.what() << std::endl;
        _info = SingularMatrix;  // Most common error in sparse factorization
        _is_factorized = false;
    } catch (...) {
        std::cerr << "Unknown error during Pardiso factorization" << std::endl;
        _info = MemoryAllocation;
        _is_factorized = false;
    }

    return _info;
}

PardisoSolver::VectorType PardisoSolver::solve(const VectorType& rhs) {
    VectorType solution;
    solution.resize(_dim);

    HighFM::VectorMap<float> highfm_solution(HighFM::VectorData<float>{
        .buffer = solution.data(), .size = static_cast<HighFM::index_t>(_dim)});
    try {
        // Check if factorization was performed
        if (!_is_factorized) {
            _info = NotFactorized;
            return VectorType::Zero(rhs.size());
        }

        // Check RHS vector size
        if (static_cast<size_t>(rhs.size()) != _dim) {
            _info = InvalidInput;
            return VectorType::Zero(rhs.size());
        }

        // Convert Eigen vector to HighFM vector
        HighFM::VectorMap<float> highfm_rhs(
            HighFM::VectorData<float>{.buffer = const_cast<float*>(rhs.data()),
                                      .size = static_cast<HighFM::index_t>(rhs.size())});

        // Solve the system
        highfm_solution = _solver.solve(highfm_rhs);

        _info = Success;

    } catch (const std::exception& e) {
        std::cerr << "Pardiso solve error: " << e.what() << std::endl;
        _info = SingularMatrix;
        solution = VectorType::Zero(rhs.size());
    } catch (...) {
        std::cerr << "Unknown error during Pardiso solve" << std::endl;
        _info = MemoryAllocation;
        solution = VectorType::Zero(rhs.size());
    }

    return solution;
}
