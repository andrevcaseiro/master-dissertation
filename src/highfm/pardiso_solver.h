#pragma once

#include <Eigen/Sparse>
#include <highfm/core.hpp>
#include <highfm/linsys/pardiso.hpp>

class PardisoSolver {
   public:
    typedef Eigen::SparseMatrix<float> SparseMatrixType;
    typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> DenseMatrixType;
    typedef Eigen::Matrix<float, Eigen::Dynamic, 1> VectorType;

    enum Status { Success, SingularMatrix, InvalidInput, MemoryAllocation, NotFactorized };

    PardisoSolver();
    ~PardisoSolver();

    /**
     * @brief Computes the sparse factorization of the matrix using HighFM Pardiso
     * @param matrix The sparse matrix to factorize
     * @return Status of the computation
     */
    Status compute(const SparseMatrixType& matrix);

    /**
     * @brief Solves the linear system Ax=b using the computed factorization
     * @param rhs The right hand side vector b
     * @return The solution vector x
     */
    VectorType solve(const VectorType& rhs);

    /**
     * @brief Returns information about the last operation
     */
    Status info() const { return _info; }

    /**
     * @brief Check if the matrix has been factorized
     */
    bool is_factorized() const { return _is_factorized; }

   private:
    HighFM::Pardiso<float> _solver;
    Status _info = Success;
    bool _is_factorized = false;
    size_t _dim = 0;
    Eigen::SparseMatrix<float, Eigen::RowMajor, HighFM::index_t> _csr_matrix;
    HighFM::CSRMatrix<float> _A;
};
