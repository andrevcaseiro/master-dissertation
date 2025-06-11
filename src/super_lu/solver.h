#pragma once

#include <Eigen/Sparse>

extern "C" {
#include <superlu_mt/slu_mt_sdefs.h>
}

class SuperLUSolver {
   public:
    typedef Eigen::SparseMatrix<float> SparseMatrixType;
    typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> DenseMatrixType;
    typedef Eigen::Matrix<float, Eigen::Dynamic, 1> VectorType;

    enum Status { Success, SingularMatrix, InvalidInput, MemoryAllocation };

    SuperLUSolver();
    ~SuperLUSolver();

    /**
     * @brief Sets the number of threads to use
     * @param nthreads Number of threads (default: 4)
     */
    void set_threads(int nthreads);

    /**
     * @brief Computes the sparse factorization of the matrix using SuperLU_MT
     * @param matrix The sparse matrix to factorize
     * @return Status of the computation
     */
    Status compute(SparseMatrixType& matrix);

    /**
     * @brief Solves the linear system Ax=b using the computed factorization
     * @param rhs The right hand side vector b
     * @return The solution vector x
     */
    VectorType solve(VectorType& rhs);

    /**
     * @brief Returns information about the last operation
     */
    Status info() const { return _info; }

   private:
    int _nthreads;
    superlumt_options_t _options;
    Gstat_t _Gstat;
    SuperMatrix _AC, _L, _U;
    int* _perm_r = nullptr;
    int* _perm_c = nullptr;

    Status _info = Success;
};
