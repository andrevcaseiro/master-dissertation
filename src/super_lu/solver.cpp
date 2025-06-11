#include "solver.h"

#include <omp.h>

#include <cstdlib>

SuperLUSolver::SuperLUSolver() { _nthreads = omp_get_max_threads(); }

SuperLUSolver::~SuperLUSolver() {
    // Free permutation vectors
    SUPERLU_FREE(_perm_r);
    SUPERLU_FREE(_perm_c);

    // Destroy SuperMatrix structures
    if (_L.Store) Destroy_SuperNode_Matrix(&_L);
    if (_U.Store) Destroy_CompCol_Matrix(&_U);

    StatFree(&_Gstat);
    pxgstrf_finalize(&_options, &_AC);
}

void SuperLUSolver::set_threads(int nthreads) { _nthreads = nthreads; }

SuperMatrix to_super(SuperLUSolver::SparseMatrixType& matrix) {
    SuperMatrix A;

    // Get matrix dimensions and data
    int m = matrix.rows();
    int n = matrix.cols();
    int nnz = matrix.nonZeros();

    // Get raw pointers to Eigen's internal data
    float* values = matrix.valuePtr();
    int* row_indices = matrix.innerIndexPtr();
    int* col_ptrs = matrix.outerIndexPtr();

    // Create SuperMatrix in Compressed Column format (same as Eigen)
    sCreate_CompCol_Matrix(&A, m, n, nnz, values, row_indices, col_ptrs, SLU_NC, SLU_S, SLU_GE);

    return A;
}

SuperMatrix to_super(SuperLUSolver::VectorType& rhs) {
    SuperMatrix B;

    // Get vector dimensions
    int m = rhs.rows();
    int n = 1;  // Single column vector

    // Get raw pointer to Eigen's internal data
    float* values = rhs.data();

    // Create SuperMatrix in Dense format
    sCreate_Dense_Matrix(&B, m, n, values, m, SLU_DN, SLU_S, SLU_GE);

    return B;
}

SuperLUSolver::Status SuperLUSolver::compute(SparseMatrixType& matrix) {
    // Declare local variables for SuperLU
    SuperMatrix A = to_super(matrix);
    int n = A.ncol;

    fact_t fact = EQUILIBRATE;
    yes_no_t refact = NO;
    trans_t trans = NOTRANS;
    int panel_size = sp_ienv(1);
    int relax = sp_ienv(2);
    double diag_pivot_thresh = 1.0;
    yes_no_t usepr = NO;
    float drop_tol = 0.0;
    void* work = NULL;
    int lwork = 0;

    int nprocs = _nthreads;

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: natural ordering
     *   permc_spec = 1: minimum degree ordering on structure of A'*A
     *   permc_spec = 2: minimum degree ordering on structure of A'+A
     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     */
    _perm_c = intMalloc(n);
    _perm_r = intMalloc(n);
    int permc_spec = 1;
    get_perm_c(permc_spec, &A, _perm_c);

    StatAlloc(n, nprocs, panel_size, relax, &_Gstat);
    StatInit(n, nprocs, &_Gstat);

    psgstrf_init(nprocs, fact, trans, refact, panel_size, relax, diag_pivot_thresh, usepr, drop_tol,
                 _perm_c, _perm_r, work, lwork, &A, &_AC, &_options, &_Gstat);

    int info;
    psgstrf(&_options, &_AC, _perm_r, &_L, &_U, &_Gstat, &info);

    if (info == 0) {
        _info = Success;
    } else if (info < 0) {
        _info = InvalidInput;
    } else if (info >= A.ncol) {
        _info = MemoryAllocation;
    } else {
        _info = SingularMatrix;
    }

    free(A.Store);

    return _info;
}

SuperLUSolver::VectorType SuperLUSolver::solve(VectorType& rhs) {
    SuperMatrix B = to_super(rhs);

    int info;
    sgstrs(NOTRANS, &_L, &_U, _perm_r, _perm_c, &B, &_Gstat, &info);

    if (info == 0) {
        _info = Success;
    } else {
        _info = InvalidInput;
    }

    free(B.Store);

    return rhs;
}