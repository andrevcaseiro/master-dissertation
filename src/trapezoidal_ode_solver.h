#pragma once

#include <Eigen/Sparse>
#include <memory>

#include "spice/time_function.h"

class TrapezoidalODESolver {
   public:
    enum class Method {
        LU,         // Direct solver using SparseLU
        CG,         // Iterative solver using Conjugate Gradient
        SUPERLU_MT  // Parallel direct solver using SuperLU_MT
    };

   private:
    const Eigen::SparseMatrix<float>& _A;
    const std::vector<std::unique_ptr<TimeFunction>>& _b;
    const std::vector<float>& _x_0;
    float _t;
    size_t _row;
    size_t _N;

    /**
     * @brief Solve using SparseLU decomposition
     */
    std::vector<float> solve_sequence_LU() const;

    /**
     * @brief Solve using Conjugate Gradient method
     */
    std::vector<float> solve_sequence_CG() const;

    /**
     * @brief Solve using SuperLU_MT parallel direct solver
     */
    std::vector<float> solve_sequence_superlu_mt() const;

   public:
    TrapezoidalODESolver(const Eigen::SparseMatrix<float>& A,
                         const std::vector<std::unique_ptr<TimeFunction>>& b,
                         const std::vector<float>& x_0, float t, size_t row, size_t N);

    /**
     * @brief Solve the ODE using the specified method
     *
     * @param method Solving method to use (default: CG)
     * @return std::vector<float> Solution sequence
     */
    std::vector<float> solve_sequence(Method method = Method::CG) const;
};
