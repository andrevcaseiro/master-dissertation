#pragma once

#include <Eigen/Sparse>
#include <vector>

#include "mna.h"

/**
 * @brief Solves the DC operating point of a circuit
 *
 * For a circuit described by Gx = b, finds the initial solution
 * when all time-varying sources are at t=0
 */
class DCSolver {
   private:
    const Eigen::SparseMatrix<float>& _G;
    const std::vector<std::unique_ptr<TimeFunction>>& _b;
    size_t _size;

    /**
     * @brief Create b vector at t=0
     */
    Eigen::VectorXf create_b0() const;

    /**
     * @brief Calculate normalized residual error ||Gx - b0|| / ||b0||
     */
    float calculate_residual_error(const Eigen::VectorXf& x) const;

    /**
     * @brief Solve using SparseLU decomposition
     */
    Eigen::VectorXf solve_LU() const;

    /**
     * @brief Solve using Conjugate Gradient method
     */
    Eigen::VectorXf solve_CG() const;

    /**
     * @brief Solve using SuperLU_MT
     */
    Eigen::VectorXf solver_SLU() const;

    /**
     * @brief Solve using HighFM Pardiso
     */
    Eigen::VectorXf solve_PARDISO() const;

   public:
    enum class Method {
        LU,      // Direct solver using LU decomposition
        CG,      // Iterative solver using Conjugate Gradient
        SLU,     // Direct solver using super LU
        PARDISO  // Direct solver using HighFM Pardiso
    };

    /**
     * @brief Construct DC analysis from MNA matrices
     *
     * @param mna MNA matrices
     */
    DCSolver(const MNA& mna) : _G(mna.get_G()), _b(mna.get_b()), _size(mna.size()) {}

    /**
     * @brief Solve the DC operating point
     *
     * @param method Solving method to use (default: PARDISO)
     * @return Eigen::VectorXf initial node voltages
     */
    Eigen::VectorXf solve(Method method = Method::PARDISO) const;
};
