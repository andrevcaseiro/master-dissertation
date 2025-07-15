/**
 * @file highfm_ode_solver.h
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief A class to solve ODEs using the trapezoidal method through HighFM and PARDISO
 * @version 0.1
 * @date 2025-07-08
 *
 * @copyright Copyright (c) 2025
 *
 */

#pragma once

#include <highfm/core.hpp>
#include <highfm/linsys/pardiso.hpp>

#include "matrix/csrd_matrix.h"
#include "spice/time_function.h"

class HigFMODESolver {
   private:
    const Eigen::SparseMatrix<float>& _A;
    const std::vector<std::unique_ptr<TimeFunction>>& _b;
    const std::vector<float>& _x_0;
    float _t;
    size_t _row;
    size_t _N;

   public:
    /**
     * @brief ODE configuration
     *
     * @param A A
     * @param b b 
     * @param x_0 initial state
     * @param t time instant
     * @param row vector entry
     * @param N number of time steps
     */
    HigFMODESolver(const Eigen::SparseMatrix<float>& A,
                   const std::vector<std::unique_ptr<TimeFunction>>& b,
                   const std::vector<float>& x_0, float t, size_t row, size_t N);

    std::vector<float> solve_sequence() const;
};
