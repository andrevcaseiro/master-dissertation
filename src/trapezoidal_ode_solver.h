#pragma once

#include <Eigen/Sparse>
#include <memory>

#include "spice/time_function.h"

class TrapezoidalODESolver {
   private:
    const Eigen::SparseMatrix<float>& _A;
    const std::vector<std::unique_ptr<TimeFunction>>& _b;
    const std::vector<float>& _x_0;
    float _t;
    size_t _row;
    size_t _N;

   public:
    TrapezoidalODESolver(const Eigen::SparseMatrix<float>& A,
                         const std::vector<std::unique_ptr<TimeFunction>>& b, 
                         const std::vector<float>& x_0,
                         float t, size_t row, size_t N);

    std::vector<float> solve_sequence();
};
