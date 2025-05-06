#pragma once

#include <Eigen/Sparse>
#include <memory>

#include "spice/time_function.h"

class TrapezoidalODESolver {
   private:
    Eigen::SparseMatrix<float>& _A;
    std::vector<std::unique_ptr<TimeFunction>>& _b;
    std::vector<float>& _x_0;
    float _t;
    size_t _row;
    size_t _N;

   public:
    TrapezoidalODESolver(Eigen::SparseMatrix<float>& A,
                         std::vector<std::unique_ptr<TimeFunction>>& b, std::vector<float>& x_0,
                         float t, size_t row, size_t N);

    std::vector<float> solve_sequence();
};
