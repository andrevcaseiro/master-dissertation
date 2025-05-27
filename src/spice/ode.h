#pragma once

#include <Eigen/Sparse>
#include <memory>
#include <vector>

#include "mna.h"
#include "time_function.h"

/**
 * @brief Converts MNA matrices to standard ODE form x'=Ax + b
 */
class ODE {
   private:
    Eigen::SparseMatrix<float> _A;
    std::vector<std::unique_ptr<TimeFunction>> _b;
    size_t _size;

   public:
    /**
     * @brief Construct ODE from MNA matrices
     */
    ODE(const MNA& mna);

    /**
     * @brief Get the A matrix where x'=Ax + b
     */
    const Eigen::SparseMatrix<float>& A() const { return _A; }

    /**
     * @brief Get the b vector where x'=Ax + b
     */
    const std::vector<std::unique_ptr<TimeFunction>>& b() const { return _b; }

    /**
     * @brief Get system size
     */
    size_t size() const { return _size; }
};
