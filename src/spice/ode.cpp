/**
 * @file ode.cpp
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief A class to represent an ODE system from an MNA system
 * @version 0.1
 * @date 2025-10-24
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include "ode.h"

ODE::ODE(const MNA& mna) : _size(mna.size()) {
    // Get inverse of C matrix
    Eigen::VectorXf c_inv = mna.get_C().diagonal().cwiseInverse();

    // Create diagonal matrix from vector
    Eigen::DiagonalMatrix<float, Eigen::Dynamic> C_inv(c_inv);

    // Compute A = -C^-1 * G
    _A = -(C_inv * mna.get_G());

    // Convert b vector
    _b.reserve(_size);
    const auto& mna_b = mna.get_b();

    for (size_t i = 0; i < _size; i++) {
        if (mna_b[i]) {
            auto b_clone = mna_b[i]->clone();
            *b_clone *= c_inv(i);
            _b.push_back(std::move(b_clone));
        } else {
            _b.push_back(std::make_unique<ConstantFunction>(0.0f));
        }
    }
}
