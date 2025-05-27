#pragma once

#include <Eigen/Sparse>

#include "netlist.h"

class MNA {
   private:
   size_t _size;                                    // Number of nodes in the circuit
    Eigen::SparseMatrix<float> G;                    // Conductance matrix
    Eigen::DiagonalMatrix<float, Eigen::Dynamic> C;  // Capacitance matrix
    std::vector<std::unique_ptr<TimeFunction>> b;    // Independent sources vector

   public:
    /* Create MNA matrices from a netlist */
    MNA(const Netlist& netlist);

    /* Get the conductance matrix G */
    const Eigen::SparseMatrix<float>& get_G() const { return G; }

    /* Get the capacitance matrix C */
    const Eigen::DiagonalMatrix<float, Eigen::Dynamic>& get_C() const { return C; }

    /* Get the independent sources vector b */
    const std::vector<std::unique_ptr<TimeFunction>>& get_b() const { return b; }

    /* Get number of nodes */
    size_t size() const { return _size; }
};