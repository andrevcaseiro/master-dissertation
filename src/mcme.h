#pragma once
#include <random>

#include "matrix/csr_matrix.h"

class MCME {
   private:
    CSRMatrix<float> _matrix;
    std::vector<float> _x_0;
    float _t;
    int _i;
    int _M;
    int _N;
    int _seed;

    std::mt19937 _gen;

    float generate_time(int current_state);

    int generate_state(int current_state);

   public:
    /**
     * @brief Configures matrix exponential solver
     *
     * @param A Matrix
     * @param x_0 initial state
     * @param t Time instant
     * @param i Vector entry
     * @param M Number of samples
     * @param N Number of time steps
     * @return float The calculated value
     */
    MCME(CSRMatrix<float> A, std::vector<float> x_0, float t, int i, int M, int N, int seed);

    /**
     * @brief Calculates the result
     *
     * @return float
     */
    float calculate();
};
