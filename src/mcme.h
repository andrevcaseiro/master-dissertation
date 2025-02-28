#pragma once
#include <random>

#include "matrix/csrd_matrix.h"

class MCME {
   private:
    CSRMatrix<float> _matrix;
    std::vector<float> _x_0;
    float _t;
    int _i;
    int _M;
    int _N;
    int _seed;

    template <typename Generator>
    float generate_time(Generator& gen, int current_state);

    template <typename Generator>
    int generate_state(Generator& gen, int current_state);

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
