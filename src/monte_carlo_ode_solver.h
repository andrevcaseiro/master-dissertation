#pragma once

#include <memory>

#include "matrix/csrd_matrix.h"
#include "spice/time_function.h"

class MonteCarloODESolver {
   private:
    CSRMatrix<float>& _A;
    std::vector<float>& _x_0;
    std::vector<std::unique_ptr<TimeFunction>> _b;
    std::vector<float> _D;
    CSRMatrix<float> _L;
    float _t;
    size_t _row;
    size_t _M;
    size_t _N;
    long _seed;
    bool _init;

    /**
     * @brief Generate an exponentially distributed time
     *
     * @tparam Generator random generator type
     * @param gen generator
     * @param current_state current state
     * @return float exponentially distributed time
     */
    template <typename Generator>
    float generate_time(Generator& gen, size_t current_state);

    /**
     * @brief Generate a random new state
     *
     * @tparam Generator random generator type
     * @param gen generator
     * @param current_state current state
     * @return size_t new state
     */
    template <typename Generator>
    size_t generate_state(Generator& gen, size_t current_state);

    /**
     * @brief Simpson quadrature rule weights
     *
     * @param n current iteration
     * @return int 1,4,2,4,2,...,4,1
     */
    int w(size_t n) const;

   public:
    /**
     * @brief ODE configuration
     *
     * @param A A
     * @param b b
     * @param x_0 initial state
     * @param t Time instant
     * @param row Vector entry
     * @param M Number of samples
     * @param N Number of time steps
     * @return float The calculated value
     */
    MonteCarloODESolver(CSRMatrix<float>& A, std::vector<float>& b, std::vector<float>& x_0,
                        float t, size_t row, size_t M, size_t N, long seed);

    /**
     * @brief ODE configuration
     *
     * @param A A
     * @param b b
     * @param x_0 initial state
     * @param t Time instant
     * @param row Vector entry
     * @param M Number of samples
     * @param N Number of time steps
     * @return float The calculated value
     */
    MonteCarloODESolver(CSRMatrix<float>& A, std::vector<std::unique_ptr<TimeFunction>>& b,
                        std::vector<float>& x_0, float t, size_t row, size_t M, size_t N,
                        long seed);

    /**
     * @brief Prepares structures to solve
     *
     */
    void init();

    /**
     * @brief Calculates the result
     *
     * @return float
     */
    float solve();

    /**
     * @brief Calculates the value at each delta t
     *
     * @return std::vector<float> vector of size N+1 with x(n*delta t) at position n
     */
    std::vector<float> solve_sequence();
};
