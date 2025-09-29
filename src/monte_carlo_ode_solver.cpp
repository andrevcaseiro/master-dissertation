/**
 * @file MonteCarloODESolver.cpp
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief Monte Carlo Ordinary Diferential Equation solver
 * @version 0.1
 * @date 2025-03-18
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "monte_carlo_ode_solver.h"

#include <omp.h>

#include <iostream>
#include <iterator>
#include <random>
#include <stdexcept>

#include "matrix/csrd_matrix.h"
#include "utils/progress_bar.h"

template <typename Generator>
float MonteCarloODESolver::generate_time(Generator& gen, size_t current_state) {
    /* Exponential distribution with positive lambda */
    std::exponential_distribution<float> dist(-_L.diagonal(current_state));

    return dist(gen);
}

template <typename Generator>
size_t MonteCarloODESolver::generate_state(Generator& gen, size_t current_state) {
    auto row = _L.row(current_state);

    std::discrete_distribution<size_t> dist(row.begin_off_diagonal_values(),
                                            row.end_off_diagonal_values());

    size_t column_index = dist(gen);

    size_t new_state = *(row.begin_off_diagonal_columns() + column_index);

    return new_state;
}

int MonteCarloODESolver::w(size_t n) const {
    if (n == 0 || n == _N) return 1;
    if (n % 2 == 0) return 2;
    return 4;
}

MonteCarloODESolver::MonteCarloODESolver(CSRMatrix<float>& A, const std::vector<float>& b,
                                         std::vector<float>& x_0, float t, size_t row, size_t M,
                                         size_t N, long seed)
    : _A(A),
      _x_0(x_0),
      _t(t),
      _row(row),
      _M(M),
      _N(N),
      _seed(seed >= 0 ? seed : std::random_device{}()),
      _init(false) {
    for (auto value : b) {
        _b.push_back(std::make_unique<ConstantFunction>(value));
    }
}

MonteCarloODESolver::MonteCarloODESolver(CSRMatrix<float>& A,
                                         const std::vector<std::unique_ptr<TimeFunction>>& b,
                                         std::vector<float>& x_0, float t, size_t row, size_t M,
                                         size_t N, long seed)
    : _A(A),
      _x_0(x_0),
      _t(t),
      _row(row),
      _M(M),
      _N(N),
      _seed(seed >= 0 ? seed : std::random_device{}()),
      _init(false) {
    _b.reserve(b.size());
    for (const auto& value : b) {
        _b.push_back(value->clone());
    }
}

void MonteCarloODESolver::init() {
    if (_init) return;

    _L = _A;
    _D = std::vector<float>(_L.rows());

#pragma omp parallel for
    for (size_t row = 0; row < _L.rows(); ++row) {
        float sum = 0;
        for (auto it : _L.row(row)) {
            sum += it.value();
        }

        _L.diagonal(row) -= sum;
        _D[row] = sum;
    }

    _init = true;
}

float MonteCarloODESolver::solve() {
    /* Ensure N is odd for Simpson quadrature */
    if (_N % 2 != 0) throw std::runtime_error("N must be an even number");

    float delta_t = _t / _N;

    init();

    float res = 0;
#pragma omp parallel reduction(+ : res)
    {
        /* Constructing one generator per thread generates different results */
        std::mt19937 gen(_seed + omp_get_thread_num());

#pragma omp for
        for (size_t m = 0; m < _M; m++) {
            /* Samples are divided equally among threads */

            size_t state = _row;
            float exponential = 1;
            float integral = (*_b[state])(_t);
            for (size_t n = 1; n <= _N; n++) {
                exponential *= exp(_D[state] * delta_t / 2);

                float tau = generate_time(gen, state);
                while (tau < delta_t) {
                    tau += generate_time(gen, state);
                    state = generate_state(gen, state);
                }

                exponential *= exp(_D[state] * delta_t / 2);
                integral += w(n) * exponential * (*_b[state])(_t - n * delta_t);
            }

            float sample = exponential * _x_0[state] + integral * delta_t / 3;
            res += sample / _M;
        }
    }

    return res;
}

std::vector<float> MonteCarloODESolver::solve_sequence(size_t output_N) {
    double delta_t = _t / _N;

    if (output_N == 0) output_N = _N;
    if (_N % output_N != 0) {
        std::runtime_error(
            "The number of output points must evenly divide the total number of steps.");
    }
    size_t output_freq = _N / output_N;
    size_t output_size = output_N + 1;
    double output_delta_t = _t / output_N;

    init();

    std::vector<float> res(output_size, 0);

    // Progress bar - only master thread updates it
    ProgressBar progress("MC Samples", _M);

#pragma omp parallel
    {
        /* Constructing one generator per thread generates different results */
        std::mt19937 gen(_seed + omp_get_thread_num());

        std::vector<float> local_res(output_size, 0);

        std::vector<float> integrals(output_size, 0);
#pragma omp for
        for (size_t m = 0; m < _M; m++) {
            // Only master thread updates progress
            if (omp_get_thread_num() == 0) {
                // Estimate progress based on master thread's work
                size_t estimated_completed = m * omp_get_num_threads();
                progress.update(estimated_completed);
            }
            
            /* Samples are divided equally among threads */

            size_t state = _row;
            float exponential = 1;
            for (size_t output_n = 1; output_n <= output_N; ++output_n) {
                integrals[output_n] = 0.5 * (*_b[state])(output_n * output_delta_t);
            }
            
            for (size_t n = 1; n <= _N; ++n) {
                exponential *= exp(_D[state] * delta_t / 2);

                float tau = generate_time(gen, state);
                while (tau < delta_t) {
                    tau += generate_time(gen, state);
                    state = generate_state(gen, state);
                }

                exponential *= exp(_D[state] * delta_t / 2);

                if (n % output_freq == 0) {
                    size_t output_n = n / output_freq;
                    for (size_t curr_output_n = output_n; curr_output_n <= output_N; ++curr_output_n) {
                        float sample_time = (curr_output_n - output_n) * output_delta_t;
                        if (curr_output_n == output_n) {
                            integrals[curr_output_n] += exponential * 0.5 * (*_b[state])(sample_time);
                        } else {
                            integrals[curr_output_n] += exponential * (*_b[state])(sample_time);
                        }
                    }

                    float sample = exponential * _x_0[state] + integrals[output_n] * output_delta_t;
                    local_res[output_n] += sample / _M;
                }
            }
        }

#pragma omp critical
        for (size_t i = 0; i < res.size(); i++) {
            /* Reduce local results */
            res[i] += local_res[i];
        }
    }

    progress.complete();

    /* Set res[0] directly to avoid float roundings */
    res[0] = _x_0[_row];

    return res;
}
