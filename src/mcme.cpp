/* Monte Carlo matrix exponential */

#include "mcme.h"

#include <math.h>

#include <iostream>
#include <random>
#include <vector>

template <typename Generator>
inline float MCME::generate_time(Generator& gen, int current_state) {
    std::exponential_distribution<float> dist(-_matrix.diagonal(current_state));

    return dist(gen);
}

template <typename Generator>
int MCME::generate_state(Generator& gen, int current_state) {
    auto row = _matrix.row(current_state);

    std::discrete_distribution<> dist(row.begin_off_diagonal_values(),
                                      row.end_off_diagonal_values());

    int column_index = dist(gen);

    int new_state = *(row.begin_off_diagonal_columns() + column_index);

    return new_state;
}

MCME::MCME(CSRMatrix<float> A, std::vector<float> x_0, float t, int i, int M, int N, int seed)
    : _matrix(A),
      _x_0(x_0),
      _t(t),
      _i(i),
      _M(M),
      _N(N),
      _seed(seed >= 0 ? seed : std::random_device{}()) {}

float MCME::calculate() {
    float delta_t = _t / _N;

    // Calculate D and L so that A = D - L and the sum of every row is 0
    std::vector<float> D(_matrix.rows());

    float res = 0;
#pragma omp parallel reduction(+ : res)
    {
#pragma omp for
        for (int row = 0; row < _matrix.rows(); ++row) {
            float sum = 0;
            for (auto it : _matrix.row(row)) {
                sum += it.value();
            }

            _matrix.diagonal(row) -= sum;
            D[row] = sum;
        }

#pragma omp for
        for (int m = 0; m < _M; m++) {
            /* Starting a generator on each line ensures reproduceability with different numbers of
             * threads */
            std::mt19937 gen(_seed + m);

            int state = _i;

            float sample = 1;
            for (int n = 1; n < _N; n++) {
                sample *= exp(D[state] * delta_t / 2);

                float tau = generate_time(gen, state);
                while (tau < delta_t) {
                    tau += generate_time(gen, state);
                    state = generate_state(gen, state);
                }

                sample *= exp(D[state] * delta_t / 2);
            }

            res += _x_0[state] * sample / _M;
        }
    }

    return res;
}

int w(int n, int size) {
    if (n == 0 || n == size) return 1;
    if (n % 2 == 0)
        return 4;
    else
        return 2;
}

float MCME::calculate_b(std::vector<float> b) {
    float delta_t = _t / _N;

    // Calculate D and L so that A = D - L and the sum of every row is 0
    std::vector<float> D(_matrix.rows());

    if (_N % 2 == 0) throw std::runtime_error("N must be an odd number");

    float res = 0;
#pragma omp parallel reduction(+ : res)
    {
#pragma omp for
        for (int row = 0; row < _matrix.rows(); ++row) {
            float sum = 0;
            for (auto it : _matrix.row(row)) {
                sum += it.value();
            }

            _matrix.diagonal(row) -= sum;
            D[row] = sum;
        }

#pragma omp for
        for (int m = 0; m < _M; m++) {
            /* Starting a generator on each line ensures reproduceability with different numbers of
             * threads */
            std::mt19937 gen(_seed + m);

            int state = _i;

            float sample = 1;
            float integral = b[state];
            for (int n = 1; n <= _N; n++) {
                sample *= exp(D[state] * delta_t / 2);

                float tau = generate_time(gen, state);
                while (tau < delta_t) {
                    tau += generate_time(gen, state);
                    state = generate_state(gen, state);
                }

                sample *= exp(D[state] * delta_t / 2);
                integral += w(n, _N) * sample * b[state];
            }

            res += (sample * _x_0[state] + integral * delta_t / 3) / _M;
        }
    }

    return res;
}
