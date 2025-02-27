/* Monte Carlo matrix exponential */

#include "mcme.h"

#include <math.h>

#include <iostream>
#include <random>
#include <vector>

float MCME::generate_time(int current_state) {
    std::exponential_distribution<float> dist(-_matrix.diagonal(current_state));

    return dist(_gen);
}

int MCME::generate_state(int current_state) {
    auto row = _matrix.row(current_state);

    std::discrete_distribution<> dist(row.begin_off_diagonal_values(),
                                      row.end_off_diagonal_values());

    int column_index = dist(_gen);

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
      _seed(seed),
      _gen(seed >= 0 ? seed : std::random_device{}()) {}

float MCME::calculate() {
    float delta_t = _t / _N;

    // Calculate D and L so that A = D - L and the sum of every row is 0
    std::vector<float> D(_matrix.rows());
    for (int row = 0; row < _matrix.rows(); ++row) {
        float sum = 0;
        for (auto it : _matrix.row(row)) {
            sum += it.value();
        }

        _matrix.diagonal(row) -= sum;
        D[row] = sum;
    }

    float res = 0;

    for (int m = 0; m < _M; m++) {
        int state = _i;

        float sample = 1;
        for (int n = 1; n < _N; n++) {
            sample *= exp(D[state] * delta_t / 2);

            float tau = generate_time(state);
            while (tau < delta_t) {
                tau += generate_time(state);
                state = generate_state(state);
            }

            sample *= exp(D[state] * delta_t / 2);
        }

        res += _x_0[state] * sample / _M;
    }

    return res;
}
