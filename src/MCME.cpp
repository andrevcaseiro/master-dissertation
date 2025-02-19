/* Monte Carlo matrix exponential */

#include "MCME.h"

#include <math.h>

#include <iostream>
#include <random>
#include <vector>

#include "matrix/CSRMatrix.h"

float MCME::generate_time(int current_state) {
    std::exponential_distribution<float> dist(-A.diagonal(current_state));

    return dist(gen);
}

int MCME::generate_state(int current_state) {
    auto it = A.row_iterator(current_state);

    std::vector<float> weights;
    std::vector<int> indices;
    while (it.first != it.second) {
        if (it.first->column != current_state) {
            weights.emplace_back(it.first->value);
            indices.emplace_back(it.first->column);
        }

        it.first++;
    }

    // TODO: avoid creating vectors every iteration (iterator?)

    std::discrete_distribution<> dist(weights.begin(), weights.end());

    auto res = dist(gen);

    return indices[res];
}

MCME::MCME(CSRMatrix<float> A, std::vector<float> x_0, float t, int i, int M, int N, int seed)
    : A(A),
      x_0(x_0),
      t(t),
      i(i),
      M(M),
      N(N),
      seed(seed),
      gen(seed >= 0 ? seed : std::random_device{}()) {}

float MCME::calculate() {
    float delta_t = t / N;

    // TODO: dont assume D=[1]
    std::vector<float> D({-1, 0, -1});

    float res = 0;

    for (int m = 0; m < M; m++) {
        int state = i;

        float sample = 1;
        for (int n = 1; n < N; n++) {
            sample *= exp(D[state] * delta_t / 2);

            float tau = generate_time(state);
            while (tau < delta_t) {
                //std::cout << tau << " ";
                tau += generate_time(state);
                state = generate_state(state);
            }
            //std::cout << std::endl;

            sample *= exp(D[state] * delta_t / 2);
        }

        res += x_0[state] * sample / M;
    }

    return res;
}
