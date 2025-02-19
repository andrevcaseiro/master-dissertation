#include <iostream>
#include <vector>

#include "MCME.h"
#include "matrix/CSRMatrix.h"

/* Prints usage message to cout
 * Returns 0 */
int exec_help() {
    std::cout << "Usage: main command" << std::endl;
    std::cout << "\nCommands:" << std::endl;
    std::cout << "help: displays this help message" << std::endl;
    std::cout << "test: tests matrix exponential" << std::endl;
    return 0;
}

/**
 * @brief Calculates the matrix exponential
 *
 * @param argc contains cli arguments
 * @param argv contains the number of cli arguments
 * @return int 0 on success
 */
int exec_test(int argc, char* argv[]) {
    if (argc != 5) {
        std::cout << "Usage: matrix-filepath entry M N seed" << std::endl;
        return 1;
    }

    std::string matrix_filepath(argv[0]);
    int entry = atoi(argv[1]);
    int M = atoi(argv[2]);
    int N = atoi(argv[3]);
    int seed = atoi(argv[4]);

    CSRMatrix<float> matrix(matrix_filepath);
    matrix.print();
    matrix.print_csr();

    for (int row = 0; row < matrix.get_rows(); row++) {
        std::vector<float> v(matrix.get_columns(), 0);
        v[row] = 1;

        /* for(auto it: v) {
            std::cout << it;
        }
        std::cout << std::endl; */

        for (int col = 0; col < matrix.get_columns(); col++) {
            MCME mcme(matrix, v, 1, col, M, N, seed);
            std::cout << mcme.calculate() << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        exec_help();
        return 0;
    }

    std::string command = argv[1];

    if (command == "help") {
        return exec_help();
    } else if (command == "test") {
        return exec_test(argc - 2, &argv[2]);
    }

    return 0;
}
