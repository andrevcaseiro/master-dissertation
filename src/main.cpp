#include <iomanip>
#include <iostream>
#include <vector>

#include "matrix/csr_matrix.h"
#include "mcme.h"

/* Prints usage message to cout
 * Returns 0 */
int exec_help() {
    std::cout << "Usage: main command" << std::endl;
    std::cout << "\nCommands:" << std::endl;
    std::cout << "  help: displays this help message" << std::endl;
    std::cout << "  expm: tests matrix exponential" << std::endl;
    return 0;
}

/**
 * @brief Calculates the matrix exponential
 *
 * @param argc contains cli arguments
 * @param argv contains the number of cli arguments
 * @return int 0 on success
 */
int exec_expm(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout << "Usage: expm matrix-filepath M N seed" << std::endl;
        return 1;
    }

    std::string matrix_filepath(argv[0]);
    int M = atoi(argv[1]);
    int N = atoi(argv[2]);
    int seed = atoi(argv[3]);

    CSRMatrix<float> matrix(matrix_filepath);
    std::cout << "Full matrix:" << std::endl;
    matrix.print();
    std::cout << std::endl;
    std::cout << "Sparce matrix:" << std::endl;
    matrix.print_csr();
    std::cout << std::endl;

    std::cout << "Matrix exponential:" << std::endl;
    for (int row = 0; row < matrix.rows(); row++) {
        for (int col = 0; col < matrix.columns(); col++) {
            std::vector<float> v(matrix.columns(), 0);
            v[col] = 1;

            MCME mcme(matrix, v, 1, row, M, N, seed);
            std::cout << std::setw(10) << std::fixed
                      << std::setprecision(5)  // Adjust width & precision
                      << mcme.calculate() << " ";
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
    } else if (command == "expm") {
        return exec_expm(argc - 2, &argv[2]);
    } else {
        exec_help();
        return 1;
    }

    return 0;
}
