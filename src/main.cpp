#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <vector>

#include "matrix/csrd_matrix.h"
#include "mcme.h"

/* Prints usage message to cout
 * Returns 0 */
int exec_help() {
    std::cout << "Usage: main command" << std::endl;
    std::cout << "\nCommands:" << std::endl;
    std::cout << "  help: displays this help message" << std::endl;
    std::cout << "  expm: calculates the matrix exponential" << std::endl;
    std::cout << "  expm-test: tests the accuracy of the matrix exponential on a set of tests"
              << std::endl;
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
        std::cout << "Usage: expm directory-path M N seed" << std::endl;
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

/**
 * @brief Calculates average accuracy of matrix exponential for a
 *
 * @param argc must be 4
 * @param argv (test directory path, M, N, seed)
 * @return int
 */
int exec_expm_test(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout << "Usage: expm-test matrix-filepath M N seed" << std::endl;
        return 1;
    }

    std::string directory_path(argv[0]);
    int M = atoi(argv[1]);
    int N = atoi(argv[2]);
    int seed = atoi(argv[3]);

    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::directory_iterator(directory_path)) {
        if (entry.is_regular_file()) {
            files.push_back(entry.path());
        }
    }

    std::sort(files.begin(), files.end());

    for (const auto& filepath : files) {
        std::string filename = filepath.filename().string();
        std::string stem = filepath.stem().string();
        std::string extension = filepath.extension().string();

        std::string solution_filename = stem + "_exp" + extension;
        std::filesystem::path solution_path = filepath.parent_path() / solution_filename;

        if (std::filesystem::exists(solution_path)) {
            CSRMatrix<float> matrix(filepath);
            CSRMatrix<float> exp_matrix(solution_path);

            int row = 0;
            int col = 0;
            std::vector<float> v(matrix.columns(), 0);
            v[col] = 1;

            MCME mcme(matrix, v, 1, row, M, N, seed);
            float res = mcme.calculate();
            float sol = exp_matrix.at(row, col);

            float error = std::abs(res - sol) / std::abs(sol);
            std::cout << filepath.string() << ": " << std::setprecision(5) << error * 100 << "%"
                      << std::endl;
        }
    }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        exec_help();
        return 1;
    }

    std::string command = argv[1];

    if (command == "help") {
        return exec_help();
    } else if (command == "expm") {
        return exec_expm(argc - 2, &argv[2]);
    } else if (command == "expm-test") {
        return exec_expm_test(argc - 2, &argv[2]);
    } else {
        exec_help();
        return 1;
    }

    return 0;
}
