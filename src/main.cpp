#include <omp.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "matrix/csrd_matrix.h"
#include "mcme.h"

std::vector<float> read_vector(std::string filepath) {
    std::vector<float> v;
    std::ifstream file(filepath);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open the file.");
    }

    std::string line;

    // Read the length
    if (!std::getline(file, line)) {
        throw std::runtime_error("File is empty.");
    }
    int expected_length = 0;
    std::istringstream(line) >> expected_length;
    v.reserve(expected_length);

    // Read the second line for the comma-separated values
    if (!std::getline(file, line)) {
        throw std::runtime_error("Vector data missing.");
    }

    std::istringstream ss(line);
    std::string value;
    while (std::getline(ss, value, ',')) {
        try {
            v.push_back(std::stof(value));  // Convert string to float and add to vector
        } catch (const std::invalid_argument& e) {
            throw std::runtime_error("Invalid float value in the CSV file.");
        }
    }

    // Check if the number of elements matches the expected length
    if (v.size() != (size_t) expected_length) {
        throw std::runtime_error("Vector data does not match the expected length.");
    }

    return v;
}

/* Prints usage message to cout
 * Returns 0 */
int exec_help() {
    std::cout << "Usage: main command" << std::endl;
    std::cout << std::endl;
    std::cout << "Commands:" << std::endl;
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

/**
 * @brief Calculates average accuracy of matrix exponential for a
 *
 * @param argc must be 1
 * @param argv (test directory path, M, N, seed)
 * @return int
 */
int exec_expm_test_coo(int argc, char* argv[]) {
    if (argc < 6 || argc > 8) {
        std::cout << "Usage: expm-test-coo matrix-filepath exp-matrix-filepath M N seed [row [col]]"
                  << std::endl;
        return 1;
    }

    std::string filepath(argv[0]);
    std::string solution_path(argv[1]);
    int M = atoi(argv[2]);
    int N = atoi(argv[3]);
    int seed = atoi(argv[4]);
    int row = argc >= 5 ? atoi(argv[5]) : 0;
    int col = argc >= 6 ? atoi(argv[6]) : 0;

    CSRMatrix<float> matrix = CSRMatrix<float>::from_coo(filepath);
    CSRMatrix<float> exp_matrix = CSRMatrix<float>::from_coo(solution_path);

    std::vector<float> v(matrix.columns(), 0);
    v[col] = 1;

    MCME mcme(matrix, v, 1, row, M, N, seed);
    float res = mcme.calculate();
    float sol = exp_matrix.at(row, col);
    float error = std::abs(res - sol) / std::abs(sol);

    std::cout << std::fixed                                  // Precision
              << std::setprecision(2) << error * 100 << "%"  // Error
              << std::setprecision(3) << " (" << res << ")"  // Result
              << std::endl;
    return 0;
}

/**
 * @brief Calculates average accuracy of matrix exponential for a
 *
 * @param argc must be 1
 * @param argv (test directory path, M, N, seed)
 * @return int
 */
int exec_expm_time_coo(int argc, char* argv[]) {
    if (argc < 4 || argc > 6) {
        std::cout << "Usage: expm-time-coo matrix-filepath M N seed [row [col]]" << std::endl;
        return 1;
    }

    std::string filepath(argv[0]);
    int M = atoi(argv[1]);
    int N = atoi(argv[2]);
    int seed = atoi(argv[3]);
    int row = argc >= 4 ? atoi(argv[4]) : 0;
    int col = argc >= 5 ? atoi(argv[5]) : 0;

    CSRMatrix<float> matrix = CSRMatrix<float>::from_coo(filepath);

    std::vector<float> v(matrix.columns(), 0);
    v[col] = 1;

    MCME mcme(matrix, v, 1, row, M, N, seed);

    double time = -omp_get_wtime();
    std::cout << std::fixed << std::setprecision(5)  // Adjust width & precision
              << mcme.calculate() << std::endl;
    time += omp_get_wtime();
    std::cout << "Executed in " << time << " seconds." << std::endl;

    return 0;
}

int exec_expm_b_test(int argc, char* argv[]) {
    if (argc < 6 || argc > 7) {
        std::cout << "Usage: expm-b-test matrix-filepath x-0-filepath b-filepath M N t seed [row]"
                  << std::endl;
        return 1;
    }

    std::string filepath(argv[0]);
    std::string x_0_filepath(argv[1]);
    std::string b_filepath(argv[2]);
    int M = atoi(argv[3]);
    int N = atoi(argv[4]);
    int seed = atoi(argv[5]);
    int row = argc > 6 ? atoi(argv[6]) : 0;

    CSRMatrix<float> matrix = CSRMatrix<float>::from_coo(filepath);

    std::vector<float> x_0 = read_vector(x_0_filepath);
    std::vector<float> b = read_vector(b_filepath);

    MCME mcme(matrix, x_0, 1, row, M, N, seed);

    std::cout << std::fixed << std::setprecision(5)  // Adjust width & precision
              << mcme.calculate_b(b) << std::endl;

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
    } else if (command == "expm-test-coo") {
        return exec_expm_test_coo(argc - 2, &argv[2]);
    } else if (command == "expm-time-coo") {
        return exec_expm_time_coo(argc - 2, &argv[2]);
    } else if (command == "expm-b-test") {
        return exec_expm_b_test(argc - 2, &argv[2]);
    } else {
        exec_help();
        return 1;
    }

    return 0;
}
