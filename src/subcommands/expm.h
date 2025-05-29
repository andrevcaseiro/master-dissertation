/**
 * @file expm.h
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief Command for calculating matrix exponentials
 * @version 0.1
 * @date 2025-05-29
 *
 * @copyright Copyright (c) 2025
 *
 */

#pragma once

#include "common.h"

/**
 * @brief A struct defining the matrix exponential command
 */
struct ExpM {
    std::string A_filename;
    std::string res_filename;
    size_t M;
    size_t N;
    size_t row;
    size_t col;
    long seed;

    ExpM(CLI::App& app) : res_filename(""), row(0), col(0), seed(-1) {
        auto cmd = app.add_subcommand("exp", "Calculate exp(A).");

        cmd->add_option("A_file", A_filename, "Filename of matrix A")->required();
        cmd->add_option("-r,--res_file", res_filename, "Filename of matrix exp(A)");

        cmd->add_option("M", M, "Number of samples")->required();
        cmd->add_option("N", N, "Splitting parameter")->required();
        cmd->add_option("row", row, "Solution row")->capture_default_str();
        cmd->add_option("col", col, "Solution column")->capture_default_str();
        cmd->add_option("seed", seed, "Seed (-1 for random)")->capture_default_str();

        cmd->callback([this]() { execute(); });
    }

    void execute() {
        CSRMatrix<float> A = CSRMatrix<float>::from_coo(A_filename);

        float expected_res = 0;
        if (!res_filename.empty()) {
            expected_res = CSRMatrix<float>::from_coo(res_filename).at(row, col);
        }

        std::vector<float> x_0(A.rows(), 0);
        x_0[col] = 1;
        std::vector<float> b(A.rows(), 0);
        float t = 1;

        MonteCarloODESolver solver(A, b, x_0, t, row, M, N, seed);

        double time = -omp_get_wtime();
        float res = solver.solve();
        time += omp_get_wtime();

        std::cout << std::setw(8) << std::left << "Result:" << std::setw(10) << std::right
                  << std::fixed << std::setprecision(5) << res << std::endl;
        std::cout << std::setw(8) << std::left << "Time:" << std::setw(10) << std::right
                  << std::fixed << std::setprecision(3) << time << " s" << std::endl;

        if (!res_filename.empty()) {
            float error = std::abs(res - expected_res) / std::abs(expected_res);

            std::cout << std::setw(8) << std::left << "Error:" << std::setw(10) << std::right
                      << std::fixed << std::setprecision(3) << error * 100 << " %" << std::endl;
        }
    }
};
