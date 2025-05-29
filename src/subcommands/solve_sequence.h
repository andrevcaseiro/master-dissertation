/**
 * @file solve_sequence.h
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief Command for solving ODEs across a time sequence
 * @version 0.1
 * @date 2025-05-29
 *
 * @copyright Copyright (c) 2025
 *
 */

#pragma once

#include <utils/read_vector.h>

#include "common.h"

/**
 * @brief A struct defining the command to solve an ODE with intermidiary results
 */
struct SolveSequence {
    std::string A_filename;
    std::string b_filename;
    std::string x_0_filename;
    std::string res_filename;
    size_t M;
    size_t N;
    size_t row;
    float t;
    long seed;

    SolveSequence(CLI::App& app) : b_filename(""), x_0_filename(""), row(0), t(1), seed(-1) {
        auto cmd =
            app.add_subcommand("solve-sequence", "Solve x'=Ax+b with x(0)=x_0 at times 0 to T.");

        cmd->add_option("A_filename", A_filename, "Filename of matrix A")->required();
        cmd->add_option("-b,--b_file", b_filename, "Filename of vector b");
        cmd->add_option("--x0,--x0_file", x_0_filename, "Filename of vector x_0");

        cmd->add_option("M", M, "Number of samples")->required();
        cmd->add_option("N", N, "Splitting parameter")->required();
        cmd->add_option("row", row, "Solution row")->capture_default_str();
        cmd->add_option("t", t, "Solution time")->capture_default_str();
        cmd->add_option("seed", seed, "Seed (-1 for random)")->capture_default_str();

        cmd->callback([this]() { execute(); });
    }

    void execute() {
        CSRMatrix<float> A = CSRMatrix<float>::from_coo(A_filename);

        std::vector<float> x_0;
        if (!x_0_filename.empty()) {
            x_0 = read_vector(x_0_filename);
        } else {
            x_0 = std::vector<float>(A.rows(), 0);
        }

        std::vector<float> b;
        if (!b_filename.empty()) {
            b = read_vector(b_filename);
        } else {
            b = std::vector<float>(A.rows(), 0);
        }

        MonteCarloODESolver solver(A, b, x_0, t, row, M, N, seed);

        std::vector<float> res = solver.solve_sequence();

        std::cout << std::setw(10) << std::fixed << std::setprecision(5);
        for (size_t i = 0; i < res.size(); i++) {
            std::cout << res[i] << std::endl;
        }
    }
};
