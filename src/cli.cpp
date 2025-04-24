/**
 * @file cli.cpp
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief CLI interface for the Monte Carlo ODE solver
 * @version 0.1
 * @date 2025-03-24
 *
 * @copyright Copyright (c) 2025
 *
 */

#include <omp.h>

#include <CLI11.hpp>
#include <Eigen/Sparse>

#include "matrix/csrd_matrix.h"
#include "monte_carlo_ode_solver.h"
#include "spice/spice_parser.h"
#include "utils/read_vector.h"

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

/**
 * @brief A struct defining the general case ODE solver command
 */
struct Solve {
    std::string A_filename;
    std::string b_filename;
    std::string x_0_filename;
    std::string res_filename;
    size_t M;
    size_t N;
    size_t row;
    float t;
    long seed;

    Solve(CLI::App& app)
        : b_filename(""), x_0_filename(""), res_filename(""), row(0), t(1), seed(-1) {
        auto cmd = app.add_subcommand("solve", "Solve x'=Ax+b with x(0)=x_0 at time T.");

        cmd->add_option("A_file", A_filename, "Filename of matrix A")->required();
        cmd->add_option("-b,--b_file", b_filename, "Filename of vector b");
        cmd->add_option("--x0,--x0_file", x_0_filename, "Filename of vector x_0");
        cmd->add_option("-r,--res_file", res_filename, "Filename of vector with expected result");

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

        float expected_res = 0;
        if (!res_filename.empty()) {
            expected_res = read_vector(res_filename)[row];
        }

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

/**
 * @brief A struct defining the general case ODE solver command
 */
struct SolveSpice {
    std::string spice_filename;
    std::string res_filename;
    size_t M;
    size_t N;
    std::string node;
    float t;
    long seed;
    bool solve_sequence;
    bool verbose;

    SolveSpice(CLI::App& app)
        : res_filename(""), N(0), node("-"), t(0), seed(-1), solve_sequence(false), verbose(false) {
        auto cmd = app.add_subcommand("solve-spice", "Solve transient analysis from a spice file.");

        cmd->add_option("spice_file", spice_filename, "Filename")->required();
        cmd->add_option("-r,--res_file", res_filename, "Filename of vector with expected result");

        cmd->add_option("M", M, "Number of samples")->required();
        cmd->add_option("N", N, "Splitting parameter (0 to read from file)")->capture_default_str();
        cmd->add_option("node", node, "Solution node (- to read from file)")->capture_default_str();
        cmd->add_option("t", t, "Solution time (-1 to read from  file)")->capture_default_str();
        cmd->add_option("seed", seed, "Seed (-1 for random)")->capture_default_str();

        cmd->add_flag("-s", solve_sequence, "Print the full sequence of values")
            ->capture_default_str();
        cmd->add_flag("-v,--verbose", verbose, "Print matrices and vectors")->capture_default_str();

        cmd->callback([this]() { execute(); });
    }

    void execute() {
        SpiceParser sp(spice_filename);
        auto& C = sp.get_C();
        CSRMatrix<float>& G = sp.get_G();
        auto& b = sp.get_b();

        if (t == 0) t = sp.get_time();
        if (N == 0) N = sp.get_time() / sp.get_timestep();
        if (node == "-") node = sp.get_print_node();

        size_t row = node.empty() ? 0 : sp.get_node_index(node);

        /* Print initial matrices */
        if (verbose) {
            std::cout << "time: " << t << std::endl;
            std::cout << "N: " << N << std::endl;
            std::cout << "node: " << node << std::endl;

            std::cout << "C:" << std::endl;
            for (auto it : C) {
                std::cout << it << std::endl;
            }

            std::cout << std::endl << "G:" << std::endl;
            G.print();

            std::cout << std::endl << "Bu:" << std::endl;
            for (const auto& it : b) {
                std::cout << it->to_string() << std::endl;
            }
        }

        std::vector<Eigen::Triplet<float>> triplets;
        triplets.reserve(G.nnz()); // If you know or can estimate it
        
        for (size_t row = 0; row < G.rows(); ++row) {
            for (auto entry : G.row(row)) {
                triplets.emplace_back(entry.row(), entry.column(), entry.value());
            }
        }
        
        Eigen::SparseMatrix<float> eigen_G(G.rows(), G.columns());
        eigen_G.setFromTriplets(triplets.begin(), triplets.end());        

        double initial_time = -omp_get_wtime();

        Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower | Eigen::Upper>
            eigen_solver;
        eigen_solver.compute(eigen_G);
        if (eigen_solver.info() != Eigen::Success) {
            std::cout << "Failed to compute the initial conditions" << std::endl;
            return;
        }

        Eigen::VectorXf eigen_b(b.size());
        for (size_t i = 0; i < b.size(); ++i) {
            eigen_b(i) = (*b[i])(0);
        }

        Eigen::VectorXf x = eigen_solver.solve(eigen_b);
        if (eigen_solver.info() != Eigen::Success) {
            std::cout << "Failed to compute the initial conditions" << std::endl;
            return;
        }

        std::vector<float> x_0(C.size(), 0);
        for (size_t i = 0; i < b.size(); ++i) {
            x_0[i] = x[i];
        }

        initial_time += omp_get_wtime();
        if (!solve_sequence) {
            std::cout << std::setw(20) << std::left << "DC analysis time:" << std::setw(10)
                      << std::right << std::fixed << std::setprecision(3) << initial_time << " s"
                      << std::endl;
        }

        if (verbose) {
            std::cout << std::endl << "Initial Transient Solution:" << std::endl;
            for (auto v : x_0) {
                std::cout << std::defaultfloat;
                std::cout << v << std::endl;
            }
        }

        /* Compute A = -C^-1 G and b = C^-1 Bu */
        for (size_t i = 0; i < C.size(); ++i) {
            float c = C[i];
            c = c != 0 ? 1 / c : 1;

            for (auto entry : G.row(i)) {
                entry.value() *= -c;
            }

            *b[i] *= c;
        }

        /* Print final matrices */
        if (verbose) {
            std::cout << std::endl << "A:" << std::endl;
            G.print();

            std::cout << std::endl << "Bu:" << std::endl;
            for (const auto& it : b) {
                std::cout << it->to_string() << std::endl;
            }
        }

        MonteCarloODESolver solver(G, b, x_0, t, row, M, N, seed);
        if (solve_sequence) {
            std::vector<float> res = solver.solve_sequence();

            float delta_t = 2 * t / N;

            std::cout << std::setw(10) << std::setprecision(5);
            for (size_t i = 0; i < res.size(); i++) {
                std::cout << delta_t * i << " " << res[i] << std::endl;
            }
        } else {
            double time = -omp_get_wtime();
            float res = solver.solve();
            time += omp_get_wtime();

            std::cout << std::setw(20) << std::left << "Result:" << std::setw(10) << std::right
                      << std::setprecision(7) << res << std::endl;
            std::cout << std::setw(20) << std::left << "Time:" << std::setw(10) << std::right
                      << std::setprecision(7) << time << " s" << std::endl;
        }
    }
};

/**
 * @brief Entry point, defines parses CLI arguments and executes command
 *
 * @param argc number of arguments
 * @param argv list of arguments
 * @return int exit status
 */
int main(int argc, char** argv) {
    CLI::App app{"A cli tool to solve systems of DAEs."};

    Solve solveCmd{app};
    SolveSequence solveSequenceCmd{app};
    ExpM expmCmd{app};
    SolveSpice solveSpiceCmd{app};

    CLI11_PARSE(app, argc, argv);
    return 0;
}
