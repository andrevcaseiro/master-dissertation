#include <omp.h>

#include <CLI11.hpp>

#include "matrix/csrd_matrix.h"
#include "monte_carlo_ode_solver.h"
#include "read_vector.h"

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
        cmd->add_option("seed", seed, "Seed")->capture_default_str();

        cmd->callback([this]() { execute(); });
    }

    void execute() {
        CSRMatrix<float> A = CSRMatrix<float>::from_coo(A_filename);

        float expected_res;
        if (!res_filename.empty()) {
            expected_res = CSRMatrix<float>::from_coo(res_filename).at(row, col);
        }

        std::vector<float> x_0(A.rows(), 0);
        x_0[col] = 1;
        std::vector<float> b(A.rows(), 0);
        float t = 1;

        MonteCarloODESolver solver(A, x_0, b, t, row, M, N, seed);

        double time = -omp_get_wtime();
        float res = solver.solve();
        time += omp_get_wtime();

        if (!res_filename.empty()) {
            float error = std::abs(res - expected_res) / std::abs(expected_res);

            std::cout << std::fixed << std::setprecision(5) << res    // Result
                      << ", " << std::setprecision(3) << error * 100  // Error
                      << std::endl;
        } else {
            std::cout << std::fixed << std::setprecision(5) << res  // Result
                      << std::endl;
        }
    }
};

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
        cmd->add_option("seed", seed, "Seed")->capture_default_str();

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

        float expected_res;
        if (!res_filename.empty()) {
            expected_res = read_vector(res_filename)[row];
        }

        MonteCarloODESolver solver(A, x_0, b, t, row, M, N, seed);

        double time = -omp_get_wtime();
        float res = solver.solve();
        time += omp_get_wtime();

        if (!res_filename.empty()) {
            float error = std::abs(res - expected_res) / std::abs(expected_res);

            std::cout << std::fixed << std::setprecision(5) << res    // Result
                      << ", " << std::setprecision(3) << error * 100  // Error
                      << std::endl;
        } else {
            std::cout << std::fixed << std::setprecision(5) << res  // Result
                      << std::endl;
        }
    }
};

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
        cmd->add_option("seed", seed, "Seed")->capture_default_str();

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

        MonteCarloODESolver solver(A, x_0, b, t, row, M, N, seed);

        std::vector<float> res = solver.solve_sequence();

        std::cout << std::fixed << std::setprecision(5);  // Adjust width & precision
        for (size_t i = 0; i < res.size(); i++) {
            std::cout << res[i] << std::endl;
        }
    }
};

int main(int argc, char** argv) {
    CLI::App app{"A cli tool to solve systems of DAEs."};

    Solve solveCmd{app};
    SolveSequence solveSequence{app};
    ExpM expm{app};

    CLI11_PARSE(app, argc, argv);
    return 0;
}
