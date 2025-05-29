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

#include <CLI11.hpp>

#include "subcommands/expm.h"
#include "subcommands/solve.h"
#include "subcommands/solve_sequence.h"
#include "subcommands/transient_analysis.h"

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
    TransientAnalysis tran{app};

    app.require_subcommand(1);  // Require exactly one subcommand

    CLI11_PARSE(app, argc, argv);
    return 0;
}
