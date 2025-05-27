#pragma once

#include <monte_carlo_ode_solver.h>
#include <spice/mna.h>
#include <spice/netlist.h>
#include <spice/ode.h>

#include "common.h"

/**
 * @brief A struct defining a testing command
 */
struct Test {
    bool verbose = false;
    size_t samples = 1000;
    size_t steps = 0;
    float time = 0;
    long seed = -1;

    Test(CLI::App& app) {
        auto cmd = app.add_subcommand("test", "This is a development command.");
        cmd->add_flag("-v,--verbose", verbose, "Print detailed information");
        cmd->add_option("-M,--samples", samples, "Number of Monte Carlo samples")
            ->capture_default_str();
        cmd->add_option("-N,--steps", steps, "Number of time steps")->capture_default_str();
        cmd->add_option("-t,--time", time, "Final time")->capture_default_str();
        cmd->add_option("-s,--seed", seed, "Random seed (-1 for random)")->capture_default_str();
        cmd->callback([this]() { execute(); });
    }

    void execute() {
        Netlist netlist("./test/spice/line1.spice");
        if (verbose) {
            std::cout << "Original netlist:" << std::endl;
            netlist.save(std::cout);
        }
        netlist.remove_voltage_sources();
        if (verbose) {
            std::cout << "Reduced netlist:" << std::endl;
            netlist.save(std::cout);
        }

        MNA mna(netlist);
        if (verbose) {
            if (mna.size() <= 100) {
                std::cout << "\nMNA matrices:\n";
                std::cout << "\nConductance matrix G:\n"
                          << Eigen::MatrixXf(mna.get_G()) << std::endl;
                std::cout << "\nCapacitance matrix C:\n" << mna.get_C().diagonal() << std::endl;

                std::cout << "\nMNA sources vector b:\n";
                const auto& mna_b = mna.get_b();
                for (size_t i = 0; i < mna_b.size(); i++) {
                    if (mna_b[i]) {
                        std::cout << "b[" << i << "] = " << mna_b[i]->to_string() << std::endl;
                    }
                }
            } else {
                std::cout << "\nMNA matrices size: " << mna.size() << " (too large to display)"
                          << std::endl;
            }
        }

        ODE ode(mna);
        if (verbose) {
            if (ode.size() <= 100) {
                std::cout << "\nODE matrices:\n";
                std::cout << "\nA matrix:\n" << Eigen::MatrixXf(ode.A()) << std::endl;

                std::cout << "\nODE sources vector b:\n";
                const auto& b = ode.b();
                for (size_t i = 0; i < b.size(); i++) {
                    if (b[i]) {
                        std::cout << "b[" << i << "] = " << b[i]->to_string() << std::endl;
                    }
                }
            } else {
                std::cout << "\nODE matrices size: " << ode.size() << " (too large to display)"
                          << std::endl;
            }
        }

        // Convert sparse matrix to CSR format
        CSRMatrix<float> A_csr(ode.A());
        std::vector<float> x0(ode.size(), 0.0f);
        x0[0] = 1.0f;  // Set initial condition

        if (time == 0) time = netlist.get_tstop();
        if (steps == 0) steps = time / netlist.get_tstep();

        if (verbose) {
            std::cout << "\nSolving ODE with Monte Carlo method:" << std::endl;
            std::cout << "Samples: " << samples << std::endl;
            std::cout << "Steps: " << steps << std::endl;
            std::cout << "Time: " << time << std::endl;
            std::cout << "Seed: " << seed << std::endl;
        }

        A_csr.print();
        A_csr.print_csr();

        MonteCarloODESolver solver(A_csr, ode.b(), x0, time, 0, samples, steps, seed);
        auto result = solver.solve_sequence();

        std::cout << "\nResults:" << std::endl;
        for (size_t i = 0; i < result.size(); i++) {
            std::cout << (time * i / (result.size() - 1)) << " " << result[i] << std::endl;
        }
    }
};
