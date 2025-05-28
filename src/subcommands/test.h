/**
 * @file test.h
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief
 * @version 0.1
 * @date 2025-05-28
 *
 * @copyright Copyright (c) 2025
 *
 */

#pragma once

#include <monte_carlo_ode_solver.h>
#include <omp.h>
#include <spice/dc_solver.h>
#include <spice/mna.h>
#include <spice/netlist.h>
#include <spice/ode.h>

#include <iomanip>

#include "common.h"

/**
 * @brief A struct defining a testing command
 */
struct Test {
   private:
    void print_netlist(const Netlist& netlist, bool original = true) const {
        std::cout << (original ? "Original" : "Reduced") << " netlist:" << std::endl;
        netlist.write(std::cout);
        std::cout << std::endl;
    }

    void print_mna_matrices(const MNA& mna) const {
        if (mna.size() <= 100) {
            std::cout << "MNA matrices:\n";
            std::cout << "\nConductance matrix G:\n" << Eigen::MatrixXf(mna.get_G()) << std::endl;
            std::cout << "\nCapacitance matrix C:\n" << mna.get_C().diagonal() << std::endl;

            std::cout << "\nMNA sources vector b:\n";
            const auto& mna_b = mna.get_b();
            for (size_t i = 0; i < mna_b.size(); i++) {
                if (mna_b[i]) {
                    std::cout << i << ": " << mna_b[i]->to_string() << std::endl;
                }
            }
            std::cout << std::endl;
        } else {
            std::cout << "MNA matrices size: " << mna.size() << " (too large to display)"
                      << std::endl;
        }
    }

    void print_ode_matrices(const ODE& ode) const {
        if (ode.size() <= 100) {
            std::cout << "ODE matrices:";
            std::cout << "\nA matrix:\n" << Eigen::MatrixXf(ode.A()) << std::endl;

            std::cout << "\nODE sources vector b:\n";
            const auto& b = ode.b();
            for (size_t i = 0; i < b.size(); i++) {
                if (b[i]) {
                    std::cout << i << ": " << b[i]->to_string() << std::endl;
                }
            }
            std::cout << std::endl;
        } else {
            std::cout << "\nODE matrices size: " << ode.size() << " (too large to display)"
                      << std::endl;
        }
    }

    void print_initial_conditions(const Eigen::VectorXf& x0) const {
        std::cout << "Initial conditions:\n";
        for (Eigen::Index i = 0; i < x0.size(); i++) {
            std::cout << i << ": " << x0[i] << std::endl;
        }
        std::cout << std::endl;
    }

    void print_solver_params() const {
        std::cout << "Solving ODE with Monte Carlo method:" << std::endl;
        std::cout << "Samples: " << samples << std::endl;
        std::cout << "Steps: " << steps << std::endl;
        std::cout << "Time: " << time << std::endl;
        std::cout << "Seed: " << seed << std::endl;
        std::cout << std::endl;
    }

    void print_results(const std::vector<float>& result) const {
        float delta_t = time / (result.size() - 1);

        std::cout << "Results:" << std::endl;
        for (size_t i = 0; i < result.size(); i++) {
            std::cout << std::scientific << std::setprecision(6) << i * delta_t << " "
                      << std::scientific << std::setprecision(7) << result[i] << std::endl;
        }
        std::cout << std::endl;
    }

    void print_execution_time(const std::string& step, double start_time) const {
        double end_time = omp_get_wtime();
        std::cout << step << " time: " << std::fixed << std::setprecision(4)
                  << (end_time - start_time) << " s" << std::endl
                  << std::endl;
    }

   public:
    std::string filepath;
    bool verbose = false;
    size_t samples = 1000;
    size_t steps = 0;
    float time = 0;
    long seed = -1;

    Test(CLI::App& app) {
        auto cmd = app.add_subcommand("test", "This is a development command.");
        cmd->add_option("filepath", filepath, "Path to the SPICE netlist file")->required();
        cmd->add_flag("-v,--verbose", verbose, "Print detailed information");
        cmd->add_option("-M,--samples", samples, "Number of Monte Carlo samples")
            ->capture_default_str();
        cmd->add_option("-N,--steps", steps, "Number of time steps")->capture_default_str();
        cmd->add_option("-t,--time", time, "Final time")->capture_default_str();
        cmd->add_option("-s,--seed", seed, "Random seed (-1 for random)")->capture_default_str();
        cmd->callback([this]() { execute(); });
    }

    void execute() {
        double start_time;

        // Load and process netlist
        start_time = omp_get_wtime();
        Netlist netlist(filepath);
        if (verbose) print_netlist(netlist, true);
        netlist.remove_voltage_sources();
        if (verbose) print_netlist(netlist, false);
        print_execution_time("Netlist processing", start_time);

        // Create MNA system
        start_time = omp_get_wtime();
        MNA mna(netlist);
        if (verbose) print_mna_matrices(mna);
        print_execution_time("MNA creation", start_time);

        // Create ODE system
        start_time = omp_get_wtime();
        ODE ode(mna);
        if (verbose) print_ode_matrices(ode);
        print_execution_time("ODE creation", start_time);

        // Calculate initial conditions using DC analysis
        start_time = omp_get_wtime();
        DCSolver dc(mna);
        Eigen::VectorXf x0 = dc.solve();
        std::vector<float> x0_vec(x0.data(), x0.data() + x0.size());
        if (verbose) print_initial_conditions(x0);
        print_execution_time("DC analysis", start_time);

        // Convert to CSR format and set up solver parameters
        start_time = omp_get_wtime();
        CSRMatrix<float> A_csr(ode.A());
        print_execution_time("CSR conversion", start_time);

        if (time == 0) time = netlist.get_tstop();
        if (steps == 0) steps = time / netlist.get_tstep();

        if (verbose) print_solver_params();

        // Solve ODE using Monte Carlo method
        start_time = omp_get_wtime();
        MonteCarloODESolver solver(A_csr, ode.b(), x0_vec, time, 0, samples, steps, seed);
        auto result = solver.solve_sequence();
        print_execution_time("Monte Carlo simulation", start_time);

        print_results(result);
    }
};
