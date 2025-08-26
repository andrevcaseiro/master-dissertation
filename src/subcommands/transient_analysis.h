/**
 * @file transient_analysis.h
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief
 * @version 0.1
 * @date 2025-05-28
 *
 * @copyright Copyright (c) 2025
 *
 */

#pragma once

#include <spice/dc_solver.h>
#include <spice/mna.h>
#include <spice/netlist.h>
#include <spice/ode.h>

#include "common.h"
#include "highfm_ode_solver.h"
#include "trapezoidal_ode_solver.h"

/**
 * @brief A struct defining a testing command
 */
struct TransientAnalysis {
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

    void print_initial_conditions(const Eigen::VectorXf& x0, const MNA& mna) const {
        std::cout << "Initial conditions:\n";
        for (Eigen::Index i = 0; i < x0.size(); i++) {
            std::string node_name = mna.get_node_name(i);

            // Format with left-aligned fixed width for node name
            std::cout << std::setw(12) << std::left << node_name + ": " << x0[i] << std::endl;
        }
        std::cout << std::endl;
    }

    void print_solver_params(const std::string& node, int mna_index) const {
        std::cout << "ODE solver: " << solver << std::endl;
        std::cout << "Steps: " << steps << std::endl;
        std::cout << "Time: " << time << std::endl;
        if (solver == "monte-carlo") {
            std::cout << "Samples: " << samples << std::endl;
            std::cout << "Seed: " << seed << std::endl;
        } else if (solver == "trapezoidal") {
            std::cout << "Trapezoidal linear solver: " << method << std::endl;
        }
        std::cout << "DC solver: " << dc_method << std::endl;
        std::cout << "Node: " << node << std::endl;
        std::cout << "MNA index: " << mna_index << std::endl;
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
    size_t print_step = 1000;
    float time = 0;
    long seed = -1;
    std::string solver = "monte-carlo";
    std::string method = "highfm-cg";
    std::string dc_method = "";

    TransientAnalysis(CLI::App& app) {
        auto cmd = app.add_subcommand("tran", "Transient analysis on a circuit.");
        cmd->add_option("filepath", filepath, "Path to the SPICE netlist file")->required();
        cmd->add_flag("-v,--verbose", verbose, "Print detailed information");

        std::vector<std::string> allowed_solvers = {"monte-carlo", "trapezoidal"};
        cmd->add_option("--solver", solver, "Solver method")
            ->check(CLI::IsMember(allowed_solvers))
            ->capture_default_str();

        std::vector<std::string> allowed_methods = {"lu", "cg", "slu", "pardiso", "highfm-cg"};
        cmd->add_option("--method", method, "Linear solver method for trapezoidal integration")
            ->check(CLI::IsMember(allowed_methods))
            ->capture_default_str();

        std::vector<std::string> allowed_dc_methods = {"lu", "cg", "slu", "pardiso"};
        cmd->add_option("--dc-method", dc_method,
                        "Linear solver method for DC analysis (defaults to method)")
            ->check(CLI::IsMember(allowed_dc_methods))
            ->capture_default_str();

        cmd->add_option("-M,--samples", samples, "Number of Monte Carlo samples")
            ->capture_default_str();
        cmd->add_option("-N,--steps", steps, "Number of time steps")->capture_default_str();
        cmd->add_option("-t,--time", time, "Final time")->capture_default_str();
        cmd->add_option("-s,--seed", seed, "Random seed (-1 for random)")->capture_default_str();
        cmd->add_option("-p,--print-step", print_step, "Number of output steps")
            ->capture_default_str();
        cmd->callback([this]() { execute(); });
    }

    void execute() {
        double start_time;

        // Override dc_method if not specified
        if (dc_method.empty()) {
            dc_method = method;
        }

        // Load and process netlist
        start_time = omp_get_wtime();
        Netlist netlist(filepath);
        if (verbose) print_netlist(netlist, true);
        netlist.remove_voltage_sources_v2();
        if (verbose) print_netlist(netlist, false);
        print_execution_time("Netlist processing", start_time);

        // Get the default tstop and tstep from netlist
        if (time == 0) time = netlist.get_tstop();
        if (steps == 0) steps = time / netlist.get_tstep();

        // Create MNA system
        start_time = omp_get_wtime();
        MNA mna(netlist);
        print_execution_time("MNA creation", start_time);
        if (verbose) print_mna_matrices(mna);

        // Get the node to analyze from netlist print list
        const auto& print_nodes = netlist.get_print_nodes();
        if (print_nodes.empty()) {
            throw std::runtime_error("No print nodes specified in netlist");
        }
        int mna_index = mna.get_mna_index(print_nodes[0]);
        if (mna_index < 0) {
            throw std::runtime_error("Cannot analyze ground node");
        }

        // Create ODE system
        start_time = omp_get_wtime();
        ODE ode(mna);
        print_execution_time("ODE creation", start_time);
        if (verbose) print_ode_matrices(ode);

        // Calculate initial conditions using DC analysis
        start_time = omp_get_wtime();
        DCSolver dc(mna);

        // Select DC solver method
        DCSolver::Method dc_solver_method;
        if (dc_method == "lu")
            dc_solver_method = DCSolver::Method::LU;
        else if (dc_method == "cg")
            dc_solver_method = DCSolver::Method::CG;
        else if (dc_method == "pardiso")
            dc_solver_method = DCSolver::Method::PARDISO;
        else
            dc_solver_method = DCSolver::Method::SLU;

        Eigen::VectorXf x0 = dc.solve(dc_solver_method);
        print_execution_time("DC analysis", start_time);
        if (verbose) print_initial_conditions(x0, mna);

        if (verbose) print_solver_params(print_nodes[0], mna_index);

        // Convert initial conditions to vector format
        std::vector<float> x0_vec(x0.data(), x0.data() + x0.size());

        std::vector<float> result;
        if (solver == "monte-carlo") {
            // Convert to CSR format for Monte Carlo solver
            start_time = omp_get_wtime();
            CSRMatrix<float> A_csr(ode.A());
            print_execution_time("CSR conversion", start_time);

            // Solve ODE using Monte Carlo method
            start_time = omp_get_wtime();
            MonteCarloODESolver mc_solver(A_csr, ode.b(), x0_vec, time, mna_index, samples, steps,
                                          seed);
            result = mc_solver.solve_sequence(print_step);
        } else if (solver == "trapezoidal") {
            start_time = omp_get_wtime();

            if (method == "pardiso") {
                // Solve ODE using HighFM with PARDISO method
                HigFMODESolver highfm_solver(ode.A(), ode.b(), x0_vec, time, mna_index, steps);
                result = highfm_solver.solve_sequence();
            } else if (method == "highfm-cg") {
                // Solve ODE using HighFM with CG method
                HigFMODESolver highfm_solver(ode.A(), ode.b(), x0_vec, time, mna_index, steps);
                result = highfm_solver.solve_sequence_cg();
            } else {
                // Solve ODE using traditional Trapezoidal method
                TrapezoidalODESolver trap_solver(ode.A(), ode.b(), x0_vec, time, mna_index, steps);

                // Select solver method
                TrapezoidalODESolver::Method solve_method;
                if (method == "lu")
                    solve_method = TrapezoidalODESolver::Method::LU;
                else if (method == "cg")
                    solve_method = TrapezoidalODESolver::Method::CG;
                else
                    solve_method = TrapezoidalODESolver::Method::SUPERLU_MT;
                result = trap_solver.solve_sequence(solve_method);
            }
        }
        print_execution_time("ODE solver", start_time);

        print_results(result);
    }
};
