#pragma once

#include <optional>
#include <unordered_map>
#include <vector>

#include "time_function.h"
#include "utils/union_find.h"

/**
 * @brief A class representing a SPICE netlist
 *
 */
class Netlist {
   private:
    struct PulseParams {
        float v1, v2, td, tr, tf, pw, per;
    };

    struct Component {
        std::string name;
        int node1, node2;
        float value;
        std::optional<PulseParams> pulse;

        Component(std::string name, int node1, int node2, float value)
            : name(name), node1(node1), node2(node2), value(value), pulse(std::nullopt) {}

        Component(std::string name, int node1, int node2, float value, float v1, float v2, float td,
                  float tr, float tf, float pw, float per)
            : name(name),
              node1(node1),
              node2(node2),
              value(value),
              pulse(PulseParams{v1, v2, td, tr, tf, pw, per}) {}
    };

    std::unordered_map<std::string, int> node_indexes;
    std::vector<std::reference_wrapper<const std::string>> node_names;
    int node_index(const std::string& node);

    std::string title;
    std::vector<Component> vsources;
    std::vector<Component> isources;
    std::vector<Component> resistors;
    std::vector<Component> capacitors;

    float tstep;
    float tstop;

    std::vector<std::string> print;

    void parse_vsource(std::string& line);
    void parse_isource(std::string& line);
    void parse_resistor(std::string& line);
    void parse_capacitor(std::string& line);
    void parse_command(std::string& line);
    
    void handle_zero_voltage_sources(UnionFind& uf);
    void handle_nonzero_voltage_sources(UnionFind& uf);
    void compact_and_update_nodes(UnionFind& uf);

   public:
    /**
     * @brief Construct a new Netlist object from a spice file
     *
     * @param spice_filepath filepath to spice netlist
     */
    Netlist(std::string spice_filepath);

    /**
     * @brief Write the netlist representation to an output stream
     *
     * @param os output stream to write to
     */
    void write(std::ostream& os) const;

    /**
     * @brief Transform the netlist to remove voltage sources
     *
     */
    void remove_voltage_sources();

    /**
     * @brief Transform the netlist to remove voltage sources (v2 implementation)
     *
     */
    void remove_voltage_sources_v2();

    /**
     * @brief Get the vector of voltage sources
     *
     * @return const std::vector<Component>&
     */
    const std::vector<Component>& get_voltage_sources() const { return vsources; }

    /**
     * @brief Get the vector of current sources
     *
     * @return const std::vector<Component>&
     */
    const std::vector<Component>& get_current_sources() const { return isources; }

    /**
     * @brief Get the vector of resistors
     *
     * @return const std::vector<Component>&
     */
    const std::vector<Component>& get_resistors() const { return resistors; }

    /**
     * @brief Get the vector of capacitors
     *
     * @return const std::vector<Component>&
     */
    const std::vector<Component>& get_capacitors() const { return capacitors; }

    /**
     * @brief Get the index from a node name
     *
     * @param name node name
     * @return int node index
     */
    int get_node_index(const std::string& name) const { return node_indexes.at(name); }

    /**
     * @brief Get the node name from a node index
     *
     * @param index node index
     * @return std::string node name
     */
    std::string get_node_name(int index) const { return node_names[index]; }

    /**
     * @brief Get the simulation time step
     *
     * @return float simulation time step
     */
    float get_tstep() const { return tstep; }

    /**
     * @brief Get the simulation stop time
     *
     * @return float simulation stop time
     */
    float get_tstop() const { return tstop; }

    /**
     * @brief Return the number of nodes, including GND
     *
     * @return size_t number of nodes
     */
    size_t nodes_size() const { return node_names.size(); }

    /**
     * @brief Get nodes to print from the .print command
     *
     * @return const std::vector<std::string>& list of nodes to print
     */
    const std::vector<std::string>& get_print_nodes() const { return print; }
};
