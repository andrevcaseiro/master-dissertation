#include "mna.h"

#include <unordered_map>
#include <unordered_set>

// Local hash function for pair of size_t
struct PairHash {
    size_t operator()(const std::pair<size_t, size_t>& p) const {
        return std::hash<size_t>()(p.first) ^ (std::hash<size_t>()(p.second) << 1);
    }
};

int MNA::get_mna_index(const std::string& node_name) const {
    if (node_name == "0" || node_name == "gnd" || node_name == "GND") {
        return -1;  // Ground node
    }
    return _netlist.get_node_index(node_name) - 1;  // Convert from netlist index to MNA index
}

std::string MNA::get_node_name(int mna_index) const {
    // Convert from MNA index to netlist index by adding 1
    int netlist_index = mna_index + 1;
    try {
        return _netlist.get_node_name(netlist_index);
    } catch (const std::exception&) {
        // Fallback if node name not found
        return "node_" + std::to_string(mna_index);
    }
}

MNA::MNA(const Netlist& netlist)
    : _netlist(netlist), _size(netlist.nodes_size() - 1), G(_size, _size), C(_size), b(_size) {
    // Initialize b vector with zeros
    for (size_t i = 0; i < _size; ++i) {
        b[i] = std::make_unique<ConstantFunction>(0.0f);
    }

    // Create a map to accumulate conductance values at each position
    std::unordered_map<std::pair<size_t, size_t>, float, PairHash> G_entries;

    // Process resistors to build conductance matrix
    for (const auto& r : netlist.get_resistors()) {
        float g = 1.0f / r.value;  // Conductance = 1/R
        auto n1 = r.node1 - 1;
        auto n2 = r.node2 - 1;

        if (n1 >= 0) {
            // Connected to ground at node2
            G_entries[{n1, n1}] += g;
        }
        if (n2 >= 0) {
            // Connected to ground at node2
            G_entries[{n2, n2}] += g;
        }
        if (n1 >= 0 && n2 >= 0) {
            G_entries[{n1, n2}] -= g;
            G_entries[{n2, n1}] -= g;
        }
    }

    // Convert map entries to triplets and build sparse matrix
    std::vector<Eigen::Triplet<float>> triplets;
    triplets.reserve(G_entries.size());

    for (const auto& entry : G_entries) {
        triplets.emplace_back(entry.first.first, entry.first.second, entry.second);
    }

    G.setFromTriplets(triplets.begin(), triplets.end());

    // Process current sources for RHS vector
    for (const auto& is : netlist.get_current_sources()) {
        auto n1 = is.node1 - 1;
        auto n2 = is.node2 - 1;

        if (is.pulse) {
            const auto& p = is.pulse.value();
            if (n2 >= 0) {
                PulseFunction pf(p.v1, p.v2, p.td, p.tr, p.tf, p.pw, p.per);
                b[n2] = *b[n2] + pf;
            }
            if (n1 >= 0) {
                PulseFunction pf(-p.v1, -p.v2, p.td, p.tr, p.tf, p.pw, p.per);
                b[n1] = *b[n1] + pf;
            }
        } else {
            if (n2 >= 0) *b[n2] += is.value;
            if (n1 >= 0) *b[n1] += -is.value;
        }
    }

    // Process voltage sources
    for (const auto& vs : netlist.get_voltage_sources()) {
        throw std::runtime_error("Voltage sources are not supported: " + vs.name);
    }

    // Process capacitors
    for (const auto& c : netlist.get_capacitors()) {
        if (c.node1 == 0) {
            C.diagonal()[c.node2 - 1] = c.value;
        } else if (c.node2 == 0) {
            C.diagonal()[c.node1 - 1] = c.value;
        } else {
            throw std::runtime_error("Capacitors are only supported when connected to GND.");
        }
    }
}
