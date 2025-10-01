#include "netlist.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>

#include "utils/union_find.h"

static void to_lowercase(std::string& s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
}

static std::vector<std::string_view> tokenize(const std::string& line) {
    /* Spice uses spaces, equal sign, comma or parenthesis as delimiters */
    static const std::regex delim(R"([\s=,()]+)");
    std::vector<std::string_view> tokens;

    for (std::sregex_token_iterator it(line.begin(), line.end(), delim, -1), end; it != end; ++it) {
        const auto& match = *it;
        if (!match.matched || match.length() == 0) continue;
        tokens.emplace_back(&*match.first, match.length());
    }

    return tokens;
}

int Netlist::node_index(const std::string& node) {
    if (node_indexes.find(node) == node_indexes.end()) {
        node_indexes[node] = node_indexes.size();  // Assign a new index
    }
    return node_indexes[node];
}

void Netlist::parse_isource(std::string& line) {
    auto tokens = tokenize(line);

    std::string name = std::string(tokens[0]);
    int node1 = node_index(std::string(tokens[1]));
    int node2 = node_index(std::string(tokens[2]));
    float value = std::stof(std::string(tokens[3]));

    if (tokens.size() > 4 && tokens[4] == "pulse") {
        float v1 = std::stof(std::string(tokens[5]));
        float v2 = std::stof(std::string(tokens[6]));
        float td = std::stof(std::string(tokens[7]));
        float tr = std::stof(std::string(tokens[8]));
        float tf = std::stof(std::string(tokens[9]));
        float pw = std::stof(std::string(tokens[10]));
        float per = std::stof(std::string(tokens[11]));

        isources.emplace_back(name, node1, node2, value, v1, v2, td, tr, tf, pw, per);
    } else {
        isources.emplace_back(name, node1, node2, value);
    }
}

void Netlist::parse_vsource(std::string& line) {
    auto tokens = tokenize(line);

    std::string name = std::string(tokens[0]);
    int node1 = node_index(std::string(tokens[1]));
    int node2 = node_index(std::string(tokens[2]));
    float value = std::stof(std::string(tokens[3]));

    if (tokens.size() > 4 && tokens[4] == "pulse") {
        float v1 = std::stof(std::string(tokens[5]));
        float v2 = std::stof(std::string(tokens[6]));
        float td = std::stof(std::string(tokens[7]));
        float tr = std::stof(std::string(tokens[8]));
        float tf = std::stof(std::string(tokens[9]));
        float pw = std::stof(std::string(tokens[10]));
        float per = std::stof(std::string(tokens[11]));

        vsources.emplace_back(name, node1, node2, value, v1, v2, td, tr, tf, pw, per);
    } else {
        vsources.emplace_back(name, node1, node2, value);
    }
}

void Netlist::parse_resistor(std::string& line) {
    auto tokens = tokenize(line);

    std::string name = std::string(tokens[0]);
    int node1 = node_index(std::string(tokens[1]));
    int node2 = node_index(std::string(tokens[2]));
    float value = std::stof(std::string(tokens[3]));

    resistors.emplace_back(name, node1, node2, value);
}

void Netlist::parse_capacitor(std::string& line) {
    auto tokens = tokenize(line);

    std::string name = std::string(tokens[0]);
    int node1 = node_index(std::string(tokens[1]));
    int node2 = node_index(std::string(tokens[2]));
    float value = std::stof(std::string(tokens[3]));

    capacitors.emplace_back(name, node1, node2, value);
}

void Netlist::parse_inductor(std::string& line) {
    auto tokens = tokenize(line);

    std::string name = std::string(tokens[0]);
    int node1 = node_index(std::string(tokens[1]));
    int node2 = node_index(std::string(tokens[2]));
    std::stof(std::string(tokens[3]));

    // Static variable to track if warning has been shown
    static bool warning_shown = false;
    if (!warning_shown) {
        std::cerr << "warning: inductors are treated as short circuits (0V voltage sources)"
                  << std::endl;
        warning_shown = true;
    }

    // Treat inductor as a 0V voltage source
    vsources.emplace_back("v_" + name, node1, node2, 0.0f);
}

void Netlist::parse_command(std::string& line) {
    auto tokens = tokenize(line);

    if (tokens[0] == ".tran") {
        tstep = std::stof(std::string(tokens[1]));
        tstop = std::stof(std::string(tokens[2]));
    } else if (tokens[0] == ".print") {
        if (tokens[1] == "tran") {
            for (size_t i = 2; i + 1 < tokens.size(); i += 2) {
                if (tokens[i] == "v") {
                    print.emplace_back(tokens[i + 1]);
                }
            }
        } else {
            std::cerr << "warning: the following print will be ignored: " << std::endl
                      << line << std::endl;
        }
    } else if (tokens[0] == ".end") {
        /* Ignore */
    } else {
        std::cerr << "warning: unknown command " << std::endl << line << std::endl;
    }
}

Netlist::Netlist(std::string spice_filepath) {
    std::ifstream file(spice_filepath);
    if (!file) {
        throw std::runtime_error("Error: Could not open file " + spice_filepath);
    }

    node_indexes["0"] = 0;

    /* Store title line */
    std::getline(file, title);

    std::string line;
    size_t line_n = 1;
    while (std::getline(file, line)) {
        ++line_n;

        /* Spice is case insensitive */
        to_lowercase(line);
        if (line.size() == 0) continue;

        switch (line[0]) {
            case 'v': /* Voltage source */
                parse_vsource(line);
                break;
            case 'i': /* Current source */
                parse_isource(line);
                break;
            case 'r': /* Resistor */
                parse_resistor(line);
                break;
            case 'c': /* Capacitor */
                parse_capacitor(line);
                break;
            case 'l': /* Inductor */
                parse_inductor(line);
                break;
            case '.': /* Command */
                parse_command(line);
                break;
            case '\n': /* Empty */
            case '*':  /* Comment */
                break;
            default:
                std::cerr << spice_filepath << ": warning: unknown line" << std::endl
                          << line << std::endl;
        }
    }

    std::string dummy = "";
    node_names = std::vector<std::reference_wrapper<const std::string>>{node_indexes.size(), dummy};
    for (auto& [name, index] : node_indexes) {
        if (index >= 0) {  // Include only non-removed nodes
            node_names[index] = std::ref(name);
        }
    }
}

void Netlist::write(std::ostream& os) const {
    os << title << std::endl;

    auto print_component = [&](const Component& c) {
        os << c.name << " " << (c.node1 >= 0 ? node_names[c.node1].get() : "0") << " "
           << (c.node2 >= 0 ? node_names[c.node2].get() : "0") << " " << c.value;
        if (c.pulse) {
            auto p = c.pulse.value();
            os << " " << PulseFunction(p.v1, p.v2, p.td, p.tr, p.tf, p.pw, p.per).to_string();
        }
        os << std::endl;
    };

    for (const auto& c : vsources) print_component(c);
    for (const auto& c : isources) print_component(c);
    for (const auto& c : resistors) print_component(c);
    for (const auto& c : capacitors) print_component(c);

    os << ".tran " << tstep << " " << tstop << std::endl;
    os << ".print tran";
    for (const auto& node : print) os << " v(" << node << ")";
    os << std::endl;

    os << ".end" << std::endl;
}

void Netlist::remove_voltage_sources() {
    std::vector<int> index_updates(node_names.size());

    // Initialize node map with identity mapping
    for (size_t i = 0; i < index_updates.size(); i++) {
        index_updates[i] = i;
    }

    // Find zero voltage sources and create node mapping
    for (const auto& vs : vsources) {
        if (std::abs(vs.value) < 1e-12) {  // Consider values close to 0

            // Always map to the smaller node number
            int min_node = std::min(index_updates[vs.node1], index_updates[vs.node2]);
            int max_node = std::max(index_updates[vs.node1], index_updates[vs.node2]);

            // Update all nodes that mapped to max_node to now map to min_node
            for (int& target : index_updates) {
                if (target == max_node) {
                    target = min_node;
                }
            }

            continue;
        }

        for (auto& r : resistors) {
            std::optional<int> shared_node;
            int pos = -1, neg = -1;

            if (vs.node1 == r.node1) {
                shared_node = vs.node1;
                pos = r.node2;
                neg = vs.node2;
            } else if (vs.node1 == r.node2) {
                shared_node = vs.node1;
                pos = r.node1;
                neg = vs.node2;
            } else if (vs.node2 == r.node1) {
                shared_node = vs.node2;
                pos = vs.node1;
                neg = r.node2;
            } else if (vs.node2 == r.node2) {
                shared_node = vs.node2;
                pos = vs.node1;
                neg = r.node1;
            }

            if (shared_node && shared_node.value() > 0) {
                float current = vs.value / r.value;
                std::string iname = "i_" + vs.name + "_" + r.name;

                isources.emplace_back(iname, neg, pos, current);

                r.node1 = pos;
                r.node2 = neg;

                // Update all nodes that mapped to max_node to now map to -1
                for (int& target : index_updates) {
                    if (target == shared_node.value()) {
                        target = -1;
                    }
                }
            }
        }
    }

    vsources.clear();

    // Compact indices
    std::map<int, int> new_index_map;
    size_t count = 0;
    for (size_t i = 0; i < index_updates.size(); ++i) {
        if (index_updates[i] >= 0) {
            if (new_index_map.find(index_updates[i]) == new_index_map.end()) {
                new_index_map[index_updates[i]] = count++;
            }
            index_updates[i] = new_index_map[index_updates[i]];
        }
    }

    // Update all component nodes
    auto update_nodes = [&](Component& comp) {
        comp.node1 = index_updates[comp.node1];
        comp.node2 = index_updates[comp.node2];
    };

    for (auto& comp : vsources) update_nodes(comp);
    for (auto& comp : isources) update_nodes(comp);
    for (auto& comp : resistors) update_nodes(comp);
    for (auto& comp : capacitors) update_nodes(comp);

    // Update node_indexes and node_names
    std::string dummy = "";
    node_names = std::vector<std::reference_wrapper<const std::string>>{count, dummy};
    for (auto& [name, old_index] : node_indexes) {
        if (old_index < 0) continue;
        int new_index = index_updates[old_index];

        node_indexes[name] = new_index;
        if (new_index >= 0) {
            node_names[new_index] = std::ref(name);
        }
    }
}

void Netlist::handle_zero_voltage_sources(UnionFind& uf) {
    // Union nodes connected by zero voltage sources
    for (const auto& vs : vsources) {
        if (vs.value == 0) {
            uf.unite(vs.node1, vs.node2);
        }
    }
}

void Netlist::handle_capacitors_with_ground_resistors(UnionFind& uf) {
    // Find the representative of ground node (0)
    int ground_rep = uf.find(0);

    // Build a mapping of representative node -> resistors connected to that node
    std::unordered_multimap<int, size_t> node_to_resistors;
    for (size_t i = 0; i < resistors.size(); ++i) {
        const auto& r = resistors[i];
        int rep1 = uf.find(r.node1);
        int rep2 = uf.find(r.node2);
        node_to_resistors.emplace(rep1, i);
        node_to_resistors.emplace(rep2, i);
    }

    for (auto& cap : capacitors) {
        int rep1 = uf.find(cap.node1);
        int rep2 = uf.find(cap.node2);

        // Check if capacitor connects to ground
        if (rep1 == ground_rep || rep2 == ground_rep) {
            continue;
        }

        // Check if capacitor has a neighboring resistor that connects to ground
        auto range1 = node_to_resistors.equal_range(rep1);
        auto range2 = node_to_resistors.equal_range(rep2);

        bool has_gnd_resistor = true;

        // Check resistors connected to node1
        for (auto it = range1.first; it != range1.second; ++it) {
            const auto& r = resistors[it->second];
            int r_rep1 = uf.find(r.node1);
            int r_rep2 = uf.find(r.node2);
            if (r_rep1 != ground_rep && r_rep2 != ground_rep) {
                has_gnd_resistor = false;
                break;
            }
        }

        auto gnd_resistor_range = range1;
        // Check resistors connected to node2 if not already found
        if (!has_gnd_resistor) {
            has_gnd_resistor = true;

            for (auto it = range2.first; it != range2.second; ++it) {
                const auto& r = resistors[it->second];
                int r_rep1 = uf.find(r.node1);
                int r_rep2 = uf.find(r.node2);
                if (r_rep1 != ground_rep && r_rep2 != ground_rep) {
                    has_gnd_resistor = false;
                    break;
                }
            }

            gnd_resistor_range = range2;
            if (!has_gnd_resistor) {
                /* throw std::runtime_error */ std::cerr
                    << ("Cannot eliminate capacitor " + cap.name +
                        ": doesn't connect to gnd, and has no connecting resistor connecting to "
                        "gnd")
                    << std::endl;
            }
        }

        for (auto it = gnd_resistor_range.first; it != gnd_resistor_range.second; ++it) {
            auto& r = resistors[it->second];
            int r_rep1 = uf.find(r.node1);
            int r_rep2 = uf.find(r.node2);

            if (r_rep1 == ground_rep) {
                // Swap r.node 1
                if (r_rep2 == rep1) {
                    r.node1 = rep2;
                    cap.node2 = ground_rep;
                } else {
                    r.node1 = rep1;
                    cap.node1 = ground_rep;
                }
            } else {
                // Swap r.node 2
                if (r_rep1 == rep1) {
                    r.node2 = rep2;
                    cap.node2 = ground_rep;
                } else {
                    r.node2 = rep1;
                    cap.node1 = ground_rep;
                }
            }
        }
    }
}

void Netlist::handle_nonzero_voltage_sources(UnionFind& uf) {
    // Build a mapping of representative node -> resistors connected to that node
    std::unordered_multimap<int, size_t> node_to_resistors;
    for (size_t i = 0; i < resistors.size(); ++i) {
        const auto& r = resistors[i];
        int rep1 = uf.find(r.node1);
        int rep2 = uf.find(r.node2);
        node_to_resistors.emplace(rep1, i);
        node_to_resistors.emplace(rep2, i);
    }

    // Find the representative of ground node (0)
    int ground_rep = uf.find(0);

    for (const auto& vs : vsources) {
        if (vs.value != 0) {  // Non-zero voltage sources

            // Check resistors connected to representative nodes
            int rep1 = uf.find(vs.node1);
            int rep2 = uf.find(vs.node2);
            auto range1 = node_to_resistors.equal_range(rep1);
            auto range2 = node_to_resistors.equal_range(rep2);

            // Choose range as one that is not empty and does not refer to ground
            auto range = range1;
            if (rep1 == ground_rep || range1.first == range1.second) {
                // Use range2 if rep1 is ground or has no resistors
                range = range2;
                if (rep2 == ground_rep || range2.first == range2.second) {
                    throw std::runtime_error(
                        "Cannot eliminate voltage source " + vs.name +
                        ": both representative nodes are ground or have no connected resistors");
                }
            }

            for (auto it = range.first; it != range.second; ++it) {
                auto& r = resistors[it->second];

                // Get representative nodes for comparison
                int r_rep1 = uf.find(r.node1);
                int r_rep2 = uf.find(r.node2);

                std::optional<int> shared_node;
                int pos = -1, neg = -1;

                if (rep1 == r_rep1) {
                    shared_node = vs.node1;
                    pos = r.node2;
                    neg = vs.node2;
                } else if (rep1 == r_rep2) {
                    shared_node = vs.node1;
                    pos = r.node1;
                    neg = vs.node2;
                } else if (rep2 == r_rep1) {
                    shared_node = vs.node2;
                    pos = vs.node1;
                    neg = r.node2;
                } else if (rep2 == r_rep2) {
                    shared_node = vs.node2;
                    pos = vs.node1;
                    neg = r.node1;
                }

                if (shared_node && shared_node.value() > 0) {
                    float current = vs.value / r.value;
                    std::string iname = "I_" + vs.name + "_" + r.name;

                    isources.emplace_back(iname, neg, pos, current);

                    r.node1 = pos;
                    r.node2 = neg;

                    // Unite the shared node with the "removed" node (at index node_names.size())
                    uf.unite(shared_node.value(), node_names.size());
                }
            }
        }
    }
}

void Netlist::compact_and_update_nodes(UnionFind& uf) {
    // Get the mapping from UnionFind and create index updates
    auto index_updates = uf.get_mapping();
    index_updates.resize(node_names.size());  // Trim the extra element
    int removed_representative = uf.find(node_names.size());

    // Convert nodes that map to the "removed" node and compact indices in one pass
    std::map<int, int> new_index_map;
    size_t count = 0;

    // Ensure ground node representative gets index 0
    int ground_rep = uf.find(0);
    if (ground_rep != removed_representative) {
        new_index_map[ground_rep] = count++;
    } else {
        throw std::runtime_error("Ground node was marked to be removed");
    }

    for (auto& update : index_updates) {
        if (update == removed_representative) {
            update = -1;
        } else if (update >= 0) {
            if (new_index_map.find(update) == new_index_map.end()) {
                new_index_map[update] = count++;
            }
            update = new_index_map[update];
        }
    }

    // Update all component nodes
    auto update_nodes = [&](Component& comp) {
        comp.node1 = index_updates[comp.node1];
        comp.node2 = index_updates[comp.node2];

        if (comp.node1 < 0 || comp.node2 < 0 || comp.node1 == comp.node2) {
            throw std::runtime_error("Component " + comp.name +
                                     " has invalid nodes after removing voltage sources.");
        }
    };

    for (auto& comp : isources) update_nodes(comp);
    for (auto& comp : resistors) update_nodes(comp);
    for (auto& comp : capacitors) update_nodes(comp);

    // Update node_indexes and node_names
    std::string dummy = "";
    node_names = std::vector<std::reference_wrapper<const std::string>>{count, dummy};

    // Handle ground node first to ensure it gets index 0
    auto ground_it = node_indexes.find("0");
    node_names[0] = std::ref(ground_it->first);

    for (auto& [name, old_index] : node_indexes) {
        if (old_index <= 0) continue;  // Skip removed nodes and ground
        int new_index = index_updates[old_index];

        node_indexes[name] = new_index;
        if (new_index > 0) {
            node_names[new_index] = std::ref(name);
        }
    }
}

void Netlist::merge_parallel_components() {
    // Merge parallel resistors
    std::vector<Component> merged_resistors;
    std::map<std::pair<int, int>, std::vector<size_t>> resistor_groups;

    // Group resistors by their node pairs (considering both (a,b) and (b,a) as same)
    for (size_t i = 0; i < resistors.size(); ++i) {
        const auto& r = resistors[i];
        std::pair<int, int> nodes = {std::min(r.node1, r.node2), std::max(r.node1, r.node2)};
        resistor_groups[nodes].push_back(i);
    }

    // Merge parallel resistors: 1/R_eq = 1/R1 + 1/R2 + ...
    for (const auto& [nodes, indices] : resistor_groups) {
        if (indices.size() == 1) {
            // Single resistor, keep as is
            merged_resistors.push_back(resistors[indices[0]]);
        } else {
            // Multiple parallel resistors, merge them
            float inverse_sum = 0.0f;
            std::string merged_name = "r_merged";

            for (size_t idx : indices) {
                const auto& r = resistors[idx];
                inverse_sum += 1.0f / r.value;
                merged_name += "_" + r.name;
            }

            float equivalent_resistance = 1.0f / inverse_sum;
            merged_resistors.emplace_back(merged_name, nodes.first, nodes.second,
                                          equivalent_resistance);
        }
    }

    resistors = std::move(merged_resistors);

    // Merge parallel capacitors
    std::vector<Component> merged_capacitors;
    std::map<std::pair<int, int>, std::vector<size_t>> capacitor_groups;

    // Group capacitors by their node pairs (considering both (a,b) and (b,a) as same)
    for (size_t i = 0; i < capacitors.size(); ++i) {
        const auto& c = capacitors[i];
        std::pair<int, int> nodes = {std::min(c.node1, c.node2), std::max(c.node1, c.node2)};
        capacitor_groups[nodes].push_back(i);
    }

    // Merge parallel capacitors: C_eq = C1 + C2 + ...
    for (const auto& [nodes, indices] : capacitor_groups) {
        if (indices.size() == 1) {
            // Single capacitor, keep as is
            merged_capacitors.push_back(capacitors[indices[0]]);
        } else {
            // Multiple parallel capacitors, merge them
            float total_capacitance = 0.0f;
            std::string merged_name = "c_merged";

            for (size_t idx : indices) {
                const auto& c = capacitors[idx];
                total_capacitance += c.value;
                merged_name += "_" + c.name;
            }

            merged_capacitors.emplace_back(merged_name, nodes.first, nodes.second,
                                           total_capacitance);
        }
    }

    capacitors = std::move(merged_capacitors);
}

void Netlist::add_missing_ground_capacitors() {
    // Calculate average capacitance and find nodes with capacitors in one pass
    if (capacitors.empty()) {
        return;  // No capacitors exist, nothing to do
    }

    float total_capacitance = 0.0f;
    std::set<int> nodes_with_capacitors;

    for (const auto& cap : capacitors) {
        total_capacitance += cap.value;
        nodes_with_capacitors.insert(cap.node1);
        nodes_with_capacitors.insert(cap.node2);
    }

    float average_capacitance = total_capacitance / capacitors.size();

    // Find ground node by looking it up in node_indexes
    int ground_node = node_indexes["0"];

    // Add capacitors to ground for nodes that don't have any capacitors
    for (int i = 0; i < static_cast<int>(node_names.size()); ++i) {
        if (i != ground_node && nodes_with_capacitors.find(i) == nodes_with_capacitors.end()) {
            // This node doesn't have any capacitor connected and is not ground
            std::string cap_name = "c_gnd_" + std::to_string(i);
            capacitors.emplace_back(cap_name, i, ground_node, average_capacitance);
        }
    }
}

void Netlist::process_netlist() {
    // Merge parallel components after voltage source removal
    merge_parallel_components();

    UnionFind uf(node_names.size() + 1);  // +1 for the "removed" node

    // Handle zero voltage sources
    handle_zero_voltage_sources(uf);

    // Handle capacitors connecting indirectly to the ground
    handle_capacitors_with_ground_resistors(uf);

    // Handle positive voltage sources
    handle_nonzero_voltage_sources(uf);

    vsources.clear();

    // Compact indices and update all node references
    compact_and_update_nodes(uf);

    // Add missing ground capacitors
    add_missing_ground_capacitors();
}
