#include "netlist.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <string>

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
        std::cerr << "warning: unknown command: " << std::endl << line << std::endl;
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
            case '.': /* Command */
                parse_command(line);
                break;
            case '\n': /* Empty */
            case '*':  /* Comment */
                break;
            default:
                std::cerr << spice_filepath << ":" << line << "warning: unknown line" << std::endl
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
    os << ".print";
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
                std::string iname = "I_" + vs.name + "_" + r.name;

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
