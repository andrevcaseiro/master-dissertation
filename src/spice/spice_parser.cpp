#include "spice_parser.h"

#include <fstream>
#include <regex>
#include <sstream>

#include "../matrix/coo_matrix.h"
#include "time_function.h"

#define FLOAT_PATTERN R"(([\d.eE+-]+))"  // Matches numbers, including scientific notation
#define COMMA_PATTERN R"(\s*)"       // Matches a comma with optional spaces
#define PULSE_START R"(pulse\(\s*)"
#define PULSE_END R"(\s*\))"

#define PULSE_REGEX                                                                            \
    std::string(PULSE_START) + FLOAT_PATTERN + COMMA_PATTERN + FLOAT_PATTERN + COMMA_PATTERN + \
        FLOAT_PATTERN + COMMA_PATTERN + FLOAT_PATTERN + COMMA_PATTERN + FLOAT_PATTERN +        \
        COMMA_PATTERN + FLOAT_PATTERN + COMMA_PATTERN + FLOAT_PATTERN + PULSE_END

SpiceParser::SpiceParser(std::string filepath) {
    parse_file(filepath);
    gen_mna();
}

void SpiceParser::parse_file(std::string filepath) {
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Error: Could not open file");
    }

    // Set ground to index 0
    node_map["0"] = -1;

    std::string line;

    /* Discard title line */
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string name;
        iss >> name;

        char c = std::tolower(name[0]);

        std::string node1_str, node2_str;
        int node1, node2;
        float value;
        switch (c) {
            case 'v': /* Voltage source */
                iss >> node1_str >> node2_str >> value;
                node1 = get_node_index(node1_str);
                node2 = get_node_index(node2_str);

                {
                    std::string rest_of_line;
                    std::getline(iss, rest_of_line);

                    static const std::regex pulse_regex(PULSE_REGEX);
                    std::smatch match;

                    if (std::regex_search(rest_of_line, match, pulse_regex)) {
                        // Extract pulse parameters from the regex match
                        float v1 = std::stof(match[1]);
                        float v2 = std::stof(match[2]);
                        float td = std::stof(match[3]);
                        float tr = std::stof(match[4]);
                        float tf = std::stof(match[5]);
                        float pw = std::stof(match[6]);
                        float per = std::stof(match[7]);

                        vsources.emplace_back(name, node1, node2, value, v1, v2, td, tr, tf, pw,
                                              per);
                    } else {
                        // Regular voltage source (no pulse function)
                        vsources.emplace_back(name, node1, node2, value);
                    }
                }
                break;
            case 'i': /* Current source */
                iss >> node1_str >> node2_str >> value;
                node1 = get_node_index(node1_str);
                node2 = get_node_index(node2_str);
                {
                    std::string rest_of_line;
                    std::getline(iss, rest_of_line);

                    static const std::regex pulse_regex(PULSE_REGEX);
                    std::smatch match;

                    if (std::regex_search(rest_of_line, match, pulse_regex)) {
                        // Extract pulse parameters from the regex match
                        float v1 = std::stof(match[1]);
                        float v2 = std::stof(match[2]);
                        float td = std::stof(match[3]);
                        float tr = std::stof(match[4]);
                        float tf = std::stof(match[5]);
                        float pw = std::stof(match[6]);
                        float per = std::stof(match[7]);

                        isources.emplace_back(name, node1, node2, value, v1, v2, td, tr, tf, pw,
                                              per);
                    } else {
                        // Regular current source (no pulse function)
                        isources.emplace_back(name, node1, node2, value);
                    }
                }

                break;
            case 'r': /* Resistor */
                iss >> node1_str >> node2_str >> value;
                node1 = get_node_index(node1_str);
                node2 = get_node_index(node2_str);
                resistors.emplace_back(name, node1, node2, value);
                break;
            case 'c': /* Capacitor */
                iss >> node1_str >> node2_str >> value;
                node1 = get_node_index(node1_str);
                node2 = get_node_index(node2_str);
                capacitors.emplace_back(name, node1, node2, value);
                break;
            case '*': /* Comment */
                break;
            case '.': /* Command */
                      /* TODO */
                break;
            default:
                break;
        }
    }
}

void SpiceParser::gen_mna() {
    size_t size = node_map.size() - 1;

    C = std::vector<float>(size, 0);
    b = std::vector<std::unique_ptr<TimeFunction>>(size);
    for (auto& ptr : b) {
        ptr = std::make_unique<ConstantFunction>(0.0);
    }

    COOMatrix<float> G_coo(size, size);
    for (Component c : resistors) {
        float conductance = 1.0f / c.value;

        if (c.node1 >= 0) G_coo.insert_or_add(c.node1, c.node1, conductance);
        if (c.node2 >= 0) G_coo.insert_or_add(c.node2, c.node2, conductance);

        if (c.node1 >= 0 && c.node2 >= 0) {
            G_coo.insert_or_add(c.node2, c.node1, -conductance);
            G_coo.insert_or_add(c.node1, c.node2, -conductance);
        }
    }

    for (Component c : capacitors) {
        if (c.node1 != -1 && c.node2 != -1) {
            throw std::runtime_error("Capacitor " + c.name + " must be connected to GND.");
        }

        int node = c.node1 == -1 ? c.node2 : c.node1;

        C[node] += c.value;
    }

    for (Component c : vsources) {
        size_t new_row = G_coo.increase_size() - 1;
        
        C.push_back(0);

        if (c.pulse) {
            const PulseParams& p = c.pulse.value();
            b.push_back(std::make_unique<PulseFunction>(p.v1, p.v2, p.td, p.tr, p.tf, p.pw, p.per));
        } else {
            b.push_back(std::make_unique<ConstantFunction>(c.value));
        }

        /* Set signal of i */
        if (c.node1 >= 0) {
            G_coo.insert_or_add(new_row, c.node1, -1.0f);
            G_coo.insert_or_add(c.node1, new_row, -1.0f);
            *b[new_row] *= -1;
        }
        if (c.node2 >= 0) {
            G_coo.insert_or_add(new_row, c.node2, -1.0f);
            G_coo.insert_or_add(c.node2, new_row, -1.0f);
        }
    }

    for (Component c : isources) {
        if (c.pulse) {
            const PulseParams& p = c.pulse.value();
            if (c.node2 >= 0)
                b[c.node2] = *b[c.node2] + PulseFunction(p.v1, p.v2, p.td, p.tr, p.tf, p.pw, p.per);
            if (c.node1 >= 0)
                b[c.node1] =
                    *b[c.node1] + PulseFunction(-p.v1, -p.v2, p.td, p.tr, p.tf, p.pw, p.per);
        } else {
            if (c.node2 >= 0) *b[c.node2] += c.value;
            if (c.node1 >= 0) *b[c.node1] += -c.value;
        }
    }

    G = CSRMatrix<float>::from_coo(G_coo);
}

int SpiceParser::get_node_index(const std::string& node) {
    if (node_map.find(node) == node_map.end()) {
        node_map[node] = node_map.size() - 1;  // Assign a new index, GND is -1
    }
    return node_map[node];
}
