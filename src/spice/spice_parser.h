#pragma once

#include <cmath>
#include <list>
#include <string>
#include <unordered_map>

#include "../matrix/csrd_matrix.h"
#include "time_function.h"

class SpiceParser {
   public:
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

   private:
    std::vector<float> C;
    CSRMatrix<float> G;
    std::vector<std::unique_ptr<TimeFunction>> b;
    std::list<Component> vsources;
    std::list<Component> isources;
    std::list<Component> resistors;
    std::list<Component> capacitors;
    std::unordered_map<std::string, int> node_map;

    void parse_file(std::string filepath);
    void gen_mna();
    int get_node_index(const std::string& node);

   public:
    SpiceParser(std::string filepath);

    std::vector<float>& get_C() { return C; }
    CSRMatrix<float>& get_G() { return G; }
    std::vector<std::unique_ptr<TimeFunction>>& get_b() { return b; }
};
