/**
 * @file read_vector.cpp
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief Read vector from file
 * @version 0.1
 * @date 2025-03-24
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include "read_vector.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

std::vector<float> read_vector(std::string filepath) {
    std::vector<float> v;
    std::ifstream file(filepath);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open the file.");
    }

    std::string line;

    // Read the length
    if (!std::getline(file, line)) {
        throw std::runtime_error("File is empty.");
    }
    int expected_length = 0;
    std::istringstream(line) >> expected_length;
    v.reserve(expected_length);

    // Read the second line for the comma-separated values
    if (!std::getline(file, line)) {
        throw std::runtime_error("Vector data missing.");
    }

    std::istringstream ss(line);
    std::string value;
    while (std::getline(ss, value, ',')) {
        try {
            v.push_back(std::stof(value));  // Convert string to float and add to vector
        } catch (const std::invalid_argument& e) {
            throw std::runtime_error("Invalid float value in the CSV file.");
        }
    }

    // Check if the number of elements matches the expected length
    if (v.size() != (size_t)expected_length) {
        throw std::runtime_error("Vector data does not match the expected length.");
    }

    return v;
}
