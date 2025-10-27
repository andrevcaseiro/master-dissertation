/**
 * @file read_vector.cpp
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief Utility function to read a vector from a file
 * @version 0.1
 * @date 2025-03-24
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#pragma once

#include <string>
#include <vector>

/**
 * @brief Reads a float vector from a file with size on first line
 * and comma-separated values on second line
 * 
 * @param filepath filepath to read
 * @return std::vector<float> vector from the file
 */
std::vector<float> read_vector(std::string filepath);
