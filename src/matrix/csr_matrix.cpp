#include "csr_matrix.h"

#include <fstream>
#include <iostream>
#include <sstream>

template <typename T>
CSRMatrix<T>::CSRMatrix(const std::string& filepath, bool force_diagonal) {
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Error: Could not open file");
    }

    std::string line;
    getline(file, line);

    try {
        _rows = _columns = stoi(line);
    } catch (const std::exception& e) {
        throw std::runtime_error("Error: Failed to read matrix size");
    }

    _row_pointers.push_back(0);

    for (int row = 0; row < _rows; row++) {
        getline(file, line);

        std::stringstream line_stream(line);
        std::string cell;

        for (int col = 0; col < _columns; col++) {
            getline(line_stream, cell, ',');

            std::stringstream cell_stream(cell);

            T value;

            cell_stream >> value;

            if (value != 0 || (force_diagonal && row == col)) {
                _column_indexes.push_back(col);
                _values.push_back(value);
            }
        }

        _row_pointers.push_back(_values.size());
    }
}

template <typename T>
T& CSRMatrix<T>::at(int row, int column) {
    auto row_begin = _column_indexes.begin() + _row_pointers[row];
    auto row_end = _column_indexes.begin() + _row_pointers[row + 1];

    auto res = std::lower_bound(row_begin, row_end, column);

    if (*res != column) {
        throw std::domain_error("Attempted to access a null entry of a sparse matrix.");
    }

    return _values[std::distance(_column_indexes.begin(), res)];
}

template <typename T>
T& CSRMatrix<T>::diagonal(int i) {
    return this->at(i, i);
}

template <typename T>
CSRRow<T> CSRMatrix<T>::row(int row) {
    return CSRRow<T>(*this, row);
}

template <typename T>
void CSRMatrix<T>::print() {
    int curr_index = 0;
    for (int row = 0; row < rows(); row++) {
        int last_row_index = _row_pointers[row + 1];

        for (int col = 0; col < columns(); col++) {
            T value;
            if (curr_index < last_row_index && _column_indexes[curr_index] == col) {
                value = _values[curr_index];
                curr_index++;
            } else {
                value = 0;
            }
            std::cout << value;

            if (col != columns() - 1) {
                std::cout << ",";
            } else {
                std::cout << std::endl;
            }
        }
    }
}

template <typename T>
void CSRMatrix<T>::print_csr() {
    std::cout << "Row Pointers:   ";
    for (auto& ptr : _row_pointers) {
        std::cout << ptr << " ";
    }
    std::cout << std::endl;

    std::cout << "Column Indices: ";
    for (auto& idx : _column_indexes) {
        std::cout << idx << " ";
    }
    std::cout << std::endl;

    std::cout << "Values:         ";
    for (auto& value : _values) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
}

template class CSRMatrix<float>;
