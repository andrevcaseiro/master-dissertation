#include "csrd_matrix.h"

#include <fstream>
#include <iostream>
#include <sstream>

template <typename T>
CSRMatrix<T>::CSRMatrix(const std::string& filepath) {
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

        /* Add diagonal first */
        _column_indexes.push_back(row);
        _values.push_back(0);

        for (int col = 0; col < _columns; col++) {
            getline(line_stream, cell, ',');

            std::stringstream cell_stream(cell);

            T value;

            cell_stream >> value;

            if (row == col) {
                int diagonal_pos = _row_pointers.back();
                _values[diagonal_pos] = value;
                continue;
            }

            if (value != 0) {
                _column_indexes.push_back(col);
                _values.push_back(value);
            }
        }

        _row_pointers.push_back(_values.size());
    }
}

template <typename T>
T& CSRMatrix<T>::at(int row, int column) {
    if (row == column) {
        return diagonal(row);
    }

    auto row_begin = _column_indexes.begin() + _row_pointers[row] + 1;
    auto row_end = _column_indexes.begin() + _row_pointers[row + 1];

    auto res = std::lower_bound(row_begin, row_end, column);

    if (*res != column) {
        throw std::domain_error("Attempted to access a null entry of a sparse matrix at " +
                                std::to_string(row) + "," + std::to_string(column) + ".");
    }

    return _values[std::distance(_column_indexes.begin(), res)];
}

template <typename T>
T& CSRMatrix<T>::diagonal(int i) {
    return _values[_row_pointers[i]];
}

template <typename T>
CSRRow<T> CSRMatrix<T>::row(int row) {
    return CSRRow<T>(*this, row);
}

template <typename T>
void CSRMatrix<T>::print() {
    for (int row = 0; row < rows(); row++) {
        for (int col = 0; col < columns(); col++) {
            try {
                std::cout << at(row, col);
            } catch (const std::domain_error& e) {
                std::cout << 0;
            }

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
