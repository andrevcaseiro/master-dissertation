#include "csrd_matrix.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>

template <typename T>
CSRMatrix<T>::CSRMatrix() : _rows(0), _columns(0) {}

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

    for (size_t row = 0; row < _rows; row++) {
        getline(file, line);

        std::stringstream line_stream(line);
        std::string cell;

        /* Add diagonal first */
        _column_indexes.push_back(row);
        _values.push_back(0);

        for (size_t col = 0; col < _columns; col++) {
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
CSRMatrix<T> CSRMatrix<T>::from_coo(const std::string& filepath) {
    CSRMatrix<T> m;

    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Error: Could not open file");
    }

    std::string line;
    getline(file, line);
    std::stringstream ss(line);
    std::string token;

    size_t nnz;
    try {
        std::getline(ss, token, ',');
        m._rows = stol(token);

        std::getline(ss, token, ',');
        m._columns = stol(token);

        std::getline(ss, token, ',');
        nnz = std::stol(token);
    } catch (const std::exception& e) {
        throw std::runtime_error("Error: Failed to read matrix size");
    }

    std::vector<std::list<std::pair<size_t, T>>> data(m._rows, std::list<std::pair<size_t, T>>());

    std::vector<T> diagonal(m._rows, 0);

    for (size_t i = 0; i < nnz; ++i) {
        getline(file, line);

        std::stringstream ss(line);

        size_t row, col;
        T value;

        try {
            std::getline(ss, token, ',');
            row = std::stoi(token);

            std::getline(ss, token, ',');
            col = std::stoi(token);

            std::getline(ss, token, ',');
            value = static_cast<T>(std::stod(token));
        } catch (const std::exception& e) {
            throw std::runtime_error("Error: Failed to parse matrix entries");
        }

        // Correct matlab indexes starting at 1
        --row;
        --col;
        if (row == col) {
            diagonal[row] = value;
            continue;
        }

        data[row].emplace_back(col, value);
    }

    for (size_t row = 0; row < m._rows; ++row) {
        data[row].sort();
    }

    m._row_pointers.reserve(m._rows + 1);
    m._column_indexes.reserve(nnz);
    m._values.reserve(nnz);

    m._row_pointers.emplace_back(0);
    for (size_t row = 0; row < m._rows; ++row) {
        m._column_indexes.emplace_back(row);
        m._values.emplace_back(diagonal[row]);
        for (auto entry : data[row]) {
            m._column_indexes.emplace_back(entry.first);
            m._values.emplace_back(entry.second);
        }
        m._row_pointers.emplace_back(m._values.size());
    }

    return m;
}

template <typename T>
CSRMatrix<T> CSRMatrix<T>::from_coo(COOMatrix<T>& coo_matrix) {
    CSRMatrix<T> m;

    m._rows = coo_matrix.rows();
    m._columns = coo_matrix.cols();

    std::vector<std::list<std::pair<size_t, T>>> data(coo_matrix.rows(),
                                                      std::list<std::pair<size_t, T>>());
    std::vector<T> diagonal(coo_matrix.rows(), 0);

    size_t nnz = coo_matrix.nnz();
    for (auto it : coo_matrix) {
        if (it.row == it.col) {
            diagonal[it.row] = it.value;
            continue;
        }

        data[it.row].emplace_back(it.col, it.value);
    }

    for (size_t row = 0; row < m._rows; ++row) {
        data[row].sort();
    }

    m._row_pointers.reserve(m._rows + 1);
    m._column_indexes.reserve(nnz);
    m._values.reserve(nnz);

    m._row_pointers.emplace_back(0);
    for (size_t row = 0; row < m._rows; ++row) {
        m._column_indexes.emplace_back(row);
        m._values.emplace_back(diagonal[row]);
        for (auto entry : data[row]) {
            m._column_indexes.emplace_back(entry.first);
            m._values.emplace_back(entry.second);
        }
        m._row_pointers.emplace_back(m._values.size());
    }

    return m;
}

template <typename T>
T& CSRMatrix<T>::at(size_t row, size_t column) {
    if (row == column) {
        return diagonal(row);
    }

    auto row_begin = _column_indexes.begin() + _row_pointers[row] + 1;
    auto row_end = _column_indexes.begin() + _row_pointers[row + 1];

    auto res = std::lower_bound(row_begin, row_end, column);

    if (res == row_end || *res != column) {
        throw std::domain_error("Attempted to access a null entry of a sparse matrix at " +
                                std::to_string(row) + "," + std::to_string(column) + ".");
    }

    return _values[std::distance(_column_indexes.begin(), res)];
}

template <typename T>
T& CSRMatrix<T>::diagonal(size_t i) {
    return _values[_row_pointers[i]];
}

template <typename T>
CSRRow<T> CSRMatrix<T>::row(size_t row) {
    return CSRRow<T>(*this, row);
}

template <typename T>
void CSRMatrix<T>::print(int prec) {
    std::cout << std::fixed << std::setprecision(prec);
    for (size_t row = 0; row < rows(); row++) {
        for (size_t col = 0; col < columns(); col++) {
            try {
                std::cout << std::setw(prec + 3) << at(row, col);
            } catch (const std::domain_error& e) {
                std::cout << std::setw(prec + 3) << 0.0;
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
