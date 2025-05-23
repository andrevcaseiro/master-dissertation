#include "coo_matrix.h"

#include <iostream>

template <typename T>
COOMatrix<T>::COOMatrix(size_t rows, size_t cols) : _rows(rows), _cols(cols) {}

template <typename T>
void COOMatrix<T>::insert_or_add(size_t row, size_t col, T value) {
    if (row >= _rows || col >= _cols) {
        throw std::out_of_range("Index out of bounds: " + std::to_string(row) + ", " +
                                std::to_string(col));
    }
    Coordinate key = {row, col};
    data[key] += value;
}

template <typename T>
size_t COOMatrix<T>::increase_size() {
    _rows += 1;
    _cols += 1;
    return _rows;
}

template <typename T>
size_t COOMatrix<T>::nnz() {
    return data.size();
}

template <typename T>
size_t COOMatrix<T>::rows() {
    return _rows;
}

template <typename T>
size_t COOMatrix<T>::cols() {
    return _cols;
}

template <typename T>
void COOMatrix<T>::print() {
    for (auto it : *this) {
        std::cout << it.row << ", " << it.col << ", " << it.value << std::endl;
    }
}

template class COOMatrix<float>;
