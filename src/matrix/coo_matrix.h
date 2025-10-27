/**
 * @file coo_matrix.h
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief A class for sparse matrices in coordinate format
 * @version 0.1
 * @date 2025-10-24
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#pragma once

#include <list>
#include <unordered_map>
#include <cstddef>

/**
 * @brief Matrix coordinate
 * 
 */
struct Coordinate {
    size_t row, col;

    bool operator==(const Coordinate& other) const { return row == other.row && col == other.col; }
};

/**
 * @brief Hash function for matrix coordinates
 * 
 */
struct CoordinateHash {
    std::size_t operator()(const Coordinate& coord) const {
        return std::hash<size_t>{}(coord.row) ^ (std::hash<size_t>{}(coord.col) << 1);
    }
};

/**
 * @brief Entry of a matrix in coordinate format
 * 
 * @tparam T value type
 */
template <typename T>
struct COOMatrixEntry {
    size_t row;
    size_t col;
    T value;
};

/**
 * @brief Iterator over the entries of a matrix in coordinate format
 * 
 * @tparam T value type
 */
template <typename T>
class COOMatrixIterator {
   private:
    typename std::unordered_map<Coordinate, T, CoordinateHash>::iterator it;

   public:
    COOMatrixIterator(typename std::unordered_map<Coordinate, T, CoordinateHash>::iterator it)
        : it(it) {}

    COOMatrixEntry<T> operator*() const {
        return COOMatrixEntry<T>{it->first.row, it->first.col, it->second};
    }

    bool operator!=(const COOMatrixIterator& other) const { return it != other.it; }

    COOMatrixIterator& operator++() {
        ++it;
        return *this;
    }
};

/**
 * @brief Matrix in coordinate format
 * 
 * @tparam T value type
 */
template <typename T>
class COOMatrix {
   private:
    size_t _rows, _cols;
    std::unordered_map<Coordinate, T, CoordinateHash> data;

   public:
    /**
     * @brief Construct a new empty COOMatrix object
     * 
     * @param rows number of rows
     * @param cols number of columns
     */
    COOMatrix(size_t rows, size_t cols);

    /**
     * @brief Add value to an entry, creating it if it doesnt exist 
     * 
     * @param row 
     * @param col 
     * @param value 
     */
    void insert_or_add(size_t row, size_t col, T value);

    /**
     * @brief Increment the number of rows and collumns by one
     * 
     * @return size_t new matrix size
     */
    size_t increase_size();

    /**
     * @brief Get the number of non-zero entries
     * 
     * @return size_t number of non-zero entries
     */
    size_t nnz();

    /**
     * @brief Get the number of rows
     * 
     * @return size_t number of rows
     */
    size_t rows();

    /**
     * @brief Get the number of columns
     * 
     * @return size_t number of columns
     */
    size_t cols();

    /**
     * @brief Iterator to the begin of the unordered list of entries
     * 
     * @return COOMatrixIterator<T> 
     */
    COOMatrixIterator<T> begin() { return COOMatrixIterator<T>(data.begin()); }

    /**
     * @brief Iterator to the end of the unordered list of entries
     * 
     * @return COOMatrixIterator<T> 
     */
    COOMatrixIterator<T> end() { return COOMatrixIterator<T>(data.end()); }

    /**
     * @brief Print the matrix to stdout as a list of triplets
     * 
     */
    void print();
};
