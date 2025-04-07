#pragma once
#include <string>
#include <vector>

#include "coo_matrix.h"

template <typename T>
class CSRRow;

template <typename T>
class CSRRowIterator;

template <typename T>
class CSREntry;

template <typename T>
/**
 * @brief Compressed Sparce Row matrix, with diagonals stored as the first element from each row,
 * allowing access to diagonals in O(1) and returning iterators to non diagonal entries
 *
 */
class CSRMatrix {
   private:
    size_t _rows, _columns;
    std::vector<size_t> _row_pointers;
    std::vector<size_t> _column_indexes;
    std::vector<T> _values;

   public:
    CSRMatrix();

    /**
     * @brief Construct a new CSRMatrix object
     *
     * @param filepath filepath to a csv containing size on first line and values on the following
     * lines
     * @param force_diagonal if true, all diagonal entries are stored even if their value is nul
     */
    CSRMatrix(const std::string& filepath);

    /**
     * @brief Constructs a new CSRMatrix from a file describing the matrix in coordinate format
     *
     * The first line of the file should contain:
     *
     * - rows, columns, nnz
     *
     * Followed by nnz lines containing:
     *
     * - row, column, value
     *
     * @param filepath
     * @return CSRMatrix<T>
     */
    static CSRMatrix<T> from_coo(const std::string& filepath);

    /**
     * @brief Constructs a new CSRMatrix from a COOMatrix
     *
     * @param matrix sparse matrix in coordinate format
     * @return CSRMatrix<T>
     */
    static CSRMatrix<T> from_coo(COOMatrix<T>& matrix);

    /**
     * @brief Number of rows
     *
     * @return int number of rows
     */
    size_t rows() { return _rows; }

    /**
     * @brief Number of columns
     *
     * @return int number of columns
     */
    size_t columns() { return _columns; }

    /**
     * @brief Returns a reference to the value at position (row, column)
     *
     * @param row row index
     * @param column column index
     * @return T& reference to the value at (row, column)
     */
    T& at(size_t row, size_t column);

    /**
     * @brief Returns the ith diagonal element
     *
     * @param i index
     * @return T the ith diagonal element
     */
    T& diagonal(size_t i);

    /**
     * @brief Iterator through row entries
     *
     * @param row index
     * @return std::pair<std::vector<Entry<T>>::iterator, std::vector<Entry<T>>::iterator> Iterator
     */
    CSRRow<T> row(size_t row);

    /**
     * @brief Prints the matrix to cout
     */
    void print(int prec = 3);

    /**
     * @brief Prints the CSR structures to cout
     */
    void print_csr();

    friend class CSRRow<T>;
    friend class CSRRowIterator<T>;
    friend class CSREntry<T>;
};

template <typename T>
class CSRRow {
    CSRMatrix<T>& _matrix;
    size_t _row;

   public:
    CSRRow(CSRMatrix<T>& matrix, size_t row) : _matrix(matrix), _row(row) {}

    size_t row() const { return row; }

    CSRRowIterator<T> begin() const { return {_matrix, _row, _matrix._row_pointers[_row]}; }

    CSRRowIterator<T> end() const { return {_matrix, _row, _matrix._row_pointers[_row + 1]}; }

    auto begin_off_diagonal_values() const {
        return _matrix._values.begin() + _matrix._row_pointers[_row] + 1;
    }
    auto end_off_diagonal_values() const {
        return _matrix._values.begin() + _matrix._row_pointers[_row + 1];
    }
    auto begin_off_diagonal_columns() const {
        return _matrix._column_indexes.begin() + _matrix._row_pointers[_row] + 1;
    }
    auto end_off_diagonal_columns() const {
        return _matrix._column_indexes.begin() + _matrix._row_pointers[_row + 1];
    }
};

template <typename T>
class CSRRowIterator {
   private:
    CSRMatrix<T>& _matrix;
    size_t _row;
    size_t _index;

   public:
    CSRRowIterator(CSRMatrix<T>& matrix, size_t row, size_t index)
        : _matrix(matrix), _row(row), _index(index) {}

    size_t row() const { return _row; }

    size_t column() const { return this->_matrix._column_indexes[_index]; }

    size_t value() const { return this->_matrix._values[_index]; }

    CSREntry<T> operator*() const { return CSREntry<T>(_matrix, _row, _index); }

    CSRRowIterator& operator++() {
        ++_index;
        return *this;
    }

    bool operator!=(const CSRRowIterator& other) const { return _index != other._index; }
};

template <typename T>
class CSREntry {
   private:
    CSRMatrix<T>& _matrix;
    size_t _row;
    size_t _index;

   public:
    CSREntry(CSRMatrix<T>& matrix, size_t row, size_t index)
        : _matrix(matrix), _row(row), _index(index) {}

    size_t row() const { return _row; }
    size_t column() const { return this->_matrix._column_indexes[_index]; }
    T& value() const { return this->_matrix._values[_index]; }
};
