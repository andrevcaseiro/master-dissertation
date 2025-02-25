#pragma once
#include <string>
#include <vector>

template <typename T>
class CSRRow;

template <typename T>
class CSRRowIterator;

template <typename T>
class CSREntry;

template <typename T>
class CSRMatrix {
   private:
    int _rows, _columns;
    std::vector<int> _row_pointers;
    std::vector<int> _column_indexes;
    std::vector<T> _values;

   public:
    /**
     * @brief Construct a new CSRMatrix object
     *
     * @param filepath filepath to a csv containing size on first line and values on the following
     * lines
     * @param force_diagonal if true, all diagonal entries are stored even if their value is nul
     */
    CSRMatrix(const std::string& filepath, bool force_diagonal = false);

    /**
     * @brief Number of rows
     *
     * @return int number of rows
     */
    int rows() { return _rows; }

    /**
     * @brief Number of columns
     *
     * @return int number of columns
     */
    int columns() { return _columns; }

    /**
     * @brief Returns a reference to the value at position (row, column)
     *
     * @param row row index
     * @param column column index
     * @return T& reference to the value at (row, column)
     */
    T& at(int row, int column);

    /**
     * @brief Returns the ith diagonal element
     *
     * @param i index
     * @return T the ith diagonal element
     */
    T& diagonal(int i);

    /**
     * @brief Iterator through row entries
     *
     * @param row index
     * @return std::pair<std::vector<Entry<T>>::iterator, std::vector<Entry<T>>::iterator> Iterator
     */
    CSRRow<T> row(int row);

    /**
     * @brief Prints the matrix to cout
     */
    void print();

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
    int _row;

   public:
    CSRRow(CSRMatrix<T>& matrix, int row) : _matrix(matrix), _row(row) {}

    int row() const { return row; }

    CSRRowIterator<T> begin() const { return {_matrix, _row, _matrix._row_pointers[_row]}; }

    CSRRowIterator<T> end() const { return {_matrix, _row, _matrix._row_pointers[_row + 1]}; }
};

template <typename T>
class CSRRowIterator {
   private:
    CSRMatrix<T>& _matrix;
    int _row;
    int _index;

   public:
    CSRRowIterator(CSRMatrix<T>& matrix, int row, int index)
        : _matrix(matrix), _row(row), _index(index) {}

    int row() const { return _row; }

    int column() const { return this->_matrix._column_indexes[_index]; }

    int value() const { return this->_matrix._values[_index]; }

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
    int _row;
    int _index;

   public:
    CSREntry(CSRMatrix<T>& matrix, int row, int index)
        : _matrix(matrix), _row(row), _index(index) {}

    int row() const { return _row; }
    int column() const { return this->_matrix._column_indexes[_index]; }
    T& value() const { return this->_matrix._values[_index]; }
};
