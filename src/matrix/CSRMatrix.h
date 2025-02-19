#pragma once
#include <string>
#include <vector>

template <typename T>
struct Entry {
    int column;
    T value;

    Entry(int column, T value) : column(column), value(value) {}
};

template <typename T>
class CSRMatrix {
   private:
    int columns, rows;
    std::vector<T> row_pointers;
    std::vector<Entry<T>> data;

   public:
    /**
     * @brief Construct a new CSRMatrix object
     *
     * @param filepath filepath to a csv containing size on first line and values on the following
     * lines
     */
    CSRMatrix(std::string filepath);

    /**
     * @brief Number of rows
     * 
     * @return int number of rows
     */
    int get_rows() { return rows; }

    /**
     * @brief Number of columns
     * 
     * @return int number of columns
     */
    int get_columns() { return columns; }

    /**
     * @brief Returns the ith diagonal element
     *
     * @param i index
     * @return T the ith diagonal element
     */
    T diagonal(int i);

    /**
     * @brief Iterator through row entries
     *
     * @param row index
     * @return std::pair<std::vector<Entry<T>>::iterator, std::vector<Entry<T>>::iterator> Iterator
     */
    std::pair<typename std::vector<Entry<T>>::iterator, typename std::vector<Entry<T>>::iterator>
    row_iterator(int row);

    /**
     * @brief Prints the matrix to cout
     */
    void print();

    /**
     * @brief Prints the CSR structures to cout
     */
    void print_csr();
};
