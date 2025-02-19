#include "CSRMatrix.h"

#include <fstream>
#include <iostream>
#include <sstream>

template <typename T>
CSRMatrix<T>::CSRMatrix(std::string filepath) {
    std::ifstream file(filepath);

    std::string line;
    getline(file, line);
    columns = stoi(line);
    rows = columns;

    row_pointers.push_back(0);

    for (int i = 0; i < rows; i++) {
        getline(file, line);

        std::stringstream line_stream(line);
        std::string cell;

        for (int j = 0; j < columns; j++) {
            getline(line_stream, cell, ',');

            std::stringstream cell_stream(cell);

            T value;

            cell_stream >> value;

            if (value != 0) {
                data.emplace_back(j, value);
            }
        }

        row_pointers.push_back(data.size());
    }
}

template <typename T>
T CSRMatrix<T>::diagonal(int i) {
    int j = row_pointers[i];
    int end = row_pointers[i + 1];

    while (j < end && data[j].column < i) {
        j++;
    }

    if (data[j].column == i) return data[j].value;
    return 0;
}

template <typename T>
std::pair<typename std::vector<Entry<T>>::iterator, typename std::vector<Entry<T>>::iterator>
CSRMatrix<T>::row_iterator(int row) {
    int start = row_pointers[row];
    int end = row_pointers[row + 1];

    return {data.begin() + start, data.begin() + end};
}

template <typename T>
void CSRMatrix<T>::print() {
    int data_index = 0;
    for (int i = 0; i < rows; i++) {
        int last = row_pointers[i + 1];

        T value;
        int j;
        for (j = 0; j < columns - 1; j++) {
            if (data_index < last && data[data_index].column == j) {
                value = data[data_index].value;
                data_index++;
            } else {
                value = 0;
            }
            std::cout << value << ", ";
        }

        // use \n instead of , on last cell
        if (data_index < last && data[data_index].column == j) {
            value = data[data_index].value;
            data_index++;
        } else {
            value = 0;
        }
        std::cout << value << std::endl;
    }
}

template <typename T>
void CSRMatrix<T>::print_csr() {
    std::cout << "Row Pointers:   ";
    for (auto& ptr : row_pointers) {
        std::cout << ptr << " ";  // Printing row pointers
    }
    std::cout << std::endl;

    std::cout << "Values:         ";
    for (auto& entry : data) {
        std::cout << entry.value << " ";  // Printing values
    }
    std::cout << std::endl;

    std::cout << "Column Indices: ";
    for (auto& entry : data) {
        std::cout << entry.column << " ";  // Printing column indices
    }
    std::cout << std::endl;
}

template class CSRMatrix<float>;
