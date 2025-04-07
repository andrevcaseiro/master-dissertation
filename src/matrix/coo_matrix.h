#pragma once

#include <list>
#include <unordered_map>

struct Coordinate {
    size_t row, col;

    bool operator==(const Coordinate& other) const { return row == other.row && col == other.col; }
};

struct CoordinateHash {
    std::size_t operator()(const Coordinate& coord) const {
        return std::hash<size_t>{}(coord.row) ^ (std::hash<size_t>{}(coord.col) << 1);
    }
};

template <typename T>
struct COOMatrixEntry {
    size_t row;
    size_t col;
    T value;
};

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

template <typename T>
class COOMatrix {
   private:
    size_t _rows, _cols;
    std::unordered_map<Coordinate, T, CoordinateHash> data;

   public:
    COOMatrix(size_t rows, size_t cols);

    void insert_or_add(size_t row, size_t col, T value);
    size_t increase_size();

    size_t nnz();
    size_t rows();
    size_t cols();

    COOMatrixIterator<T> begin() { return COOMatrixIterator<T>(data.begin()); }
    COOMatrixIterator<T> end() { return COOMatrixIterator<T>(data.end()); }

    void print();
};
