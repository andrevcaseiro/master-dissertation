/**
 * @file union_find.h
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief Union find data structure for creating and merging nodes in a graph
 * @version 0.1
 * @date 2025-10-24
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#pragma once

#include <vector>

/**
 * @brief Simple Union-Find (Disjoint Set Union) data structure for efficient node merging
 *
 * This implementation uses path compression and union by rank optimizations
 * to achieve nearly constant time operations for practical purposes.
 */
class UnionFind {
   private:
    mutable std::vector<int> parent;
    std::vector<int> rank;

   public:
    /**
     * @brief Construct a new Union Find object
     *
     * @param n Number of elements to initialize
     */
    explicit UnionFind(size_t n);

    /**
     * @brief Find the representative (root) of the set containing x
     *
     * @param x Element to find the representative for
     * @return int Representative of the set containing x
     */
    int find(int x) const;

    /**
     * @brief Unite the sets containing x and y
     *
     * @param x First element
     * @param y Second element
     */
    void unite(int x, int y);

    /**
     * @brief Get the mapping of all elements to their representatives
     *
     * @return std::vector<int> Vector where index i contains the representative of element i
     */
    std::vector<int> get_mapping() const;
};
