#include "union_find.h"

#include <numeric>

UnionFind::UnionFind(size_t n) : parent(n), rank(n, 0) {
    std::iota(parent.begin(), parent.end(), 0);
}

int UnionFind::find(int x) const {
    if (parent[x] != x) {
        parent[x] = find(parent[x]);  // Path compression
    }
    return parent[x];
}

void UnionFind::unite(int x, int y) {
    int px = find(x);
    int py = find(y);

    if (px == py) return;

    // Union by rank: attach smaller rank tree under root of higher rank tree
    if (rank[px] < rank[py]) {
        parent[px] = py;
    } else if (rank[px] > rank[py]) {
        parent[py] = px;
    } else {
        parent[py] = px;
        rank[px]++;
    }
}

std::vector<int> UnionFind::get_mapping() const {
    std::vector<int> mapping(parent.size());
    for (size_t i = 0; i < parent.size(); ++i) {
        mapping[i] = find(i);
    }
    return mapping;
}
