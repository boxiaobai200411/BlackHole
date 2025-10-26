#ifndef TYPE_H
#define TYPE_H

#include <algorithm>
#include <cstdint>
#include <map>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
typedef unsigned int ui;
typedef ui LabelID;
typedef uint32_t VertexID;
typedef uint64_t ul;
#include <cmath>
using namespace std;

class BloomFilter {
  private:
    size_t m;
    size_t k;
    std::vector<uint64_t> bits;

    inline bool getBit(size_t idx) const {
        return (bits[idx >> 6] >> (idx & 63)) & 1ULL;
    }

    inline void setBit(size_t idx) { bits[idx >> 6] |= (1ULL << (idx & 63)); }

    static uint64_t hash1(uint64_t x) {
        x ^= (x >> 33);
        x *= 0xff51afd7ed558ccdULL;
        x ^= (x >> 33);
        x *= 0xc4ceb9fe1a85ec53ULL;
        x ^= (x >> 33);
        return x;
    }

    static uint64_t hash2(uint64_t x) {
        x = (~x) + (x << 21);
        x ^= (x >> 24);
        x = (x + (x << 3)) + (x << 8);
        x ^= (x >> 14);
        x = (x + (x << 2)) + (x << 4);
        x ^= (x >> 28);
        x += (x << 31);
        return x;
    }

    std::pair<uint64_t, uint64_t> aggregate(const VertexID *arr, ui n) const {
        uint64_t h1 = 0, h2 = 0;
        for (ui i = 0; i < n; i++) {
            uint64_t v = arr[i];
            h1 ^= hash1(v);
            h2 += hash2(v);
        }
        return {h1, h2};
    }

  public:
    BloomFilter(size_t m_bits, size_t k_hash) : m(m_bits), k(k_hash) {
        bits.resize((m + 63) / 64);
    }

    BloomFilter(size_t expected_elements, double false_positive_rate = 0.01) {
        if (expected_elements == 0)
            expected_elements = 100000000000000ll;
        if (false_positive_rate <= 0 || false_positive_rate >= 1)
            false_positive_rate = 0.01;

        m = static_cast<size_t>(
            ceil(-(double)expected_elements * log(false_positive_rate) /
                 (log(2) * log(2))));
        k = static_cast<size_t>(round((double)m / expected_elements * log(2)));
        if (k < 1)
            k = 1;
        if (k > 30)
            k = 30;

        bits.resize((m + 63) / 64);
    }

    BloomFilter() : BloomFilter(1000000, 0.01) {}

    bool possiblyContainsArray(const VertexID *arr, ui n) const {
        auto [h1, h2] = aggregate(arr, n);
        for (size_t i = 0; i < k; i++) {
            size_t pos = (h1 + i * h2) % m;
            if (!getBit(pos))
                return false;
        }
        return true;
    }

    void insertArray(const VertexID *arr, ui n) {
        auto [h1, h2] = aggregate(arr, n);
        for (size_t i = 0; i < k; i++) {
            size_t pos = (h1 + i * h2) % m;
            setBit(pos);
        }
    }

    bool checkAndInsertArray(const VertexID *arr, ui n) {
        auto [h1, h2] = aggregate(arr, n);
        bool exists = true;

        for (size_t i = 0; i < k; i++) {
            size_t pos = (h1 + i * h2) % m;
            if (!getBit(pos)) {
                exists = false;
                break;
            }
        }

        if (!exists) {
            for (size_t i = 0; i < k; i++) {
                size_t pos = (h1 + i * h2) % m;
                setBit(pos);
            }
        }

        return exists;
    }

    void clear() { std::fill(bits.begin(), bits.end(), 0ULL); }
};
inline int checkRelation(vector<ui> &a, vector<ui> &b) {
    size_t i = 0, j = 0;
    bool hasCommon = false;
    ui lenA = a.size();
    ui lenB = b.size();
    while (i < lenA && j < lenB) {
        if (a[i] == b[j]) {
            hasCommon = true;
            ++i;
            ++j;
        } else if (a[i] < b[j]) {
            ++i;
        } else {
            ++j;
        }
    }

    if (!hasCommon)
        return 0;

    i = j = 0;
    while (i < lenA && j < lenB) {
        if (a[i] == b[j]) {
            ++i;
            ++j;
        } else if (a[i] < b[j]) {
            ++i;
        } else {
            break;
        }
    }
    if (j == lenB)
        return 2;

    return 1;
}

class UnionFind {
  private:
    std::vector<int> parent;
    std::vector<int> rank;

  public:
    UnionFind(int n) : parent(n) {
        for (int i = 0; i < n; ++i)
            parent[i] = i;
    }

    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    }

    void unite(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);
        if (rootX == rootY)
            return;
        parent[rootY] = rootX;
    }

    bool connected(int x, int y) { return find(x) == find(y); }

    int countRoots() {
        int count = 0;
        for (int i = 0; i < parent.size(); ++i) {
            if (parent[i] == i)
                count++;
        }
        return count;
    }
};

#endif