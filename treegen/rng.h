#ifndef RNG_H
#define RNG_H

#include <iostream>
using namespace std;

static inline uint64_t init(uint64_t seed) {
    return (seed ^ 0x5deece66d) & (0xffffffffffff);
}

static inline int next(uint64_t *seed, const int bits) {
    *seed = (*seed * 0x5deece66d + 0xb) & (0xffffffffffff);
    return (int) ((int64_t)*seed >> (48 - bits));
}

static inline int nextInt(uint64_t *seed, const int n) {
    int bits, val;
    const int m = n - 1;
    if ((m & n) == 0) {
        uint64_t x = n * (uint64_t)next(seed, 31);
        return (int) ((int64_t) x >> 31);
    }
    do {
        bits = next(seed, 31);
        val = bits % n;
    } while (bits - val + m < 0);
    return val;
}

static inline float nextFloat(uint64_t *seed) {
    return next(seed, 24) / (float) (1 << 24);
}

#endif