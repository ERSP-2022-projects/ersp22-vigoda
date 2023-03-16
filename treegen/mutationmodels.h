#ifndef MUTATIONMODELS_H
#define MUTATIONMODELS_H

#include <iostream>
#include "rng.h"
using namespace std;

using mutation_model = int(*)(uint64_t*);

// note purines % 2 = 0, pyrimidines % 2 = 1
static inline char convert_nucleotide(int id) {
    if (id == 0)
        return 'A'; // pairs T, purine
    if (id == 1)
        return 'T'; // pairs A, pyrimidine
    if (id == 2)
        return 'G'; // pairs C, purine
    if (id == 3)
        return 'C'; // pairs G, pyrimidine
    return '-'; // unknown or missing data
}

static inline int roll_nucleotide(mutation_model smm, uint64_t* seed) {
    return (smm(seed));
}

static inline int jc69(uint64_t* seed) {
    return nextInt(seed, 4);
}

static inline int k2p(uint64_t* seed) {
    int a = nextInt(seed, 6);
    if (a == 4) a = 0;
    if (a == 5) a = 2;
    return a;
}

static inline mutation_model stomm(string smm) {
    if (smm == "k2p") return k2p;
    if (smm == "jc69") return jc69;
    return jc69;
}

#endif