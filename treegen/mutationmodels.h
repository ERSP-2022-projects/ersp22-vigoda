#ifndef MUTATIONMODELS_H
#define MUTATIONMODELS_H

#include <iostream>
#include "rng.h"
using namespace std;

using mutation_model = int(*)(uint64_t*);

int roll_nucleotide(mutation_model smm, uint64_t* seed) {
    return (smm(seed));
}

int jc69(uint64_t* seed) {
    return nextInt(seed, 4);
}

int k2p(uint64_t* seed) {
    int a = nextInt(seed, 6);
    if (a == 4) a = 0;
    if (a == 5) a = 2;
    return a;
}

mutation_model stomm(string smm) {
    if (smm == "k2p") return k2p;
    if (smm == "jc69") return jc69;
    return jc69;
}

#endif