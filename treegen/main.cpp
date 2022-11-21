#include "generator.h"
#include <iostream>
using namespace std;

int main(int argc, char** argv) {
    uint64_t seed = (argc >= 2) ? stoll(argv[1]) : time(0);
    Generator g(seed);
    g.generate_tree();
    return 0;
}