#include <iostream>
#include <numeric>
#include <LUDecomposition.h>
#include <norm.h>
#include <resol.h>
#include "io.h"


inline
void printError(const char *message, const char *arg = nullptr) {
    std::cerr << "error: " << message;
    if (arg) std::cerr << ' ' << arg;
    std::cerr << '\n';
}

inline
void handleArg(const char *arg) {
    auto oName = solveFile(arg);
    std::cout << "written solution of " << arg << " to " << oName << std::endl;
}

inline void handleArg2(const char *arg) {
    auto file = openFile(arg);
    Matrix A = readMatrix(file);
    const size_t n = A.size();
    Vector b = readVector(file, n);
    LUDecomposition luObj(A);

    Vector eps(101);
    Vector results(eps.size());
    for (index_t k = 0; k <= 100; ++k) {
        eps[k] = 0.01*k;
        results[k] = norm1(solveLU(luObj, b));
    }
    auto outFile = writeFile(getOutName(arg));
    printColumns(outFile, eps, results);
}

//---------- main ----------//

int main(int argc, char *argv[]) {
    if (argc <= 1) std::cout << "no files to process" << std::endl;
    for (int i = 1; i < argc; ++i) {
        try { handleArg(argv[i]); }
        catch (std::exception &e) { printError(e.what()); }
    }
}


