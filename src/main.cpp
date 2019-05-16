#include <iostream>
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

//---------- main ----------//

int main(int argc, char *argv[]) {
    if (argc <= 1) std::cout << "no files to process" << std::endl;
    for (int i = 1; i < argc; ++i) {
        try { handleArg(argv[i]); }
        catch (std::exception &e) { printError(e.what()); }
    }
}


