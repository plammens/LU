#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <exception>
#include <sstream>
#include <regex>

#include "sistema.h"
#include "errors.h"
#include "norm.h"


// Constants
const unsigned PRECISION = 12;


//---------- Helper functions ----------//

Matrix readMatrix(std::istream &is) {
    size_t n = 0, m = 0;
    is >> n >> m;
    if (m > n * n)
        throw BadFormat("number of elements larger than n^2 while reading matrix");

    Matrix mat(n);
    for (; m > 0; --m) {
        index_t i = 0, j = 0;
        is >> i >> j;
        is >> mat(i, j);
    }

    if (is.fail()) throw BadFormat("couldn't parse matrix");
    return mat;
}


Vector readVector(std::istream &is, size_t n) {
    size_t k = 0;
    is >> k;

    Vector vec(0.0, n);
    for (; k > 0; --k) {
        index_t i;
        is >> i;
        is >> vec[i];
    }

    return vec;
}

template<typename Integer>
unsigned noDigits(Integer num) {
    unsigned count = 1;
    for (num /= 10; num > 0; ++count) num /= 10;
    return count;
};


template<typename Vec>
void printVector(const Vec &vec, std::ostream &os, const char *title = "") {
    size_t n = vec.size();
    const unsigned width = noDigits(n);

    if (title and *title != '\0') os << title << ":\n";
    for (index_t i = 0; i < n; ++i)
        os << std::setw(width) << i << " \t"
           << (vec[i] >= 0 ? " " : "")
           << std::scientific << vec[i]
           << '\n';

    os << '\n';
}

void printInfoNumber(double infoNum, std::ostream &os, const char *title) {
    os << title << ": " << std::scientific << infoNum << "\n\n";
}

void printResult(const SolveResult &result, const ExtraSolveInfo &info, std::ostream &os) {
    os.precision(PRECISION);

    printVector(result.solution(), os);
    printInfoNumber(info.residue, os, "residue");
    printInfoNumber(info.cond1, os, "condition number μ_1");
    printInfoNumber(info.condInf, os, "condition number μ_Inf");
    printVector(result.getLU().perm().vector(), os, "permutation vector");

    os.flush();
}

std::string getFileName(const std::string &path) {
    static const std::regex regex(R"(^.*?[/\\]?(\w+)(\.\w+)?$)");
    std::smatch match;
    std::regex_match(path, match, regex);
    if (not match.empty()) return match[1];
    return "";
}


std::string solveFile(const char *filePath) {
    std::ifstream inputFile(filePath);
    if (not inputFile.is_open()) throw IOError(filePath);

    Matrix &&A = readMatrix(inputFile);
    Vector &&b = readVector(inputFile, A.size());
    auto result = solve(A, b);
    if (not result) throw SingularMatrixError(result.tol());
    ExtraSolveInfo info = getExtraSolveInfo(A, b, result);

    std::ostringstream oss;
    oss << "SOLUTION_" << getFileName(filePath) << ".DAT";
    std::string outName = oss.str();
    std::ofstream outFile(outName);
    if (not outFile.is_open()) throw IOError(outName.c_str());
    printResult(result, info, outFile);
    return outName;
};

inline
void printError(const char *message, const char *arg = nullptr) {
    std::cerr << "error: " << message;
    if (arg) std::cerr << ' ' << arg;
    std::cerr << '\n';
}


//---------- main ----------//

int main(int argc, char *argv[]) {
    if (argc <= 1) std::cout << "no files to process" << std::endl;
    for (int i = 1; i < argc; ++i) {
        try {
            auto iName = argv[i];
            auto oName = solveFile(iName);
            std::cout << "written solution of " << iName
                      << " to " << oName << std::endl;
        }
        catch (std::exception &e) { printError(e.what()); }
    }
}


