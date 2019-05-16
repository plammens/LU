//
// Created by Paolo on 17/05/2019.
//

#include "io.h"

#include <regex>
#include <sstream>
#include <exception>
#include <fstream>
#include <iomanip>
#include "norm.h"
#include "errors.h"
#include "sistema.h"

// Constants
const unsigned PRECISION = 12;

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
}

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
}
