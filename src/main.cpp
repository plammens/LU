#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <exception>
#include <sstream>

#define SISTEMA_DECLARATIONS_ONLY

#include "sistema.cpp"



//---------- Exceptions ----------//


class BaseException : public std::exception {
protected:
    std::string message;  ///< message to be returned by what()

    BaseException() = default;
    explicit BaseException(const char *baseMessage, const char *details) {
        std::ostringstream oss(baseMessage, std::ios::ate);
        if (details) oss << " (" << details << ')';
        message = oss.str();
    }

public:
    /// `std::exception`-compliant message getter
    const char *what() const noexcept override { return message.c_str(); }
};


class IOError : public BaseException {
public:
    static constexpr auto base = "unable to access file";
    explicit IOError(const char *specific = nullptr) : BaseException(base, specific) {}
};


class BadFormat : public BaseException {
public:
    static constexpr auto base = "bad input format";
    explicit BadFormat(const char *specific = nullptr) : BaseException(base, specific) {}
};



//---------- Helper functions ----------//

Matrix readMatrix(std::istream &is) {
    size_t n = 0, m = 0;
    is >> n >> m;
    if (m > n*n)
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

void printResult(const SolveResult &result, std::ostream &os) {
    const Vector &x = result.solution();
    const size_t n = x.size();
    const unsigned width = noDigits(n);

    for (index_t i = 0; i < n; ++i)
        os << std::setw(width) << i << " \t"
           << (x[i] >= 0 ? " ": "") << std::scientific << std::setprecision(9) << x[i]
           << '\n';

    os << '\n' << "residue: " << std::scientific << result.residue() << '\n';
    os.flush();
}

void parseFile(const char *filename) {
    static index_t id = 0;

    std::ifstream inputFile(filename);
    if (not inputFile.is_open()) throw IOError(filename);

    Matrix &&A = readMatrix(inputFile);
    Vector &&b = readVector(inputFile, A.size());
    auto result = solve(A, b);
    if (not result) throw SingularMatrixError();
    else {
        std::ostringstream oss;
        oss << "SOLUTION" << std::setw(2) << std::setfill('0') << id++ << ".DAT";
        std::string outName = oss.str();
        std::ofstream outFile(outName);
        if (not outFile.is_open()) throw IOError(outName.c_str());
        printResult(result, outFile);
    }
};

inline
void printError(const char *message, const char *arg = nullptr) {
    std::cerr << "error: " << message;
    if (arg) std::cerr << ' ' << arg;
    std::cerr << '\n';
}


//---------- main ----------//

int main(int argc, char *argv[]) {
    for (int i = 1; i < argc; ++i) {
        try { parseFile(argv[i]); }
        catch (std::exception &e) { printError(e.what()); }
    }
}


