//
// Created by Paolo on 17/05/2019.
//

#ifndef LU_IO_H
#define LU_IO_H

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "aliases.h"
#include "Matrix.h"
#include "Vector.h"


const unsigned PRECISION = 12;

std::string solveFile(const char *filePath);

template<typename Integer>
unsigned noDigits(Integer num) {
    unsigned count = 1;
    for (num /= 10; num > 0; ++count) num /= 10;
    return count;
}

template<typename Vec>
void printVector(const Vec &vec, std::ostream &os, const char *title = nullptr) {
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

void printInfoNumber(double infoNum, std::ostream &os, const char *title);

Matrix readMatrix(std::istream &is);

Vector readVector(std::istream &is, size_t n);

std::ifstream openFile(const char *path);

std::ofstream writeFile(const std::string &oname);

std::string getOutName(const char *iname);


template <typename Vec, typename... Vecs>
void printColumns(std::ostream &os, Vec vec, Vecs... vecs) {
    os.setf(std::ios::scientific);
    const size_t n = vec.size();
    for (index_t i = 0; i < n; ++i) {
        for (auto &&v : {vecs...})
            os << v[i] << ' ';
        os << '\n';
    }
    os << std::flush;
}


#endif //LU_IO_H
