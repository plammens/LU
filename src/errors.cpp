//
// Created by Paolo on 14/05/2019.
//

#include "errors.h"

BaseException::BaseException(const char *baseMessage, const std::string &details) {
    std::ostringstream oss(baseMessage, std::ios::ate);
    if (not details.empty()) oss << " (" << details << ')';
    message = oss.str();
}

SingularMatrixError::SingularMatrixError(double tol, const std::string &details)
        : BaseException(baseMessage, getDetails(details, tol)) {}

std::string SingularMatrixError::getDetails(const std::string &details, double tol) {
    std::ostringstream oss;
    if (not details.empty()) oss << details << ", ";
    oss.precision(2);
    oss << "tol = " << std::scientific << tol;
    return oss.str();
}
