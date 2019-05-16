//
// Created by Paolo on 14/05/2019.
//

#ifndef LU_ERRORS_H
#define LU_ERRORS_H

#include <exception>
#include <sstream>


/// Abstract base exception class
class BaseException : public std::exception {
protected:
    std::string message;  ///< message to be returned by what()

    BaseException() = default;
    explicit BaseException(const char *baseMessage, const std::string &details);

public:
    /// `std::exception`-compliant message getter
    const char *what() const noexcept override { return message.c_str(); }
};


/// Indicates that an algorithm encountered a singular matrix
class SingularMatrixError : public BaseException {
public:
    static constexpr auto baseMessage = "singular matrix";
    explicit SingularMatrixError(double tol, const std::string &details = "");

private:
    static std::string getDetails(const std::string &details, double tol);
};

class ValueError : public BaseException {
public:
    static constexpr auto base = "invalid value";
    explicit ValueError(const char *specific = nullptr) : BaseException(base, specific) {}
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


#endif //LU_ERRORS_H
