#pragma once
#include <map>
#include <complex>
#include <iostream>
#include <cstdint>

namespace AbAl
{
    using POperator = std::map<int, std::uint8_t>;
    using GOperator = std::map<POperator, std::complex<double> >;

    GOperator operator*(const POperator&, const POperator&);
    GOperator operator*(const std::complex<double>&, const POperator&);
    GOperator operator*(const POperator&, const std::complex<double>&);
    GOperator operator*(const double&, const POperator&);
    GOperator operator*(const POperator&, const double&);
    GOperator operator*(const int&, const POperator&);
    GOperator operator*(const POperator&, const int&);
    GOperator operator/(const POperator&, const std::complex<double>&);
    GOperator operator/(const POperator&, const double&);
    GOperator operator/(const POperator&, const int&);

    GOperator operator*(const GOperator&, const GOperator&);
    GOperator operator*(const std::complex<double>&, const GOperator&);
    GOperator operator*(const GOperator&, const std::complex<double>&);
    GOperator operator*(const double&, const GOperator&);
    GOperator operator*(const GOperator&, const double&);
    GOperator operator*(const int&, const GOperator&);
    GOperator operator*(const GOperator&, const int&);
    GOperator operator/(const GOperator&, const std::complex<double>&);
    GOperator operator/(const GOperator&, const double&);
    GOperator operator/(const GOperator&, const int&);

    GOperator operator*(const GOperator&, const POperator&);
    GOperator operator*(const POperator&, const GOperator&);

    GOperator& operator*=(GOperator&, const POperator&);
    GOperator& operator*=(GOperator&, const GOperator&);
    GOperator& operator*=(GOperator&, const std::complex<double>&);
    GOperator& operator*=(GOperator&, const double&);
    GOperator& operator*=(GOperator&, const int&);

    GOperator operator-(const POperator&);
    GOperator operator-(const GOperator&);

    GOperator operator+(const POperator&, const POperator&);
    GOperator operator+(const std::complex<double>&, const POperator&);
    GOperator operator+(const POperator&, const std::complex<double>&);
    GOperator operator+(const double&, const POperator&);
    GOperator operator+(const POperator&, const double&);
    GOperator operator+(const int&, const POperator&);
    GOperator operator+(const POperator&, const int&);
    GOperator operator-(const POperator&, const POperator&);
    GOperator operator-(const std::complex<double>&, const POperator&);
    GOperator operator-(const POperator&, const std::complex<double>&);
    GOperator operator-(const double&, const POperator&);
    GOperator operator-(const POperator&, const double&);
    GOperator operator-(const int&, const POperator&);
    GOperator operator-(const POperator&, const int&);

    GOperator operator+(const GOperator&, const GOperator&);
    GOperator operator+(const std::complex<double>&, const GOperator&);
    GOperator operator+(const GOperator&, const std::complex<double>&);
    GOperator operator+(const double&, const GOperator&);
    GOperator operator+(const GOperator&, const double&);
    GOperator operator+(const int&, const GOperator&);
    GOperator operator+(const GOperator&, const int&);
    GOperator operator-(const GOperator&, const GOperator&);
    GOperator operator-(const std::complex<double>&, const GOperator&);
    GOperator operator-(const GOperator&, const std::complex<double>&);
    GOperator operator-(const double&, const GOperator&);
    GOperator operator-(const GOperator&, const double&);
    GOperator operator-(const int&, const GOperator&);
    GOperator operator-(const GOperator&, const int&);

    GOperator operator+(const POperator&, const GOperator&);
    GOperator operator+(const GOperator&, const POperator&);
    GOperator operator-(const POperator&, const GOperator&);
    GOperator operator-(const GOperator&, const POperator&);

    GOperator& operator+=(GOperator&, const POperator&);
    GOperator& operator+=(GOperator&, const GOperator&);
    GOperator& operator+=(GOperator&, const std::complex<double>&);
    GOperator& operator+=(GOperator&, const double&);
    GOperator& operator+=(GOperator&, const int&);

    GOperator& operator-=(GOperator&, const POperator&);
    GOperator& operator-=(GOperator&, const GOperator&);
    GOperator& operator-=(GOperator&, const std::complex<double>&);
    GOperator& operator-=(GOperator&, const double&);
    GOperator& operator-=(GOperator&, const int&);

    std::ostream& operator<<(std::ostream& os, const POperator& dt);
    std::ostream& operator<<(std::ostream& os, const GOperator& dt);

    GOperator commute(const GOperator&, const GOperator&);
    GOperator commute(const GOperator&, const POperator&);
    GOperator commute(const POperator&, const GOperator&);
    GOperator commute(const POperator&, const POperator&);

    POperator a(const unsigned int&);
    POperator ak(const unsigned int&);

}
