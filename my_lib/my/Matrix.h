#pragma once

#include "vec.h"

template <typename T>
class Matrix
{
    vec elements;
    std::size_t width = 0;

public:

    Matrix(std::size_t m = 0, std::size_t n = 0);

    Matrix(vec const &);

    double operator()(std::size_t, std::size_t) const;
    std::size_t width();
    std::size_t height();

    vec<T> vals() const;
    vec<T> pos_vals() const;
};



template <typename T>
T Matrix<T>::operator()(std::size_t i, std::size_t j) const {return elements[i * width + j];}

template <typename T>
vec<T> Matrix<T>::vals() const {return elements;};

template <typename T>
vec<T> Matrix<T>::pos_vals() const
{
    vec<T> res;
    vec<T> values = vals();
    for (std::size_t i = 0; i < values.size(); ++i)
        if (values[i] > 0)
            res.push_back(values[i])
    return res;
}

template <typename T>
Matrix<T>::Matrix(std::size_t m, std::size_t n)
    : elements(m*n),
      width(n),
    {}

template <typename T>
Matrix<T>::Matrix(vec<T> vals, std::size_t width)
    : elements(vals),
      width(width),
    {}

template <typename T>
std::size_t Matrix<T>height() const {return elements.size() / width}

template <typename T>
std::size_t Matrix<T>::width() const {return width;}


template <typename T, U>
Matrix<T> &Matrix<T>::operator*=(Matrix<U> const &other)
{
    vec<T> res;
    for(std::size_t i = 0; i < height(); ++i)
        for (std::size_t j = 0; j < other.width(); ++j)
        {
            T res_ij = 0.;
            for(std::size_t k = 0; k < width; ++j)
                res_ij += (*this)(i, k) * other(k, j);
            res.push_back(res_ij);
        }
    *this = Martix(res, other.width());
    return *this;
}

template <typename T, U>
Matrix<T> &Matrix<T>::operator+=(Matrix<U> const &other)
{
    elements = elements + other.vals();
    return *this;
}

template <typename T, U>
Matrix<T> &Matrix<T>::operator-=(Matrix<U> const &other) 
{
    elements = elements - other.vals();
    return *this;
}

template <typename T>
Matrix<T> &operator*(Matrix<T> const &a, Matrix<T> const &b) {Matrix<T> c = a; return c *= b;}
template <typename T>
Matrix<T> &operator-(Matrix<T> const &a, Matrix<T> const &b) {Matrix<T> c = a; return c -= b;}
template <typename T>
Matrix<T> &operator+(Matrix<T> const &a, Matrix<T> const &b) {Matrix<T> c = a; return c += b;}
