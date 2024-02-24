#pragma once
#include <iostream>
#include "vec.h"

template <typename T>
class Matrix
{
    vec<T> elements;
    std::size_t width = 0;

public:
    
    Matrix(std::size_t m = 0, std::size_t n = 0);
    Matrix(vec<T> const &, std::size_t width = 0);


    T operator()(std::size_t, std::size_t) const;
    std::size_t get_width() const;
    std::size_t get_height() const;

    vec<T> vals() const;
    vec<T> pos_vals() const;
    
    Matrix<T> &operator*=(Matrix<T> const &);
    Matrix<T> &operator+=(Matrix<T> const &);
    Matrix<T> &operator-=(Matrix<T> const &);
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
            res.push_back(values[i]);
    return res;
}

template <typename T>
Matrix<T>::Matrix(std::size_t m, std::size_t n)
    : elements(m*n),
      width(n)
    {}

template <typename T>
Matrix<T>::Matrix(vec<T> const &vals, std::size_t width)
    : elements(vals),
      width(width)
    {}

template <typename T>
std::size_t Matrix<T>::get_height() const {return elements.size() / width;}

template <typename T>
std::size_t Matrix<T>::get_width() const {return width;}

template <typename T>
vec<T> operator*(Matrix<T> const &M, vec<T> const &v)
{
    vec<T> res;
    for(std::size_t i = 0; i < M.get_height(); ++i)
    {
        T res_j = 0;
        for(std::size_t k = 0; k < M.get_width(); ++k)
            res_j += M(i, k) * v[k];
        res.push_back(res_j);
    }
    return res;
}

template <typename T>
Matrix<T> &Matrix<T>::operator*=(Matrix<T> const &other)
{
    vec<T> res;
    for(std::size_t i = 0; i < get_height(); ++i)
        for (std::size_t j = 0; j < other.get_width(); ++j)
        {
            T res_ij = 0.;
            for(std::size_t k = 0; k < width; ++k)
                res_ij += (*this)(i, k) * other(k, j);
            res.push_back(res_ij);
        }
    *this = {res, other.get_width()};
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator+=(Matrix<T> const &other)
{
    elements = elements + other.vals();
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator-=(Matrix<T> const &other) 
{
    elements = elements - other.vals();
    return *this;
}

template <typename T>
Matrix<T> operator*(Matrix<T> const &a, Matrix<T> const &b) {Matrix<T> c = a; return c *= b;}
template <typename T>
Matrix<T> operator-(Matrix<T> const &a, Matrix<T> const &b) {Matrix<T> c = a; return c -= b;}
template <typename T>
Matrix<T> operator+(Matrix<T> const &a, Matrix<T> const &b) {Matrix<T> c = a; return c += b;}
