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
    Matrix(std::size_t const &);

    T operator()(std::size_t, std::size_t) const;
    std::size_t get_width() const;
    std::size_t get_height() const;
    
    vec<T> get_col(std::size_t const &) const;
    vec<T> get_str(std::size_t const &) const;
    vec<T> vals() const;
//    vec<T> pos_vals() const;
    /*
    Matrix<T> &operator*=(Matrix<T> const &);
    Matrix<T> &operator+=(Matrix<T> const &);
    Matrix<T> &operator-=(Matrix<T> const &);
    */
    vec<T>    operator*(vec<T> const &) const;
    Matrix<T> operator*(Matrix<T> const &) const;
    Matrix<T> transpose() const;
};



template <typename T>
T Matrix<T>::operator()(std::size_t i, std::size_t j) const {return elements[i * width + j];}

template <typename T>
vec<T> Matrix<T>::vals() const {return elements;};
/*
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
*/
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
std::size_t Matrix<T>::get_height() const {return width != 0 ? elements.size() / width : 0;}

template <typename T>
std::size_t Matrix<T>::get_width() const {return width;}

template <typename T>
vec<T> Matrix<T>::operator*(vec<T> const &v) const
{
	vec<T> res = vec<T>(get_height());
	for(std::size_t i = 0;  i < res.size(); ++i)
    {
		res[i] = 0;
		for(std::size_t j = 0; j < width; ++j)
			res[i] += elements[i * width + j] * v[j];
	}
	return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(Matrix<T> const &other) const
{
    vec<T> res_vals = vec<T>(get_height() * other.get_width());
    for(std::size_t i = 0;  i < get_height(); ++i)
        for (std::size_t j = 0; j < other.get_width(); ++j)
        {
            res_vals[i * other.get_width() + j] = 0;
            for(std::size_t k = 0; k < width; ++k)
                res_vals[i * other.get_width() + j]+= elements[i * width + k] * other.get_col(j)[k];
        }

    return Matrix<T>(res_vals, other.get_width());
}


template <typename T>
Matrix<T> add_col_H(Matrix<T> &M, vec<T> const &v)
{
    vec<T> res_vals;
    for (std::size_t i = 0; i < M.get_height(); ++i)
    {
        for (std::size_t j = 0; j < M.get_width(); ++j)
            res_vals.push_back(M.vals()[i * M.get_width() + j]);
        res_vals.push_back(v[i]);
    }

    for (std::size_t j = 0; j <= M.get_width(); ++j)
        if (j == M.get_height()) res_vals.push_back(v[j]);
        else res_vals.push_back(0);
    
    return Matrix<T>(res_vals, M.get_width()+1);
}

template <typename T>
Matrix<T> add_col(Matrix<T> &M, vec<T> const &v)
{
    vec<T> res_vals;
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        for (std::size_t j = 0; j < M.get_width(); ++j)
            res_vals.push_back(M.vals()[i * M.get_width() + j]);
        res_vals.push_back(v[i]);
    }

    return Matrix<T>(res_vals, M.get_width()+1);
}

template <typename T>
Matrix<T> Em(std::size_t const &m) 
{
    vec<T> elements(m*m);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j)
        {
            if (i == j) elements[i * m + j] = 1;
            else        elements[i * m + j] = 0;
        }
    return Matrix<T>(elements, m);
}

template <typename T>
Matrix<T> Matrix<T>::transpose() const
{
    vec<T> elementsT;
     
	for(std::size_t j = 0; j < width; ++j)
		for(std::size_t i = 0; i < get_height(); ++i)
            elementsT.push_back(elements[i * width + j]);
    return Matrix<T>(elementsT, get_height());
}

template <typename T>
vec<T> Matrix<T>::get_col(std::size_t const &n) const
{
    vec<T> res = vec<T>(get_height());
    for (std::size_t i = 0; i < res.size(); ++i)
        res[i] = elements[i * width + n];
    return res;
}   

template <typename T>
vec<T> Matrix<T>::get_str(std::size_t const &n) const
{
    vec<T> res = vec<T>(width);
    for (std::size_t j = 0; j < width; ++j)
        res[j] = elements[n * width + j];
    return res;
}

template <typename T>
Matrix<T> theta(vec<T> const &nu, Matrix<T> const &M)
{
	vec<T> res = vec<T>(M.get_height() * M.get_width());
    for(std::size_t j = 0; j < M.get_width(); ++j)
    {
		vec<T> x = M.get_col(j);
		vec<T> theta = x - 2 * dot(x, nu) / dot(nu, nu) * nu;
		for(std::size_t i = 0; i < M.get_height(); ++i)
			res[i * M.get_height() + j] = theta[i];
	}
	return Matrix<T>(res, M.get_width());
}

template <typename T>
std::pair<Matrix<T>, Matrix<T>> QR_decomp(Matrix<T> const &M)
{
    Matrix<T> E = Em<T>(M.get_width());

	Matrix<T> Q = E;
	Matrix<T> R = M;
	for(std::size_t j = 0; j < R.get_width(); ++j){
		vec<T> nu = vec<T>(R.get_width());
		for(std::size_t i = 0; i < R.get_width(); ++i)
        {
			if(i >= j) nu[i] = R(i, j);
			else       nu[i] = 0;
		}
		
        nu[j] += nu[j] > 0 ? norm(nu) : -norm(nu);
		Q = theta(nu, Q.transpose()).transpose();
		R = theta(nu, R);
	}
	return std::pair<Matrix<T>, Matrix<T>>(Q, R);
}

template <typename T>
vec<T> Inverse_Gauss_Method(Matrix<T> const &A, vec<T> const &f)
{
    vec<vec<T>> rows;
    vec<T> b = f;
    vec<T> x;
    T alpha;
    for (std::size_t i = 0; i < f.size(); ++i)
    {
        rows.push_back(A.get_str(i));
        for (std::size_t k = 0; k < i; ++k)
        {
            for (std::size_t z = 1; rows[i][i] == 0; ++z)
            { 
                rows[i] = rows[i] + rows[i-z];
                b[i] = b[i] + b[i-z];
            }
            alpha = rows[k][i] / rows[i][i];
            rows[k] = rows[k] - rows[i] * alpha;
            b[k] = b[k] - b[i] * alpha;
        }
    }
    
    for (std::size_t k = 0; k < f.size(); ++k)
            x.push_back(b[k] / rows[k][k]);
    return x;
}
