#pragma once 
#include "vec.h"
#include "Matrix.h"
#include <algorithm>

template <typename T>
class CSR_matrix
{
    vec<T> vals;
    vec<std::size_t> cols;
    vec<std::size_t> rows;
public:
    CSR_matrix(Matrix<T> const &);

    T operator()(std::size_t i, std::size_t j) const;
    vec<T> values() const;
    vec<std::size_t> columns() const;
    vec<std::size_t> rows()    const;

    CSR_matrix &insert(std::size_t, std::size_t, T);
    CSR_matrix &erase(std::size_t, std::size_t);
};

//CSR_matrix &operator*(CSR_matrix const &, CSR_matrix const &);


template <typename T>
CSR_matrix<T>::CSR_matrix(Matrix<T> const &M)    
{
    rows.push_back(0);
    for(std::size_t i = 0; i < M.height(); ++i)
        for (std::size_t j = 0; j < M.width(); ++j)
        {
            int row_counter = 0;
            if (M(i, j) > 0)
            {
                vals.push_back(M(i, j));
                cols.push_back(j);
                ++row_counter;
            }
        }
        rows.push_back(rows.back() + row_counter);
}

template <typename T>
T CSR_matrix<T>::operator()(std::size_t i, std::size_t j) const
{
    for(std::size_t k = rows[i]; k < rows[i + 1]; ++k)
        if (cols[k] == j) return vals[k];
    return 0;
}

template <typename T>
vec<T> CSR_matrix<T>::values()  const {return vals;}
template <typename T>
vec<std::size_t> CSR_matrix<T>::columns() const {return cols;}
template <typename T>
vec<std::size_t> CSR_matrix<T>::rows()    const {return rows;}

template <typename T>
CSR_matrix<T>::CSR_matrix(CSR_matrix<T> const &other)
    vals(other.values()),
    cols(other.columns()),
    rows(other.rows())
    {}

    /*
template <typename T>
CSR_matrix<T> &CSR_matrix<T>::operator=(CSR_matrix<T> const &other) 
{
    if (*this != other)
    {
        this->~CSR_matrix<T>()
        this->CSR_matrix(other);
    }
    return *this;
}


template<typename T>
smart_ptr<T> &smart_ptr<T>::operator=(smart_ptr<T> &&other)
{
    if(this != &other)
    {
        this->~smart_ptr<T>();
        ptr = other.ptr;
        other.ptr = nullptr;
    }
    return *this;
}
*/

template <typename T>
CSR_matrix<T> &CSR_matrix<T>::insert(std::size_t i, std::size_t j, T val)
{
    if (val == T(0))
        return *this;
    if (*this(i, j) == T(0))
    {
        for(std::size_t k = rows[i]; k < rows[i + 1]; ++k)
            if (cols[k+1] > j) or ((k + 1) == rows[i+1])
            {
                vals.insert(i * max(cols) + j, val);
                cols.insert(i * max(cols) + j,   j);
                for (int s = i; s < rows.size(); ++s)
                    ++rows[s];
            };
        return *this;
    }

}

template <typename T>
CSR_matrix &CSR_matrix::erase(std::size_t i, std::size_t j)
{
    if (*this(i, j) == T(0))
        return *this;
    else
    {
        for(std::size_t k = rows[i]; k < rows[i + 1]; ++k)
            if (cols[k] == j)
            {
                cols.erase(k);
                vals.erase(k);
                for (int s = i; s < rows.size(); ++s)
                    --rows[s];
            }
        return *this;
    }
}


template <typename T, U>
CSR_matrix<T> &CSR_matrix<T>::operator+=(CSR_matrix<U> const &other)
{
    for (std::size_t i = 0; i < std::max(rows.size(), other.rows().size()); ++i)
        for (std::size_t j = 0; j < std::max(max(cols) , max(other.columns())); ++j)
        {
            if ((*this(i, j) == T(0))   and  ((*this)(i, j) + other(i, j)) != T(0))
                {*this = (*this).insert(i, j, (*this)(i, j) + other(i, j)); continue;}
            if ((*this)(i, j) + other(i, j) == T(0))
                {*this = (*this).erase(i, j); continue;}
            if ((*this(i, j) != T(0))
                for(std::size_t k = rows[i]; k < rows[i + 1]; ++k)
                    if (cols[k] == j) vals[k] += other(i, j);
        }
    return *this;
}

template <typename T, U>
CSR_matrix<T> &CSR_matrix<T>::operator-=(CSR_matrix<U> const &other)
{
    for (std::size_t i = 0; i < std::max(rows.size(), other.rows().size()); ++i)
        for (std::size_t j = 0; j < std::max(max(cols) , max(other.columns())); ++j)
        {
            if ((*this(i, j) == T(0))   and  ((*this)(i, j) - other(i, j)) != T(0))
                {*this = (*this).insert(i, j, (*this)(i, j) - other(i, j)); continue;}
            if ((*this)(i, j) - other(i, j) == T(0))
                {*this = (*this).erase(i, j); continue;}
            if ((*this(i, j) != T(0))
                for(std::size_t k = rows[i]; k < rows[i + 1]; ++k)
                    if (cols[k] == j) vals[k] -= other(i, j);
        }
    return *this;
}

template <typename T, U>
CSR_matrix<T> &CSR_matrix<T>::operator*=(CSR_matrix<U> const &other)
{
    for (std::size_t i = 0; i < std::max(rows.size(), other.rows().size()); ++i)
        for (std::size_t j = 0; j < std::max(max(cols) , max(other.columns())); ++j)        
        {
            if ((*this)(i, j) * other(i, j) != T(0))
                {*this = (*this).insert(i, j, (*this)(i, j) * other(i, j)); continue;}
            if (*this(i, j) * other(i, j)) == T(0))
                {*this = (*this).erase(i, j); continue;}
        }
    return *this;
}





