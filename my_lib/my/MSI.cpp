#pragma once
#include "CSR_matrix.h"

template <typename T>
vec<T> MSI(CSR_matrix<T> const &A, vec<T> const &b, vec<T> const &x0, T const &tol std::size_t const &Nmax)
{
    T tau = 2 / lambda_max(A);
    vec<T> r(A.width());
    vec<T> x(x0);

    for (std::size_t i = 0; i < Nmax; ++i)
    {
        r = tau * (A * x - b);
        x = x - r;
        if (max(r) < tol)
            break;
    }

    return x;
}


template <typename T>
vec<T> Gauss_Zejdel(CSR_matrix<T> const &A, vec<T> const &b, vec<T> const &x0, T const &tol, std::size_t const &Nmax)
{
    vec<T> x = x0;
    T d;

    for (std::size_t it = 0; it < Nmax; ++it)
    {
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            x[i] = b[i];
            for (std::size_t k = A.get_rows()[i]; k < A.get_rows()[i+1]; ++k)
            {
                if (i != A.get_cols()[k])
                    x[i] -= A.get_vals()[k] * x[A.get_cols()[k]]; 
                else
                    d = A.get_vals()[k];
            }
            x[i] /= d;
        }

        if (max(A*x - b) < tol) break;
    }

    return x;
}


template <typename T>
vec<T> Jacobi_Method(CSR_matrix<T> const &A, vec<T> const &b, vec<T> const &x0, T const &tol std::size_t const &Nmax)
{
    vec<T> x = x0;
    vec<T> res;
    T d;

    for (std::size_t it = 0; it < Nmax; ++it)
    {
        res = b;
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            for (std::size_t k = A.get_rows()[i]; k < A.get_rows()[i+1]; ++k)
            {
                if (i != A.get_cols()[k])
                    res[i] -= A.get_vals()[k] * x[A.get_cols()[k]];
                else
                    d = A.get_vals()[k];
            }
            res[i] /= d;
        }

        x = res;
        if (max(A*x - b) < tol) break;
    }

    return x;
}






