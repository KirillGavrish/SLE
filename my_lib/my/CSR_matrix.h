#pragma once 
#include "Matrix.h"

template <typename T>
class CSR_matrix
{
    vec<T> vals;
    vec<std::size_t> cols;
    vec<std::size_t> rows;
public:
    CSR_matrix();
    CSR_matrix(Matrix<T> const &);

    vec<T>           get_vals() const;
    vec<std::size_t> get_cols() const;
    vec<std::size_t> get_rows() const;
    
    T lambda_max() const;

    vec<T> operator*(vec<T> const &) const;
    T operator()(std::size_t, std::size_t) const;
    
    bool check(vec<T> const &, vec<std::size_t> const &, vec<std::size_t> const &) const;
};


template <typename T>
CSR_matrix<T>::CSR_matrix()
    : vals(0),
      cols(0),
      rows(0)
    {}
 
template <typename T>
CSR_matrix<T>::CSR_matrix(Matrix<T> const &M)    
{
    rows.push_back(0);
    for(std::size_t i = 0; i < M.get_height(); ++i)   
    {
        std::size_t row_counter = 0;
        for (std::size_t j = 0; j < M.get_width(); ++j)
            if (M(i, j) > 0)
            {
                vals.push_back(M(i, j));
                cols.push_back(j);
                ++row_counter;
            }

        rows.push_back(rows.back() + row_counter);
    }
}

template <typename T>
bool CSR_matrix<T>::check(vec<T> const &values, vec<std::size_t> const &columns, vec<std::size_t> const &rws) const
    {return vals == values && cols == columns && rows == rws;} 

template <typename T>
T CSR_matrix<T>::operator()(std::size_t i, std::size_t j) const
{
    for(std::size_t k = rows[i]; k < rows[i + 1]; ++k)
        if (cols[k] == j) return vals[k];
    return 0;
}

template <typename T> vec<T>           CSR_matrix<T>::get_vals() const {return vals;}
template <typename T> vec<std::size_t> CSR_matrix<T>::get_cols() const {return cols;}
template <typename T> vec<std::size_t> CSR_matrix<T>::get_rows() const {return rows;}

template<typename T>
vec<T> CSR_matrix<T>::operator*(vec<T> const &v) const
{
	vec<T> res = vec<T>(rows.size() - 1);
	for(std::size_t row = 0; row < rows.size() - 1; ++row)
    {
		T sum = 0;
		for(std::size_t i = rows[row]; i < rows[row + 1]; ++i)
			sum += vals[i] * v[cols[i]];
		res[row] = sum;
	}
	return res;
}

template <typename T>
T CSR_matrix<T>::lambda_max() const
{
    T tol = 1e-20;
    vec<T> r(max(cols));
    for (std::size_t i = 0; i < r.size(); ++i)
        r[i] = 1;
    
    T mu = 0, prev;
    for (std::size_t it = 0; it < 10000; ++it)
    {
        prev = mu;
        r = ((*(this)) * r) / norm((*(this)) * r);
        mu = dot(r, (*(this)) * r) / dot(r , r);
        if (mu - prev < tol) break;
    }

    return mu;

}

template <typename T>
vec<T> Simple_Iteration_Method(CSR_matrix<T> const &A, vec<T> const &b, vec<T> const &x0, T const &tol, std::size_t const &Nmax)
{
    T tau = 1 / A.lambda_max();
    vec<T> x = x0;
    vec<T> r; 
    for (std::size_t it = 0; it < Nmax; ++it)
    {
        r = (A * x) - b;
        if (norm(r) < tol) break;

        x = x - tau * r;
    }
    return x;
}

template <typename T>
vec<T> Jacobi_Method(CSR_matrix<T> const &A, vec<T> const &b, vec<T> const &x0, T const &tol, std::size_t const &Nmax)
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

template <typename T>
vec<T> Gauss_Zejdel_Method(CSR_matrix<T> const &A, vec<T> const &b, vec<T> const &x0, T const &tol, std::size_t const &Nmax)
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
vec<T> SIM_Chebyshev(CSR_matrix<T> const &A, vec<T> const &b, vec<T> const &x0, T const &tol, std::size_t const &Nmax, std::size_t const &r, T const &lambda_min, T const &lambda_max)
{
    vec<T> tau(std::size_t(std::pow(2,r)));
    std::size_t n = tau.size();
    T xi;
    for (std::size_t i = 0; i < n; ++i)
    {
        xi = (lambda_max + lambda_min)/2 + (lambda_max - lambda_min)/2 * std::cos((2*i + 1) / (2*n) * std::numbers::pi); tau[i] = 1 / xi;
    }

    vec<std::size_t> permut(n); permut[0] = 0; permut[1] = 1;
    vec<std::size_t> permut_new(n);
    for (std::size_t numbers = 2; numbers < n; numbers *= 2)
    {
        for (std::size_t i = 0; i < numbers; ++i)
        {
            permut_new[i*2] = permut[i];
            permut_new[i*2 + 1] = numbers * 2 - 1 - permut[i];
        }
        permut = permut_new;
    }
       
    vec<T> x = x0;
    vec<T> rx = (A * x) - b;
    for (std::size_t it = 0; it < Nmax; ++it)
    {
        if (norm(rx) < tol) break;

        for (std::size_t const &i : permut)
        {
            x = x - tau[i] * rx;
            rx = (A * x) - b;
        }
    }
    return x;
}

template <typename T>
vec<T> Steepest_Descent(CSR_matrix<T> const &A, vec<T> const &b, vec<T> const &x0, T const &tol, std::size_t const &Nmax)
{
    T tau;
    vec<T> x = x0;
    vec<T> r;
    for (std::size_t it = 0; it < Nmax; ++it)
    {
        r = (A * x) - b;
        if (norm(r) < tol) break;
        
        tau = dot(r, r) / dot(r, A * r);
        x = x - tau * r;
    }
    return x;
}

template <typename T>
vec<T> Sym_Gauss_Zejdel_Method(CSR_matrix<T> const &A, vec<T> const &b, vec<T> const &x0, T const &tol, std::size_t const &Nmax)
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

        for (std::size_t i = x.size() - 1; i > 0; --i)
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


/*
template <typename T>
std::size_t  Gauss_Zejdel_Method_kr(CSR_matrix<T> const &A, vec<T> const &b, vec<T> const &x0, T const &tol, std::size_t const &Nmax)
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

    vec<T> x1 = x0;
    
    std::size_t it = 0;
    while(norm(x1 - x) > (norm(x)/1000))
    {
        for (std::size_t i = 0; i < x1.size(); ++i)
        {
            x1[i] = b[i];
            for (std::size_t k = A.get_rows()[i]; k < A.get_rows()[i+1]; ++k)
            {
                if (i != A.get_cols()[k])
                    x1[i] -= A.get_vals()[k] * x1[A.get_cols()[k]];
                else
                    d = A.get_vals()[k];
            }
            x1[i] /= d;
        }
        it++;
    }


    return it;
}


template <typename T>
auto CSR_matrix<T>::insert(std::size_t i, std::size_t j, T val)
{
    if (val == T(0))
        return vals.begin();
    if ((*this)(i, j) == T(0))
    {
        for(std::size_t k = rows[i]; k < rows[i + 1]; ++k)
            if ((cols[k+1] > j) or ((k + 1) == rows[i+1]))
            {
                std::size_t counter = 0;
                for (auto it = vals.begin(); it != vals.end(); ++it, ++counter)
                    if (counter == i * max(cols) + j)
                        {vals.insert(it, val); break;}

                counter = 0;
                for (auto it = vals.begin(); it != vals.end(); ++it, ++counter)
                    if (counter == i * max(cols) + j)
                        {cols.insert(it,   j); break;}

                for (std::size_t s = i; s < rows.size(); ++s)
                    ++rows[s];
            }
        return vals.begin();
    }
}

template <typename T>
auto CSR_matrix<T>::erase(std::size_t i, std::size_t j)
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
                for (std::size_t s = i; s < rows.size(); ++s)
                    --rows[s];
            }
        return *this;
    }
}


template <typename T>
CSR_matrix<T> &CSR_matrix<T>::operator+=(CSR_matrix<T> const &other)
{
    for (std::size_t i = 0; i < std::max(rows.size(), other.rows().size()); ++i)
        for (std::size_t j = 0; j < std::max(max(cols) , max(other.columns())); ++j)
        {
            if ((*this(i, j) == T(0))   and  ((*this)(i, j) + other(i, j)) != T(0))
                {*this = (*this).insert(i, j, (*this)(i, j) + other(i, j)); continue;}
            if ((*this)(i, j) + other(i, j) == T(0))
                {*this = (*this).erase(i, j); continue;}
            if ((*this)(i, j) != T(0))
                for (std::size_t k = rows[i]; k < rows[i + 1]; ++k)
                    if (cols[k] == j) vals[k] += other(i, j);
        }
    return *this;
}

template <typename T>
CSR_matrix<T> &CSR_matrix<T>::operator-=(CSR_matrix<T> const &other)
{
    for (std::size_t i = 0; i < std::max(rows.size(), other.rows().size()); ++i)
        for (std::size_t j = 0; j < std::max(max(cols) , max(other.columns())); ++j)
        {
            if ((*this(i, j) == T(0))   and  ((*this)(i, j) - other(i, j)) != T(0))
                {*this = (*this).insert(i, j, (*this)(i, j) - other(i, j)); continue;}
            if ((*this)(i, j) - other(i, j) == T(0))
                {*this = (*this).erase(i, j); continue;}
            if ((*this)(i, j) != T(0))
                for (std::size_t k = rows[i]; k < rows[i + 1]; ++k)
                    if (cols[k] == j) vals[k] -= other(i, j);
        }
    return *this;
}

template <typename T>
CSR_matrix<T> &CSR_matrix<T>::operator*=(CSR_matrix<T> const &other)
{
    for (std::size_t i = 0; i < std::max(rows.size(), other.rows().size()); ++i)
        for (std::size_t j = 0; j < std::max(max(cols) , max(other.columns())); ++j)
        {
            if ((*this)(i, j) * other(i, j) != T(0))
                {*this = (*this).insert(i, j, (*this)(i, j) * other(i, j)); continue;}
            if (((*this)(i, j) * other(i, j)) == T(0))
                {*this = (*this).erase(i, j); continue;}
        }
    return *this;
}
*/
