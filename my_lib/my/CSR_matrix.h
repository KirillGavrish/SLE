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
    std::size_t n = x.size(), I;

    for (std::size_t it = 0; it < Nmax; ++it)
    {
        for (std::size_t i = 0; i < n; ++i)
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

        for (std::size_t i = 0; i < n; ++i)
        {
            I = n - i - 1;
            x[I] = b[I];
            for (std::size_t k = A.get_rows()[I]; k < A.get_rows()[I+1]; ++k)
            {
                if (I != A.get_cols()[k])
                    x[I] -= A.get_vals()[k] * x[A.get_cols()[k]];
                else
                    d = A.get_vals()[k];
            }
            x[I] /= d;
        }

        if (max(A*x - b) < tol) break;
    }

    return x;
}

template <typename T, typename Method> 
vec<T> Chebyshev2(CSR_matrix<T> const &A, vec<T> const &b, vec<T> const &x0, T const &tol, std::size_t const &Nmax, T const &rho = 0.9, Method method = Simple_Iteration_Method)
{
    vec<T> x = x0;
    x = method(A, b, x, tol, 1);
    T mu_prev = 1, mu = rho, mu_tmp;
    for (std::size_t it = 0; it < Nmax; ++it)
    {
        mu_tmp = mu;
        mu = 2 / rho * mu - mu_prev;
        
        x = 2 * mu_tmp / (rho * mu) * method(A, b, x, tol, 1) - mu_prev / mu * x;

        if (norm(A*x - b) < tol) break;
        mu_prev = mu_tmp;
    }
    return x;
}

template <typename T>
vec<T> Conjugate_Gradient_Method(CSR_matrix<T> const &A, vec<T> const &b, vec<T> const &x0, T const &tol)
{
    T alpha, beta = 0;
    vec<T> x = x0;
    vec<T> r = A * x - b, r_prev;
    vec<T> d = r;
    for (std::size_t it = 0; it < x.size(); ++it)                       // вычислительная сложность O(n^2)
    {                                                                   // по памяти О(n)
        if (norm(r) < tol) break;                                       
                                                                        
        alpha = dot(r, r) / dot(d, A * d);                              
        x = x - alpha * d;                                              
        
        r_prev = r;
        r = A * x - b;
        
        beta = dot(r, r) / dot(r_prev, r_prev);
        d = r + beta * d;
    }
    return x;
}

template <typename T>
vec<T> GMRES_m(CSR_matrix<T> const &A, vec<T> const &b, vec<T> const &x0, T const &tol, std::size_t const &m, std::size_t const &Nmax)
{
    vec<T> r0, x, x0_, vk, h, h_rots, y, by;                                // по памяти O(n*m)
    vec<std::pair<T, T>> giv_rots;                                          // вычислительная сложность O(n*m)
    Matrix<T> V, H;                                                          
    x = x0;                                                                  
    for (std::size_t it = 0; it < Nmax; ++it)                                
    {
        x0_ = x;
        r0 = A * x0_ - b;
        V = Matrix<T>(r0 / norm(r0), 1);
        
        H = Matrix<T>({}, 0);
        giv_rots = vec<std::pair<T, T>>(0);
        
        for (std::size_t k = 0; k < m; ++k)
        {  
            vk = A * V.get_col(k);
            h = vec<T>(0);
            for (std::size_t i = 0; i < k+1; ++i)
            {   
                h.push_back(dot(A * V.get_col(k), V.get_col(i)));
                vk = vk - h[i] * V.get_col(i);
            }
            h.push_back(norm(vk));
            vk = vk / norm(vk);

            h_rots = {h[k], h[k+1]};
            giv_rots.push_back(std::pair(h_rots[0]/norm(h_rots), (-1)*h_rots[1]/norm(h_rots)));
            
            h = givens_rots(h, giv_rots);
            h.pop_back();
            H = add_col_H(H, h);

            std::reverse(giv_rots.begin(), giv_rots.end());
            
            by = vec<T>(0);
            by.push_back(norm(r0));
            for (std::size_t i = 0; i < k ; ++i)
                by.push_back(0);

            by = givens_rots(by, giv_rots);
            
            by = (-1.) * by;
            y = Inverse_Gauss_Method(H, by);

            x = x0_ + V * y;
            if (norm(A*x - b) < tol) break;

            V = add_col(V, vk);
        }
    }
    return x;
}

