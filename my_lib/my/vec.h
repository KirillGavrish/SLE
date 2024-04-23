#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numbers>

template<typename T>
using vec = std::vector<T>;

template<typename T>
T dot(vec<T> const &a, vec<T> const &b)
{
    T res = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
        res += a[i] * b[i];
    return res;
}

template<typename T>
vec<T> operator+(vec<T> const &a, vec<T> const &b)
{
    vec<T> res;
    for(std::size_t i = 0; i < a.size(); ++i)
        res.push_back(a[i] + b[i]);
    return res;
}

template<typename T>
vec<T> operator-(vec<T> const &a, vec<T> const &b)
{
    vec<T> res;
    for(std::size_t i = 0; i < a.size(); ++i)
        res.push_back(a[i] - b[i]);
    return res;
}

template<typename T>
vec<T> operator*(vec<T> const &v, T n)
{
	vec<T> res;
	for(std::size_t i = 0; i < v.size(); i++){
		res.push_back(v[i] * n);
	}
	return res;
}

template<typename T>
vec<T> operator/(vec<T> const &v, T n) {return v * (1/n);}

template<typename T>
vec<T> operator*(T n, vec<T> const &v) {return v * n;}

template <typename T>
T abs(T const &a)
{
    if (a > 0) return a;
    else return -a;
}

template <typename T>
T max(vec<T> const &v)
{
    T max = abs(v[0]);
    for (std::size_t i = 0; i < v.size(); ++i)
        if (max < abs(v[i])) max = abs(v[i]);
    return max;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, vec<T> const &vec) 
{
    for (std::size_t i = 0; i < vec.size(); i++)
        os << vec[i] << " ";
    std::cout << std::endl;
    return os;
}

template <typename T>
T norm(vec<T> const &v) {return sqrt(dot(v, v));}

template <typename T>
vec<T> givens_rots(vec<T> const &h, vec<std::pair<T, T>> const &rots)
{
    vec<T> res = h; 
    vec<T> resi = vec<T>(2);
    for (std::size_t i = 0; i < h.size() - 1; ++i)
    {
        resi[0] = rots[i].first * res[i] - rots[i].second * res[i+1];
        resi[1] = rots[i].second * res[i] + rots[i].first * res[i+1];
        res[i] = resi[0];
        res[i+1] = resi[1];
    }
    return res;
}

