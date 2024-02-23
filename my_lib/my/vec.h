#pragma once

#include <vector>

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
vec<T> operator*(T n, vec<T> const &v) {return v * n;}

template <typename T>
T max(vec<T> const &v)
{
    T max = v[0];
    for (std::size_t i = 0; i < v.size(); ++i)
        if (max < v[i]) max = v[i];
    return max;
}
