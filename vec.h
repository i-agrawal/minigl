#ifndef __vec__
#define __vec__

#include <cmath>
#include <iostream>
#include <cassert>
#include "mat.h"

template<class T, int n> struct vec;
template<class T, int n> T dot(const vec<T,n>& u,const vec<T,n>& v);

template<class T, int n>
struct vec: mat<T,n,1>
{
    typedef mat<T,n,1> BASE;
    using BASE::make_zero; using BASE::values;

    enum{dim=n};
    vec()
    {make_zero();}

    vec(const T& a)
    {static_assert(n == 1,"vector size does not match this constructor");values[0]=a;}

    vec(const T& a, const T& b)
    {static_assert(n == 2,"vector size does not match this constructor");values[0]=a;values[1]=b;}

    vec(const T& a, const T& b, const T& c)
    {static_assert(n == 3,"vector size does not match this constructor");values[0]=a;values[1]=b;values[2]=c;}
    
    vec(const T& a, const T& b, const T& c, const T& d)
    {static_assert(n == 4,"vector size does not match this constructor");values[0]=a;values[1]=b;values[2]=c;values[3]=d;}

    vec(const mat<T,n,1>& M)
    {for(int i=0; i<n; i++) values[i]=M.values[i];}

    int size()
    {return n;}

    const T& operator[] (int i) const
    {return values[i];}

    T& operator[] (int i)
    {return values[i];}

    T magnitude_squared() const
    {return dot(*this, *this);}

    T magnitude() const
    {return sqrt(magnitude_squared());}

    // Be careful to handle the zero vector gracefully
    vec normalized() const
    {T mag = magnitude(); if(mag) return *this / mag; vec r; r[0] = 1; return r;};

    vec<T,n-1> remove_dim(int d) const
    {vec<T,n-1> v; for(int i=0; i<d; i++) v[i]=values[i];
        for(int i=d; i<n; i++) v[i]=values[i+1];
        return v;}

    vec<T,n-1> remove_last_dim() const
    {vec<T,n-1> v; for(int i=0; i<n-1; i++) v[i]=values[i]; return v;}

};

template <class T, int n>
T dot(const vec<T,n> & u, const vec<T,n> & v)
{
    T r  =  0;
    for(int i = 0; i < n; i++) r += u.values[i] * v.values[i];
    return r;
}

template <class T >
vec<T,3> cross(const vec<T,3> & u, const vec<T,3> & v)
{
    return vec<T,3> (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0]);
}

template <class T, int n>
std::ostream& operator << (std::ostream& out, const vec<T,n> & u)
{
    for(int i = 0; i < n; i++)
    {
        if(i) out << ' ';
        out << u[i];
    }
    return out;
}

typedef vec<float,2> vec2;
typedef vec<float,3> vec3;
typedef vec<int,2> ivec2;

#endif
