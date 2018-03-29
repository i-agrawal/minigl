#ifndef __mat__
#define __mat__

#include <cmath>
#include <iostream>

template<class T, int m, int n> struct mat;
template<class T, int m, int n> T dot(const mat<T,m,n>& u,const mat<T,m,n>& v);

//Class mat: Aggregate data structure, initial values may not be 0;
//column major ordered matrix. e.g. A 3x4 matrix, M, will be stored in values as {M00,M10,M20,M01,M11,M21,M01,M12,M22,M02,M12,M22}

template<class T, int m, int n>
struct mat
{
    T values[m*n];

    int rows()
    {return m;}
    
    int cols()
    {return n;}

    void make_zero()
    {for(int i = 0; i < m*n; i++) values[i]=0;}

    mat& operator += (const mat& v)
    {for(int i = 0; i < m*n; i++) values[i] += v.values[i]; return *this;}

    mat& operator -= (const mat& v)
    {for(int i = 0; i < m*n; i++) values[i] -= v.values[i]; return *this;}

    mat& operator *= (const T& c)
    {for(int i = 0; i < m*n; i++) values[i] *= c; return *this;}

    mat& operator /= (const T& c)
    {for(int i = 0; i < m*n; i++) values[i] /= c; return *this;}

    mat operator + () const
    {return *this;}

    mat operator - () const
    {mat r; for(int i = 0; i < m*n; i++) r.values[i] = -values[i]; return r;}

    mat operator + (const mat& v) const
    {mat r; for(int i = 0; i < m*n; i++) r.values[i] = values[i] + v.values[i]; return r;}

    mat operator - (const mat& v) const
    {mat r; for(int i = 0; i < m*n; i++) r.values[i] = values[i] - v.values[i]; return r;}

    mat operator * (const T& c) const
    {mat r; for(int i = 0; i < m*n; i++) r.values[i] = values[i] * c; return r;}

    mat operator / (const T& c) const
    {mat r; for(int i = 0; i < m*n; i++) r.values[i] = values[i] / c; return r;}

    template<int p>
    mat<T,m,p> operator * (const mat<T,n,p>& M) const
    {
        mat<T,m,p> r; r.make_zero();
        for(int i = 0; i < m; i++)
            for(int j = 0; j < n; j++)
                for(int k = 0; k < p; k++)
                    r(i,k) += this->operator()(i,j) * M(j,k);
        return r;
    }

    inline const T& operator() (int i, int j) const
    {return values[i+j*m];}

    inline T& operator() (int i, int j)
    {return values[i+j*m];}
};

template <class T, int m, int n>
mat<T,m,n> operator * (const T& c, const mat<T,m,n>& v)
{return v*c;}

template <class T, int m, int n>
std::ostream& operator << (std::ostream& out, const mat<T,m,n> & M)
{
    for(int i = 0; i < m; i++)
    {
        if(i) out << " ; ";
        for(int j = 0; j < n; j++)
        {
            if(j) out << ' ';
            out << M(i,j);
        }
    }
    return out;
}

#endif
