#ifndef XMATH_POLYNOMIAL_H
#define XMATH_POLYNOMIAL_H
#include <utility>
#include <memory.h>
#include "Numerical.h"
template<typename T>
class Polynomial
{
    private:
        unsigned int n = 0;
        T* a = nullptr;
    public:
        Polynomial()
        {
            static_assert(std::is_arithmetic<T>(), "Error: polynomial coefficients must be arithmetic type");
            a = new T;
            *a = static_cast<T>(0);
        }
        Polynomial(unsigned int _n, const T* F) : n(_n)
        {
            static_assert(std::is_arithmetic<T>(), "Error: polynomial coefficients must be arithmetic type");
            a = new T[n];
            for(int i = 0; i <= n; i++)
                a[i] = F[i];
        }
        Polynomial(const Polynomial& F) : n(F.n)
        {
            delete[] a;
            a = new T[n];
            for (int i = 0; i <= n; i++)
                a[i] = F[i];
        }
        Polynomial(Polynomial&& F) noexcept
        {
            n = std::move(F.n);
            delete[] a;
            a = std::move(F.a);
        }
        Polynomial& operator=(const Polynomial& F)
        {
            if(&F == this) return *this;
            else
            {
                n = F.n;
                delete[] a;
                a = new T[n];
                for (int i = 0; i <= n; i++)
                    a[i] = F[i];
            }
        }
        Polynomial& operator=(Polynomial&& F) noexcept
        {
            n = std::move(F.n);
            delete[] a;
            a = std::move(F.a);
        }
        T& operator[](unsigned int i) const { return a[i]; }
        T operator()(const T& x)
        {
            if(!n) return a[0];
            T ret = a[0], pw = x;
            for(int i = 1; i <= n; i++)
            {
                ret += pw * a[i];
                pw += pw * x;
            }
            return ret;
        }
        int deg()
        {
            if(!n)
            {
                if(a[0] == 0) return -1;
                else return 0;
            }
            return n;
        }
        ~Polynomial() { delete[] a; }
};

Polynomial<Real> conv(std::initializer_list<Polynomial<Real> >);
#endif