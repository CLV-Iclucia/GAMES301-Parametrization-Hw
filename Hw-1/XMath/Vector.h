#ifndef XMATH_VECTOR_H
#define XMATH_VECTOR_H
#include <utility>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <cassert>
template<typename T>
class Vector
{
    private:
        unsigned int n = 0;
        T* v = nullptr;
    public:
        Vector(unsigned int _n, const T* A) : n(_n)
        {
            static_assert(std::is_arithmetic<T>(), "Error: Elements of vector must be arithmetic type!");
            v = new T[n];
            for(int i = 0; i < n; i++)
                v[i] = A[i];
        }
        explicit Vector(unsigned int _n) : n(_n)
        {
            v = new T[n];
            for(int i = 0; i < n; i++)
                v[i] = static_cast<T>(0);
        }
        Vector(const Vector& A)
        {
            assert(n <= A.n);
            n = A.n;
            delete[] v;
            v = new T[n];
            for(int i = 0; i < n; i++)
                v[i] = A[i];
        }
        Vector(Vector&& A) noexcept
        {
            assert(n <= A.n);
            n = std::move(A.n);
            delete[] v;
            v = std::move(A.v);
        }
        Vector& operator=(const Vector& A)
        {
            assert(n <= A.n);
            if(&A == this) return *this;
            else
            {
                n = A.n;
                delete[] v;
                v = new T[n];
                for(int i = 0; i < n; i++)
                    v[i] = A[i];
                return *this;
            }
        }
        Vector& operator=(Vector&& A) noexcept
        {
            assert(n <= A.n);
            n = std::move(A.n);
            delete[] v;
            v = std::move(A.v);
            return *this;
        }
        T dot(const Vector& A) const
        {
            T sum = 0;
            for (int i = 0; i < n; i++)
                sum += this->v[i] * A.v[i];
            return sum;
        }
        T& operator[](int i) { return v[i]; }
        const T& operator[](int i) const { return v[i]; }
        Vector operator*(const Vector& A) const
        {
            Vector V;
            for (int i = 0; i < n; i++)
                V.v[i] = v[i] * A.v[i];
            return V;
        }
        Vector operator+(const Vector& A) const
        {
            Vector V;
            for (int i = 0; i < n; i++)
                V.v[i] = v[i] + A.v[i];
            return V;
        }
        Vector operator-()
        {
            Vector V;
            for (int i = 0; i < n; i++)
                V.v[i] = -v[i];
            return V;
        }
        Vector operator/(const Vector& A) const
        {
            for (int i = 0; i < n; i++)
            {
                try
                {
                    if (A.v[i] == 0)throw std::runtime_error("Division by zero!");
                }
                catch (const std::exception& e)
                {
                    std::cerr << e.what() << '\n';
                    exit(-1);
                }
            }
            Vector V;
            for (int i = 0; i < n; i++)
                V.v[i] = v[i] / A.v[i];
            return V;
        }
        Vector operator-(const Vector& A) const
        {
            Vector V;
            for (int i = 0; i < n; i++)
                V.v[i] = v[i] - A.v[i];
            return V;
        }
        Vector operator/(const T& val) const
        {
            try
            {
                if (val == 0)throw std::runtime_error("Division by zero!");
            }
            catch (const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                exit(-1);
            }
            Vector V;
            for (int i = 0; i < n; i++)
                V.v[i] = v[i] / val;
            return V;
        }
        Vector operator*(T val) const
        {
            Vector V;
            for (int i = 0; i < n; i++)
                V.v[i] = val * v[i];
            return V;
        }
        Vector& operator*=(const Vector& A)
        {
            for (int i = 0; i < n; i++)
                v[i] *= A.v[i];
            return *this;
        }
        Vector& operator+=(const Vector& A)
        {
            for (int i = 0; i < n; i++)
                v[i] += A.v[i];
            return *this;
        }
        Vector& operator/=(const Vector& A)
        {
            for (int i = 0; i < n; i++)
            {
                try
                {
                    if (A.v[i] == 0)throw std::runtime_error("Division by zero!");
                    else v[i] /= A.v[i];
                }
                catch (const std::exception& e)
                {
                    std::cerr << e.what() << '\n';
                    exit(-1);
                }
            }
            return *this;
        }
        Vector& operator/=(T val)
        {
            try
            {
                if (val == 0)throw std::runtime_error("Division by zero!");
            }
            catch (const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                exit(-1);
            }
            for (int i = 0; i < n; i++)
                v[i] /= v;
            return *this;
        }
        Vector& operator-=(const Vector& A)
        {
            for (int i = 0; i < n; i++)
                v[i] -= A.v[i];
            return *this;
        }
        Vector& operator*=(T val)
        {
            for (int i = 0; i < n; i++)
                v[i] *= val;
            return *this;
        }
        void MakeZero()
        {
            for (int i = 0; i < n; i++)
                v[i] = static_cast<T>(0);
        }
        unsigned int dim() const { return n; }
        ~Vector() { delete[] v; }
};

template<typename T>
T dot(const Vector<T>& A, const Vector<T>& B)
{
    assert(A.dim() == B.dim());
    T ret = static_cast<T>(0);
    for(int i = 0; i < A.dim(); i++)
        ret += A[i] * B[i];
    return ret;
}

template<typename T>
T norm(const Vector<T>& A) { return static_cast<T>(std::sqrt(dot(A, A))); }

using Vec = Vector<double>;
using VecInt = Vector<int>;
using VecFloat = Vector<float>;
#endif //XMATH_VECTOR_H