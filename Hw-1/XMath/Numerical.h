#ifndef XMATH_NUMERICAL_H
#define XMATH_NUMERICAL_H
//#include "Polynomial.h"
using uint = unsigned int;
using Real = double;
const Real EPS = 1e-12;
#include "Matrix.h"
//void DFT(Vec&);
//void IDFT(Vec&);
void LUDecompose(const Matrix<Real>&, Matrix<Real>&, Matrix<Real>&, Matrix<Real>&);
void inverse(const Matrix<Real>&, Matrix<Real>&);
bool fequals(Real, Real);
int sign(Real);
#endif //XMATH_NUMERICAL_H