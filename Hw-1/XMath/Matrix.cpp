#include "Matrix.h"
void FirstElementaryRowOpt(Matrix<Real>& M, uint i, uint j)
{
    for(uint p = 0; p < M.getN(); p++)
        std::swap(M[i][p], M[j][p]);
}
void FirstElementaryColumnOpt(Matrix<Real>& M, uint i, uint j)
{
    for(uint p = 0; p < M.getM(); p++)
        std::swap(M[p][i], M[p][j]);
}
void SecondElementaryRowOpt(Matrix<Real>& M, uint i, Real k, uint j)
{
    try
    {
        if(i == j)throw std::runtime_error("i and j must be different in the second elementary operation");
        for(uint p = 0; p < M.getN(); p++)
            M[j][p] += k * M[i][p];
    }
    catch(std::exception& e)
    {
        std::cerr << e.what() << std::endl;
    }
}
void SecondElementaryColumnOpt(Matrix<Real>& M, uint i, Real k, uint j)
{
    try
    {
        if(i == j)throw std::runtime_error("i and j must be different in the second elementary operation");
        for(uint p = 0; p < M.getM(); p++)
            M[p][j] += k * M[p][i];
    }
    catch(std::exception& e)
    {
        std::cerr << e.what() << std::endl;
    }
}
void ThirdElementaryRowOpt(Matrix<Real>& M, uint i, Real k)
{
    for(uint j = 0; j < M.getN(); j++)
        M[i][j] *= k;
}
void ThirdElementaryColumnOpt(Matrix<Real>& M, uint i, Real k)
{
    for(uint j = 0; j < M.getM(); j++)
        M[j][i] *= k;
}
void FirstElementaryRowOpt(Matrix<Real>& M, uint i, uint j, uint l, uint r)
{
    for(uint p = l; p < r; p++)
        std::swap(M[i][p], M[j][p]);
}
void FirstElementaryColumnOpt(Matrix<Real>& M, uint i, uint j, uint l, uint r)
{
    for(uint p = l; p < r; p++)
        std::swap(M[p][i], M[p][j]);
}
void SecondElementaryRowOpt(Matrix<Real>& M, uint i, Real k, uint j, uint l, uint r)
{
    try
    {
        if(i == j)throw std::runtime_error("i and j must be different in the second elementary operation");
        for(uint p = l; p < r; p++)
            M[j][p] += k * M[i][p];
    }
    catch(std::exception& e)
    {
        std::cerr << e.what() << std::endl;
    }
}
void SecondElementaryColumnOpt(Matrix<Real>& M, uint i, Real k, uint j, uint l, uint r)
{
    try
    {
        if(i == j)throw std::runtime_error("i and j must be different in the second elementary operation");
        for(uint p = l; p < r; p++)
            M[p][j] += k * M[p][i];
    }
    catch(std::exception& e)
    {
        std::cerr << e.what() << std::endl;
    }
}
void ThirdElementaryRowOpt(Matrix<Real>& M, uint i, Real k, uint l, uint r)
{
    for(uint j = l; j < r; j++)
        M[i][j] *= k;
}
void ThirdElementaryColumnOpt(Matrix<Real>& M, uint i, Real k, uint l, uint r)
{
    for(uint j = l; j < r; j++)
        M[j][i] *= k;
}