// $Id$
//==============================================================================
//!
//! \file Chebyshev.C
//!
//! \date Jul 13 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Evaluation of Chebyshev polynomials.
//!
//==============================================================================

#include "Chebyshev.h"
#include "MatVec.h"

#include <fstream>
#include <numeric>


double Chebyshev::evalPol1 (int polnum, double xi)
{
  if (polnum <= 0)
    return 1.0;
  if (polnum == 1)
    return xi;

  return 2.0*xi*evalPol1(polnum-1, xi) - evalPol1(polnum-2, xi);
}


double Chebyshev::evalPol2 (int polnum, double xi)
{
  if (polnum <= 0)
    return 1.0;
  if (polnum == 1)
    return 2.0*xi;

  return 2.0*xi*evalPol2(polnum-1, xi) - evalPol2(polnum-2, xi);
}


double Chebyshev::evalDer1 (int polnum, double xi)
{
  if (polnum <= 0)
    return 0.0;

  return polnum * evalPol2(polnum-1, xi);
}


double Chebyshev::evalDer2 (int polnum, double xi)
{
  return ((polnum + 1)*evalPol1(polnum + 1, xi) - xi*evalPol2(polnum, xi)) / (xi*xi - 1.0);
}


double Chebyshev::eval2Der1 (int polnum, double xi)
{
  if (polnum < 2)
    return 0.0;

  if (std::abs(xi-1.0) < 1e-6)
    return (pow(polnum,4) - pow(polnum, 2)) / 3.0;
  else if (std::abs(1.0+xi) < 1e-6)
    return pow(-1, polnum) * (pow(polnum,4) - pow(polnum, 2)) / 3.0;

  return polnum * ((polnum+1)*evalPol1(polnum, xi) - evalPol2(polnum, xi)) / (xi*xi - 1.0);
}


ChebyshevFunc::ChebyshevFunc(const char* file)
{
  std::ifstream in(file);
  if (!in.good()) {
    n[0] = n[1] = n[2] = 0;
    return;
  }

  in >> n[0] >> n[1] >> n[2];
  if (n[1] == 0)
    n[1] = 1;
  if (n[2] == 0)
    n[2] = 1;
  coefs.resize(n[0]*n[1]*n[2]);
  for (double& coef : coefs)
    in >> coef;
}


Real ChebyshevFunc::evaluate(const Vec3& X) const
{
  const Vec4& X4 = static_cast<const Vec4&>(X);
  double res = 0.0;
  Vector TX(n[0]);
  Matrix T;
  for (int i = 0; i < n[0]; ++i)
    TX[i] = Chebyshev::evalPol1(i, (-1.0+2.0*X4.u[0]));
  Vector TY(n[1]);
  for (int j = 0; j < n[1]; ++j)
    TY[j] = Chebyshev::evalPol1(j, (-1.0+2.0*X4.u[1]));
  T.outer_product(TX, TY);
  if (n[2] == 1)
    res = std::inner_product(coefs.begin(), coefs.end(), T.begin(), 0.0);
  else {
    Vector TZ(n[2]);
    for (int k = 0; k < n[2]; ++k)
      TZ[k] = Chebyshev::evalPol1(k, (-1.0 + 2.0*X4.u[2]));
    Matrix T2;
    T2.outer_product(T, TZ);
    res = std::inner_product(coefs.begin(), coefs.end(), T2.begin(), 0.0);
  }

  return res;
}


ChebyshevVecFunc::ChebyshevVecFunc(const std::vector<const char*>& file, bool second)
  : secondDer(second)
{
  f[0].reset(new ChebyshevFunc(file[0]));
  if (file.size() > 1)
    f[1].reset(new ChebyshevFunc(file[1]));
  if (file.size() > 2)
    f[2].reset(new ChebyshevFunc(file[2]));
  if (f[0]->getSize()[2] == 1)
    ncmp = 2;
}


Vec3 ChebyshevVecFunc::evaluate(const Vec3& X) const
{
  const Vec4& X4 = static_cast<const Vec4&>(X);
  Vec3 res;

  // Multiple-components - no derivatives
  if (f[1]) {
    res[0] = (*f[0])(X);
    res[1] = (*f[1])(X);
    if (f[2])
      res[2] = (*f[2])(X);
    return res;
  }

  const std::array<int,3>& n = f[0]->getSize();
  const std::vector<Real>& coefs = f[0]->getCoefs();

  Vector TX(n[0]), dTX(n[0]);
  Matrix T;
  for (int i = 0; i < n[0]; ++i) {
    TX[i] = Chebyshev::evalPol1(i, (-1.0+2.0*X4.u[0]));
    if (secondDer)
      dTX[i] = 4.0*Chebyshev::eval2Der1(i, (-1.0+2.0*X4.u[0])); // 4.0 due to dxi/du twice
    else
      dTX[i] = 2.0*Chebyshev::evalDer1(i, (-1.0+2.0*X4.u[0])); // 2.0 due to dxi/du
  }
  Vector TY(n[1]), dTY(n[1]);
  for (int j = 0; j < n[1]; ++j) {
    TY[j] = Chebyshev::evalPol1(j, (-1.0+2.0*X4.u[1]));
    if (secondDer)
      dTY[j] = 4.0*Chebyshev::eval2Der1(j, (-1.0+2.0*X4.u[1])); // 4.0 due to dxi/du twice
    else
      dTY[j] = 2.0*Chebyshev::evalDer1(j, (-1.0+2.0*X4.u[1])); // 2.0 due to dxi/du
  }
  if (n[2] == 1) {
    T.outer_product(dTX, TY);
    res[0] = std::inner_product(coefs.begin(), coefs.end(), T.begin(), 0.0);
    T.outer_product(TX, dTY);
    res[1] = std::inner_product(coefs.begin(), coefs.end(), T.begin(), 0.0);
  } else {
    Vector TZ(n[2]), dTZ(n[2]);
    for (int k = 0; k < n[2]; ++k) {
      TZ[k] = Chebyshev::evalPol1(k, (-1.0 + 2.0*X4.u[2]));
      if (secondDer)
        dTZ[k] = Chebyshev::eval2Der1(k, (-1.0 + 2.0*X4.u[2])); // 4.0 due to dxi/du twice
      else
        dTZ[k] = Chebyshev::evalDer1(k, (-1.0 + 2.0*X4.u[2])); // 2.0 due to dxi/du
    }
    Matrix T2;
    T.outer_product(dTX, TY);
    T2.outer_product(T, TZ);
    res[0] = std::inner_product(coefs.begin(), coefs.end(), T2.begin(), 0.0);
    T.outer_product(TX, dTY);
    T2.outer_product(T, TZ);
    res[1] = std::inner_product(coefs.begin(), coefs.end(), T2.begin(), 0.0);
    T.outer_product(TX, TY);
    T2.outer_product(T, dTZ);
    res[2] = std::inner_product(coefs.begin(), coefs.end(), T2.begin(), 0.0);
  }

  return res;
}


ChebyshevTensorFunc::ChebyshevTensorFunc(const std::vector<const char*>& file, bool second)
{
  if (file.size() < 4) {
    for (size_t i = 0; i < file.size(); ++i)
      f[i].reset(new ChebyshevVecFunc({file[i]}, second));
  } else {
    if (file.size() == 4) {
      f[0].reset(new ChebyshevVecFunc({file[0], file[1]}));
      f[1].reset(new ChebyshevVecFunc({file[2], file[3]}));
    } else {
      f[0].reset(new ChebyshevVecFunc({file[0], file[1], file[2]}));
      f[1].reset(new ChebyshevVecFunc({file[3], file[4], file[5]}));
      f[2].reset(new ChebyshevVecFunc({file[6], file[7], file[8]}));
    }
  }
  ncmp = f[0]->dim();
}


Tensor ChebyshevTensorFunc::evaluate(const Vec3& X) const
{
  Tensor result(ncmp);
  Vec3 r1 = (*f[0])(X);
  Vec3 r2, r3;
  if (f[1])
    r2 = (*f[1])(X);
  if (f[2])
    r3 = (*f[2])(X);
  for (size_t i = 0; i < ncmp; ++i)
    result(1, i+1) = r1[i];
  if (f[1])
    for (size_t i = 0; i < ncmp; ++i)
      result(2, 1+i) = r2[i];
  if (f[2])
    for (size_t i = 0; i < ncmp; ++i)
      result(3, 1+i) = r3[i];

  return result;
}
