// $Id$
//==============================================================================
//!
//! \file SplineFields3D.C
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for spline-based finite element vector fields in 3D.
//!
//==============================================================================

#include "GoTools/trivariate/SplineVolume.h"

#include "SplineFields3D.h"
#include "ASMs3D.h"
#include "ItgPoint.h"
#include "CoordinateMapping.h"
#include "Utilities.h"
#include "Vec3.h"


SplineFields3D::SplineFields3D (const ASMs3D* patch,
                                const RealArray& v, char nbasis,
                                int nnf, const char* name)
  : Fields(name), basis(patch->getBasis(nbasis)), vol(patch->getVolume())
{
  const int n1 = basis->numCoefs(0);
  const int n2 = basis->numCoefs(1);
  const int n3 = basis->numCoefs(2);
  nno = n1*n2*n3;

  const int p1 = basis->order(0);
  const int p2 = basis->order(1);
  const int p3 = basis->order(2);
  nelm = (n1-p1+1)*(n2-p2+1)*(n3-p3+1);

  nsd = patch->getNoSpaceDim();

  size_t ofs = 0;
  for (char i = 1; i < nbasis; ++i)
    ofs += patch->getNoNodes(i)*patch->getNoFields(i);
  auto vit = v.begin()+ofs;

  if (nnf == 0)
    nnf = 3;

  nf = nnf;
  int nfc = patch->getNoFields(nbasis);
  values.resize(nno*nf);
  int ndof = nfc*nno;
  auto end = v.size() > ofs+ndof ? vit+ndof : v.end();
  if (nfc == nf)
    std::copy(vit,end,values.begin());
  else
    for (size_t i = 0; i < nno && vit != end; ++i, vit += nfc)
      for (size_t j = 0; j < nf; ++j)
        values[nf*i+j] = *(vit+j);
}


SplineFields3D::SplineFields3D (const Go::SplineVolume* svol,
                                const RealArray& v, int cmp, const char* name)
  : Fields(name), basis(svol), vol(svol)
{
  nsd = 3;
  nf = cmp;
  values = v;
}


bool SplineFields3D::valueNode (size_t node, Vector& vals) const
{
  if (node < 1 || node > nno) return false;

  vals.resize(nf);
  vals.fill(values.ptr()+(node-1)*nf);
  return true;
}


bool SplineFields3D::valueFE (const ItgPoint& x, Vector& vals) const
{
  if (!basis) return false;

  // Evaluate the basis functions at the given point
  Go::BasisPts spline;
#pragma omp critical
  basis->computeBasis(x.u,x.v,x.w,spline);

  // Evaluate the solution field at the given point
  std::vector<int> ip;
  ASMs3D::scatterInd(basis->numCoefs(0),basis->numCoefs(1),basis->numCoefs(2),
		     basis->order(0),basis->order(1),basis->order(2),
		     spline.left_idx,ip);

  Matrix Vnod;
  utl::gather(ip,nf,values,Vnod);
  Vnod.multiply(spline.basisValues,vals); // vals = Vnod * basisValues

  return true;
}


bool SplineFields3D::valueCoor (const Vec4& x, Vector& vals) const
{
  if (x.u)
    return this->valueFE(ItgPoint(x.u[0],x.u[1],x.u[2]),vals);

  // Use with caution, very slow!
  Go::Point pt(x.x,x.y,x.z), clopt(3);
  double clo_u, clo_v, clo_w, dist;
#pragma omp critical
  vol->closestPoint(pt, clo_u, clo_v, clo_w, clopt, dist, 1.0e-5);

  return this->valueFE(ItgPoint(clo_u,clo_v,clo_w),vals);
}


bool SplineFields3D::gradFE (const ItgPoint& x, Matrix& grad) const
{
  if (!basis) return false;
  if (!vol)   return false;

  // Evaluate the basis functions at the given point
  Go::BasisDerivs spline;
#pragma omp critical
  vol->computeBasis(x.u,x.v,x.w,spline);

  const int uorder = vol->order(0);
  const int vorder = vol->order(1);
  const int worder = vol->order(2);
  const size_t nen = uorder*vorder*worder;

  Matrix dNdu(nen,3), dNdX;
  for (size_t n = 1; n <= nen; n++)
  {
    dNdu(n,1) = spline.basisDerivs_u[n-1];
    dNdu(n,2) = spline.basisDerivs_v[n-1];
    dNdu(n,3) = spline.basisDerivs_w[n-1];
  }

  std::vector<int> ip;
  ASMs3D::scatterInd(vol->numCoefs(0),vol->numCoefs(1),vol->numCoefs(2),
		     uorder,vorder,worder,spline.left_idx,ip);

  // Evaluate the Jacobian inverse
  Matrix Xnod(nsd,ip.size()), Jac;
  for (size_t i = 0; i < ip.size(); i++)
    Xnod.fillColumn(1+i,&(*vol->coefs_begin())+vol->dimension()*ip[i]);
  if (!utl::Jacobian(Jac,dNdX,Xnod,dNdu))
    return false; // Singular Jacobian

  // Evaluate the gradient of the solution field at the given point
  if (basis != vol)
  {
    // Mixed formulation, the solution uses a different basis than the geometry
#pragma omp critical
    basis->computeBasis(x.u,x.v,x.w,spline);

    const size_t nbf = basis->order(0)*basis->order(1)*basis->order(2);
    dNdu.resize(nbf,3);
    for (size_t n = 1; n <= nbf; n++)
    {
      dNdu(n,1) = spline.basisDerivs_u[n-1];
      dNdu(n,2) = spline.basisDerivs_v[n-1];
      dNdu(n,3) = spline.basisDerivs_w[n-1];
    }
    dNdX.multiply(dNdu,Jac); // dNdX = dNdu * Jac

    ip.clear();
    ASMs3D::scatterInd(basis->numCoefs(0),basis->numCoefs(1),basis->numCoefs(2),
		       basis->order(0),basis->order(1),basis->order(2),
		       spline.left_idx,ip);
  }

  utl::gather(ip,nf,values,Xnod);
  return !grad.multiply(Xnod,dNdX).empty(); // grad = Xnod * dNdX
}


bool SplineFields3D::hessianFE (const ItgPoint& x, Matrix3D& H) const
{
  if (!basis) return false;
  if (!vol)  return false;

  // Evaluate the basis functions at the given point
  Go::BasisDerivs2 spline2;
  Matrix3D d2Ndu2;
  IntVec ip;
  if (vol == basis) {
#pragma omp critical
    vol->computeBasis(x.u,x.v,x.w,spline2);

    const size_t nen = vol->order(0)*vol->order(1)*vol->order(2);
    d2Ndu2.resize(nen,3,3);
    for (size_t n = 1; n <= nen; n++) {
      d2Ndu2(n,1,1) = spline2.basisDerivs_uu[n-1];
      d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline2.basisDerivs_uv[n-1];
      d2Ndu2(n,1,3) = d2Ndu2(n,3,1) = spline2.basisDerivs_uw[n-1];
      d2Ndu2(n,2,2) = spline2.basisDerivs_vv[n-1];
      d2Ndu2(n,2,3) = d2Ndu2(n,3,2) = spline2.basisDerivs_vw[n-1];
      d2Ndu2(n,3,3) = spline2.basisDerivs_ww[n-1];
    }

    ASMs3D::scatterInd(vol->numCoefs(0),vol->numCoefs(1),vol->numCoefs(2),
		       vol->order(0),vol->order(1),vol->order(2),
		       spline2.left_idx,ip);
  }
  else {
    // Mixed formulation, the solution uses a different basis than the geometry
#pragma omp critical
    basis->computeBasis(x.u,x.v,x.w,spline2);

    const size_t nbf = basis->order(0)*basis->order(1)*basis->order(2);
    d2Ndu2.resize(nbf,3,3);
    for (size_t n = 1; n <= nbf; n++) {
      d2Ndu2(n,1,1) = spline2.basisDerivs_uu[n-1];
      d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline2.basisDerivs_uv[n-1];
      d2Ndu2(n,1,3) = d2Ndu2(n,3,1) = spline2.basisDerivs_uw[n-1];
      d2Ndu2(n,2,2) = spline2.basisDerivs_vv[n-1];
      d2Ndu2(n,2,3) = d2Ndu2(n,3,2) = spline2.basisDerivs_vw[n-1];
      d2Ndu2(n,3,3) = spline2.basisDerivs_ww[n-1];
    }

    ASMs3D::scatterInd(basis->numCoefs(0),basis->numCoefs(1),basis->numCoefs(2),
		       basis->order(0),basis->order(1),basis->order(2),
		       spline2.left_idx,ip);
  }

  Matrix Vnod;
  utl::gather(ip,nf,values,Vnod);
  return H.multiply(Vnod,d2Ndu2);
}
