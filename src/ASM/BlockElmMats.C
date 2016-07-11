// $Id$
//==============================================================================
//!
//! \file BlockElmMats.C
//!
//! \date Jul 10 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the block element matrices for a FEM problem.
//!
//==============================================================================

#include "BlockElmMats.h"


bool BlockElmMats::redim (size_t blkIndex, size_t nen, size_t ncmp, char basis)
{
  bool ok = true;
  if (basis < 0) // Zero diagonal block
    basis = -basis;
  else if (blkIndex <= blockInfo.size() && blkIndex < A.size())
    A[blkIndex].resize(nen*ncmp,nen*ncmp);
  else
    ok = false;

  if (blkIndex <= blockInfo.size() && blkIndex < b.size())
    b[blkIndex].resize(nen*ncmp);
  else
    ok = false;

  if (blkIndex == 0)
    return ok;
  else if (blkIndex <= blockInfo.size())
    blockInfo[blkIndex-1] = Block(nen,ncmp,basis);
  else
    ok = false;

  if (basis > 0 && (size_t)basis <= basisInfo.size())
  {
    basisInfo[basis-1].nen = nen;
    basisInfo[basis-1].ncmp += ncmp;
  }
  else
    ok = false;

  return ok;
}


bool BlockElmMats::redimOffDiag (size_t blkIndex, char symm)
{
  size_t nDiagB = blockInfo.size();
  if (blkIndex <= nDiagB || blkIndex > A.size())
    return false; // Not an off-diagonal sub-matrix

  size_t kb = nDiagB+1;
  for (size_t ib = 0; ib < nDiagB; ib++)
    for (size_t jb = 0; jb < nDiagB; jb++, kb++)
      if (jb == ib)
        kb--;
      else if (kb == blkIndex)
      {
        A[kb].resize(blockInfo[ib].nen*blockInfo[ib].ncmp,
                     blockInfo[jb].nen*blockInfo[jb].ncmp);
        kb -= nDiagB;
        if (symmetry.size() < kb)
          symmetry.resize(kb,0);
        symmetry[kb-1] = symm;
        return true;
      }

  return false;
}


void BlockElmMats::redimNewtonMat ()
{
  size_t nDiagB = blockInfo.size();
  if (A.front().empty())
  {
    // Calculate the total number of element dofs
    size_t neldof = 0;
    for (size_t i = 0; i < nDiagB; i++)
      neldof += blockInfo[i].ncmp*blockInfo[i].nen;

    // Set the element Newton matrix dimension
    A.front().resize(neldof,neldof);
    b.front().resize(neldof);
  }

  // Calculate the offset of each block sub-matrix
  size_t idof = 1;
  std::vector<size_t> blockStart(nDiagB,1);
  for (size_t j = 0; j < basisInfo.size(); j++)
  {
    size_t ndof = 0;
    size_t jdof = idof;
    for (size_t i = 0; i < nDiagB; i++)
      if (blockInfo[i].basis == (int)j+1)
      {
        blockInfo[i].idof = jdof;
        jdof += blockInfo[i].ncmp;
        ndof += blockInfo[i].ncmp*blockInfo[i].nen;
      }
    idof += ndof;
  }
}


const Matrix& BlockElmMats::getNewtonMatrix () const
{
  Matrix& N = const_cast<Matrix&>(A.front());

  size_t i, j, k;
  size_t nDiagB = blockInfo.size();
  for (i = 0, k = nDiagB; i < nDiagB; i++)
  {
    size_t ibas = blockInfo[i].basis-1;
    size_t neni = basisInfo[ibas].nen;
    size_t ndi  = basisInfo[ibas].ncmp;
    size_t nci  = blockInfo[i].ncmp;
    size_t idof = blockInfo[i].idof;

    // Insert the diagonal block matrix
    const Matrix& Aii = A[1+i];
    if (!Aii.empty())
      for (size_t in = 0; in < neni; in++)
        for (size_t jn = 0; jn < neni; jn++)
          for (size_t id = 0; id < nci; id++)
            for (size_t jd = 0; jd < nci; jd++)
              N(idof+ndi*in+id,idof+ndi*jn+jd) = Aii(nci*in+id+1,nci*jn+jd+1);

    for (j = 0; j < nDiagB; j++, k++)
    {
      if (j == i) {
        k--;
        continue;
      }

      size_t jbas = blockInfo[j].basis-1;
      size_t nenj = basisInfo[jbas].nen;
      size_t ndj  = basisInfo[jbas].ncmp;
      size_t ncj  = blockInfo[j].ncmp;
      size_t jdof = blockInfo[j].idof;

      // Insert the off-diagonal block matrix
      if (A.size() > 1+k && !A[1+k].empty()) {
        const Matrix& Aij = A[1+k];
        for (size_t in = 0; in < neni; in++)
          for (size_t jn = 0; jn < nenj; jn++)
            for (size_t id = 0; id < nci; id++)
              for (size_t jd = 0; jd < ncj; jd++)
              {
                N(idof+ndi*in+id,jdof+ndj*jn+jd) = Aij(nci*in+id+1,ncj*jn+jd+1);
                if (symmetry[k-nDiagB] > 0)
                  N(jdof+ndj*jn+jd,idof+ndi*in+id) = Aij(nci*in+id+1,ncj*jn+jd+1);
                else if (symmetry[k-nDiagB] < 0)
                  N(jdof+ndj*jn+jd,idof+ndi*in+id) = -Aij(nci*in+id+1,ncj*jn+jd+1);
              }
      }
    }
  }

  return A.front();
}


const Vector& BlockElmMats::getRHSVector () const
{
  Vector& R = const_cast<Vector&>(b.front());

  for (size_t i = 0; i < blockInfo.size(); i++)
  {
    size_t ibas = blockInfo[i].basis-1;
    size_t neni = basisInfo[ibas].nen;
    size_t ndi  = basisInfo[ibas].ncmp;
    size_t nci  = blockInfo[i].ncmp;
    size_t idof = blockInfo[i].idof;

    const Vector& bi = b[1+i];
    for (size_t in = 0; in < neni; in++)
      for (size_t id = 0; id < nci; id++)
        R(idof+ndi*in+id) = bi(nci*in+id+1);
  }

  return b.front();
}
