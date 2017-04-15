// $Id$
//==============================================================================
//!
//! \file SystemMatrix.C
//!
//! \date Jan 4 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General representation of the system matrix on different forms.
//!
//==============================================================================

#include "SystemMatrix.h"
#include "DenseMatrix.h"
#include "SPRMatrix.h"
#include "SparseMatrix.h"
#ifdef HAS_PETSC
#include "PETScMatrix.h"
#endif
#ifdef HAS_ISTL
#include "ISTLMatrix.h"
#endif
#include "LinSolParams.h"


SystemVector* SystemVector::create (const ProcessAdm* adm, Type vectorType)
{
  switch (vectorType)
    {
    case STD   : return new StdVector();
#ifdef HAS_PETSC
    case PETSC : return new PETScVector(*adm);
#endif
#ifdef HAS_ISTL
    case ISTL  : return new ISTLVector(*adm);
#endif
    default:
      std::cerr <<"SystemVector::create: Unsupported vector type "
                << vectorType << std::endl;
    }

  return nullptr;
}


SystemVector& SystemVector::copy (const SystemVector& x)
{
  this->redim(x.size());
  Real* vec = this->getPtr();
  memcpy(vec,x.getRef(),x.dim()*sizeof(Real));
  this->restore(vec);

  return *this;
}


void StdVector::dump (std::ostream& os, char format, const char* label)
{
  switch (format)
    {
    case 'M':
    case 'm':
      utl::writeMatlab(label,*this,os);
      break;

    default:
      if (label) os << label <<" =";
      os << *this;
    }
}


SystemMatrix* SystemMatrix::create (const ProcessAdm* adm, Type matrixType,
                                    const LinSolParams& spar,
                                    LinAlg::LinearSystemType ltype)
{
#ifdef HAS_PETSC
  if (matrixType == PETSC)
    return new PETScMatrix(*adm,spar,ltype);
#endif
#ifdef HAS_ISTL
  if (matrixType == ISTL)
    return new ISTLMatrix(*adm,spar,ltype);
#endif

  return SystemMatrix::create(adm,matrixType,ltype);
}


SystemMatrix* SystemMatrix::create (const ProcessAdm* adm, Type matrixType,
                                    LinAlg::LinearSystemType ltype,
                                    int num_thread_SLU)
{
#ifndef HAS_PETSC
  if (matrixType == PETSC) {
    std::cerr <<"SystemMatrix::create: PETSc not compiled in, bailing out..."
              << std::endl;
    exit(1);
  }
#endif
#ifndef HAS_ISTL
  if (matrixType == ISTL) {
    std::cerr <<"SystemMatrix::create: ISTL not compiled in, bailing out..."
              << std::endl;
    exit(1);
  }
#endif

#if defined(HAS_PETSC) || defined(HAS_ISTL)
  // Use default settings when no parameters are provided by user
  static LinSolParams defaultPar;
#endif

  switch (matrixType)
    {
    case DENSE : return new DenseMatrix();
    case SPR   : return new SPRMatrix();
    case SPARSE: return new SparseMatrix(SparseMatrix::SUPERLU,num_thread_SLU);
    case SAMG  : return new SparseMatrix(SparseMatrix::S_A_M_G);
#ifdef HAS_PETSC
    case PETSC : return new PETScMatrix(*adm,defaultPar,ltype);
#endif
#ifdef HAS_ISTL
    case ISTL  : return new ISTLMatrix(*adm,defaultPar,ltype);
#endif
    default:
      std::cerr <<"SystemMatrix::create: Unsupported matrix type "
                << matrixType << std::endl;
    }

  return 0;
}


bool SystemMatrix::assemble (const Matrix&, const SAM&,
                             SystemVector&, const std::vector<int>&)
{
  std::cerr <<"SystemMatrix::assemble(const Matrix&,const SAM&,"
            <<"SystemVector&,const std::vector<int>&): "
            <<"Not implemented for the chosen matrix type."<< std::endl;
  return false;
}


StdVector SystemMatrix::operator* (const StdVector& b) const
{
  StdVector results;
  this->multiply(b,results);
  return results;
}


StdVector SystemMatrix::operator/ (const StdVector& b)
{
  StdVector results;
  this->solve(b,results);
  return results;
}
