// $Id$
//==============================================================================
//!
//! \file ISTLMatrix.C
//!
//! \date Mar 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Representation of the system matrix in ISTL format.
//!
//==============================================================================

#include "ISTLMatrix.h"
#include "LinSolParams.h"
#include "ProcessAdm.h"
#include "SAM.h"
#include "LinAlgInit.h"


ISTLVector::ISTLVector(const ProcessAdm& padm) : adm(padm)
{
  LinAlgInit::increfs();
}


ISTLVector::ISTLVector(const ProcessAdm& padm, size_t n)
  : StdVector(n), adm(padm)
{
  x.resize(n);
  LinAlgInit::increfs();
}


ISTLVector::ISTLVector(const ProcessAdm& padm, const Real* values, size_t n)
  : StdVector(values,n), adm(padm)
{
  x.resize(n);
  this->endAssembly();
  LinAlgInit::increfs();
}


ISTLVector::ISTLVector(const ISTLVector& vec) : StdVector(vec), adm(vec.adm)
{
  x = vec.x;
  LinAlgInit::increfs();
}


ISTLVector::~ISTLVector()
{
  LinAlgInit::decrefs();
}


void ISTLVector::init(Real value)
{
  x = value;
  this->StdVector::init(value);
}


void ISTLVector::redim(size_t n)
{
  x.resize(n);
  this->StdVector::redim(n);
}


bool ISTLVector::endAssembly()
{
  for (size_t i = 0; i < x.size(); ++i)
    x[i] = (*this)[i];

  return true;
}


Real ISTLVector::L1norm() const
{
  return x.one_norm();
}


Real ISTLVector::L2norm() const
{
  return x.two_norm();
}


Real ISTLVector::Linfnorm() const
{
  return x.infinity_norm();
}


ISTLMatrix::ISTLMatrix (const ProcessAdm& padm, const LinSolParams& spar)
  : SparseMatrix(SUPERLU,1), adm(padm), solParams(spar,adm)
{
  LinAlgInit::increfs();
}


ISTLMatrix::ISTLMatrix (const ISTLMatrix& B)
  : SparseMatrix(B), adm(B.adm), solParams(B.solParams.get(),B.adm)
{
  iA = B.iA;

  LinAlgInit::increfs();
}


ISTLMatrix::~ISTLMatrix ()
{
  LinAlgInit::decrefs();
}


void ISTLMatrix::initAssembly (const SAM& sam, char)
{
  this->resize(sam.getNoEquations());
  this->preAssemble(sam,false);

  std::vector<IntSet> dofc;
  sam.getDofCouplings(dofc);

  // Set correct number of rows and columns for matrix.
  size_t sum = 0;
  for (const IntSet& dofs : dofc)
    sum += dofs.size();

  iA.setSize(rows(), cols(), sum);
  iA.setBuildMode(ISTL::Mat::random);

  for (size_t i = 0; i < dofc.size(); ++i)
    iA.setrowsize(i,dofc[i].size());
  iA.endrowsizes();

  for (size_t i = 0; i < dofc.size(); ++i)
    for (int dof : dofc[i])
      iA.addindex(i,dof-1);

  iA.endindices();

  iA = 0;
}


bool ISTLMatrix::endAssembly()
{
  for (size_t j = 0; j < cols(); ++j)
    for (int i = IA[j]; i < IA[j+1]; ++i)
      iA[JA[i]][j] = A[i];

  return this->SparseMatrix::endAssembly();
}


void ISTLMatrix::init ()
{
  this->SparseMatrix::init();

  // Set all matrix elements to zero
  iA = 0;
}


bool ISTLMatrix::solve (SystemVector& B, Real*)
{
  this->handleSolverReset();

  if (!pre)
    std::tie(solver, pre, op) = solParams.setupPC(iA);

  ISTLVector* Bptr = dynamic_cast<ISTLVector*>(&B);
  if (!Bptr || !solver || !pre)
    return false;

  try {
    Dune::InverseOperatorResult r;
    ISTL::Vec b(Bptr->getVector());
    Bptr->getVector() = 0;
    solver->apply(Bptr->getVector(), b, r);
  } catch (Dune::ISTLError& e) {
    std::cerr << "ISTL exception " << e << std::endl;
    return false;
  }

  for (size_t i = 0; i < rows(); ++i)
    (*Bptr)(i+1) = Bptr->getVector()[i];

  ++nLinSolves;

  return true;
}


bool ISTLMatrix::solve (const SystemVector& b, SystemVector& x)
{
  this->handleSolverReset();

  if (!pre)
    std::tie(solver, pre, op) = solParams.setupPC(iA);

  const ISTLVector* Bptr = dynamic_cast<const ISTLVector*>(&b);
  if (!Bptr || ! solver || !pre)
    return false;

  ISTLVector* Xptr = dynamic_cast<ISTLVector*>(&x);
  if (!Xptr)
    return false;

  try {
    Dune::InverseOperatorResult r;
    solver->apply(Xptr->getVector(),
                  const_cast<ISTL::Vec&>(Bptr->getVector()), r);
  } catch (Dune::ISTLError& e) {
    std::cerr << "ISTL exception " << e << std::endl;
    return false;
  }

  for (size_t i = 0; i < rows(); ++i)
    (*Xptr)(i+1) = Xptr->getVector()[i];

  ++nLinSolves;

  return true;
}


Real ISTLMatrix::Linfnorm () const
{
  return iA.infinity_norm();
}


void ISTLMatrix::handleSolverReset ()
{
  // Reset linear solver
  if (nLinSolves && solParams.get().hasValue("reset_pc")) {
    const std::string string_val = solParams.get().getStringValue("reset_pc");
    int val = solParams.get().getIntValue("reset_pc");
    if (string_val == "all" ||
        (string_val == "first" && nLinSolves == 1) ||
        (val > 0 && nLinSolves % val == 0)) {
      pre.reset();
      solver.reset();
      op.reset();
      adm.cout << "Resetting preconditioner" << std::endl;
    }
  }
}
