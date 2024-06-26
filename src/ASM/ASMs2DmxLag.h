// $Id$
//==============================================================================
//!
//! \file ASMs2DmxLag.h
//!
//! \date Dec 28 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 2D %Lagrange mixed FE models.
//!
//==============================================================================

#ifndef _ASM_S2D_MX_LAG_H
#define _ASM_S2D_MX_LAG_H

#include "ASMs2DLag.h"
#include "ASMmxBase.h"
#include <array>


/*!
  \brief Driver for assembly of structured 2D %Lagrange mixed FE models.
  \details This class implements a two-field mixed formulation with Lagrangian
  basis functions. The geometry and the first field are of equal order and
  is one order higher than the second field.
*/

class ASMs2DmxLag : public ASMs2DLag, private ASMmxBase
{
protected:
  //! \brief Implementation of basis function cache.
  class BasisFunctionCache : public ASMs2DLag::BasisFunctionCache
  {
  public:
    //! \brief The constructor initializes the class.
    //! \param pch Patch the cache is for
    BasisFunctionCache(const ASMs2DLag& pch);

    //! \brief Constructor reusing quadrature info from another instance.
    //! \param cache Instance holding quadrature information
    //! \param b Basis to use
    BasisFunctionCache(const ASMs2D::BasisFunctionCache& cache, int b);

    //! \brief Empty destructor.
    virtual ~BasisFunctionCache() = default;

  protected:
    //! \brief Calculates basis function info in a single integration point.
    //! \param[in] gp Integration point on element (0-indexed)
    //! \param[in] red If \e true, returns for reduced integration scheme
    BasisFunctionVals calculatePt(size_t, size_t gp, bool red) const override;

    //! \brief Calculates basis function info in all integration points.
    void calculateAll() override;
  };

public:
  //! \brief The constructor initializes the dimension of each basis.
  ASMs2DmxLag(unsigned char n_s, const CharVec& n_f);
  //! \brief Copy constructor.
  ASMs2DmxLag(const ASMs2DmxLag& patch, const CharVec& n_f = CharVec(2,0));
  //! \brief Empty destructor.
  virtual ~ASMs2DmxLag() {}


  // Methods for model generation
  // ============================

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! the nodal coordinate array, as well as global node and element numbers.
  virtual bool generateFEMTopology();

  //! \brief Clears the contents of the patch, making it empty.
  //! \param[in] retainGeometry If \e true, the spline geometry is not cleared.
  //! This is used to reinitialize the patch after it has been refined.
  virtual void clear(bool retainGeometry);

  //! \brief Returns the number of bases.
  virtual size_t getNoBasis() const { return 2; }
  //! \brief Returns the total number of nodes in this patch.
  virtual size_t getNoNodes(int basis) const;
  //! \brief Returns the number of solution fields.
  virtual unsigned char getNoFields(int basis) const;
  //! \brief Returns the number of DOFs per node.
  //! \param[in] inod 1-based node index local to current patch
  virtual unsigned char getNodalDOFs(size_t inod) const;
  //! \brief Returns the classification of a node.
  //! \param[in] inod 1-based node index local to current patch
  virtual char getNodeType(size_t inod) const;

  //! \brief Initializes the patch level MADOF array for mixed problems.
  virtual void initMADOF(const int* sysMadof);

  //! \brief Connects all matching nodes on two adjacent boundary edges.
  //! \param[in] edge Local edge index of this patch, in range [1,4]
  //! \param neighbor The neighbor patch
  //! \param[in] nedge Local edge index of neighbor patch, in range [1,4]
  //! \param[in] revers Indicates whether the two edges have opposite directions
  //! \param[in] basis The basis to connect (for mixed problems)
  //! \param[in] coordCheck False to disable coordinate checks (periodic connections)
  //! \param[in] thick Thickness of connection
  virtual bool connectPatch(int edge, ASM2D& neighbor, int nedge, bool revers,
                            int basis, bool coordCheck, int thick);

  //! \brief Makes two opposite boundary edges periodic.
  //! \param[in] dir Parameter direction defining the periodic edges
  virtual void closeBoundaries(int dir, int, int);


  // Methods for integration of finite element quantities.
  // These are the main computational methods of the ASM class hierarchy.
  // ====================================================================

  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand,
                         GlobalIntegral& glbInt, const TimeDomain& time);

  //! \brief Evaluates a boundary integral over a patch edge.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index of the boundary edge
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand, int lIndex,
                         GlobalIntegral& glbInt, const TimeDomain& time);


  // Post-processing methods
  // =======================

  //! \brief Extract the primary solution field at the specified nodes.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] nodes 1-based local node numbers to extract solution for
  virtual bool getSolution(Matrix& sField, const Vector& locSol,
                           const IntVec& nodes) const;

  using ASMs2DLag::evalSolution;
  //! \brief Evaluates the secondary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //!
  //! \details The number of visualization points is the same as the order of
  //! the %Lagrange elements by default.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
                            const int*, char = 0) const;

  //! \brief Evaluates the primary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] nf If nonzero, evaluates nf fields on first basis
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const RealArray*, bool, int, int nf) const;

  //! \brief Evaluates the secondary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] gpar Parameter values of the result sampling points
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
                            const RealArray* gpar, bool) const;

  //! \brief Extracts nodal results for this patch from the global vector.
  //! \param[in] globVec Global solution vector in DOF-order
  //! \param[out] nodeVec Nodal result vector for this patch
  //! \param[in] basis Which basis (or 0 for both) to extract nodal values for
  virtual void extractNodeVec(const RealArray& globVec, RealArray& nodeVec,
                              unsigned char, int basis) const;

protected:

  // Internal utility methods
  // ========================

  using ASMs2DLag::getSize;
  //! \brief Returns the number of nodal points in each parameter direction.
  //! \param[out] n1 Number of nodes in first (u) direction
  //! \param[out] n2 Number of nodes in second (v) direction
  //! \param[in] basis Which basis to return size parameters for
  virtual bool getSize(int& n1, int& n2, int basis) const;

private:
  std::vector<size_t> nxx; //!< Basis dimensions in 1st parameter direction
  std::vector<size_t> nyx; //!< Basis dimensions in 2nd parameter direction

  std::vector< std::array<size_t,2> > elem_sizes; //!< Size on each basis
};

#endif
