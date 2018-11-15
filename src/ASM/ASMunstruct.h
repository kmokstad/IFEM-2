// $Id$
//==============================================================================
//!
//! \file ASMunstruct.h
//!
//! \date September 2011
//!
//! \author Kjetil Andre Johannessen / SINTEF
//!
//! \brief Base class for unstructured spline-based FE assembly drivers.
//!
//==============================================================================

#ifndef _ASM_UNSTRUCT_H
#define _ASM_UNSTRUCT_H

#include "ASMbase.h"

typedef std::set<int> IntSet; //!< General integer set


namespace LR //! Utilities for LR-splines.
{
  /*!
    \brief A struct of data to control the mesh refinement.
    \details The \a options parameters have the following interpretation:
    options[0] is the beta percentage of elements to refine,
    options[1] is the knotline multiplicity (default 1),
    options[2] is the refinement scheme (default 0),
    (FULLSPAN=0, MINSPAN=1, ISOTROPIC_ELEMENTS=2, ISOTROPIC_FUNCTIONS=3),
    options[3] is nonzero if testing for linear independence at all iterations,
    options[4] is the maximum number of T-joints allowed in the model,
    options[5] is the maximum allowed parametric aspect ratio of an element,
    options[6] is one if all "gaps" are to be closed,
    options[7] is one if using "true beta".
  */
  struct RefineData
  {
    bool      refShare; //!< If \e true, force refinement of shared FE grids
    IntVec    options;  //!< Parameters used to control the refinement
    IntVec    elements; //!< 0-based indices of the elements to refine
    RealArray errors;   //!< List of error indicators for the elements

    //! \brief Default constructor.
    explicit RefineData(bool rs = false) : refShare(rs) {}
  };
}


/*!
  \brief Base class for unstructured spline-based FE assembly drivers.
  \details This class is an abstract interface for unstructured spline patches.
*/

class ASMunstruct : public ASMbase
{
protected:
  //! \brief The constructor sets the space dimensions.
  //! \param[in] n_p Number of parameter dimensions
  //! \param[in] n_s Number of spatial dimensions
  //! \param[in] n_f Number of primary solution fields
  ASMunstruct(unsigned char n_p, unsigned char n_s, unsigned char n_f)
    : ASMbase(n_p,n_s,n_f) {}

  //! \brief Special copy constructor for sharing of FE data.
  //! \param[in] patch The patch to use FE data from
  //! \param[in] n Number of primary solution fields
  ASMunstruct(const ASMunstruct& patch, unsigned char n) : ASMbase(patch,n) {}

public:
  //! \brief Empty destructor.
  virtual ~ASMunstruct() {}

  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the mesh refinement
  //! \param sol Control point results values that are transferred to new mesh
  //! \param[in] fName Optional file name for an image of the resulting mesh
  virtual bool refine(const LR::RefineData& prm, Vectors& sol,
                      const char* fName = nullptr) = 0;

  //! \brief Remaps element-wise errors from geometry mesh to refinement mesh.
  //! \param[out] errors The remapped errors
  //! \param[in] orig The element-wise errors on the geometry mesh
  virtual void remapErrors(RealArray& errors, const RealArray& orig,
                           bool = false) const { errors = orig; }

  //! \brief Refines the parametrization based on a mesh density function.
  //! \param[in] refC Mesh refinement criteria function
  //! \param[in] refTol Mesh refinement threshold
  virtual bool refine(const RealFunc& refC, double refTol) = 0;

  //! \brief Returns all boundary functions that are covered by the given nodes.
  virtual IntVec getBoundaryCovered(const IntSet&) const { return IntVec(); }
  //! \brief Extends the refinement domain with information for neighbors.
  virtual void extendRefinementDomain(IntSet&, const IntSet&) const {}
};

#endif
