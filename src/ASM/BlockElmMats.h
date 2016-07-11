// $Id$
//==============================================================================
//!
//! \file BlockElmMats.h
//!
//! \date Jul 10 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the block element matrices for a FEM problem.
//!
//==============================================================================

#ifndef _BLOCK_ELM_MATS_H
#define _BLOCK_ELM_MATS_H

#include "ElmMats.h"


/*!
  \brief Class representing the block element matrices for a FEM problem.
  \details This class can be used for multi-integrand FEM problems with strong
  coupling, where the element matrices typically consists of blocks resulting
  from the various integrands. The block sub-matrices are stiched together in a
  element-level Newton matrix corresponding to all DOFs of the element by the
  getNewtonMatrix method.
*/

class BlockElmMats : public ElmMats
{
public:
  //! \brief The constructor initializes the block information arrays.
  //! \param[in] nBlk Number of matrix blocks (in each direction, row & column)
  //! \param[in] nb Number of bases (> 1 for mixed problems)
  BlockElmMats(size_t nBlk, size_t nb = 1) : blockInfo(nBlk), basisInfo(nb) {}
  //! \brief Empty destructor.
  virtual ~BlockElmMats() {}

  using ElmMats::redim;
  //! \brief Sets the dimension of a diagonal block sub-matrix.
  //! \param[in] blkIndex Sub-matrix block index
  //! \param[in] nen Number of element nodes for this matrix block
  //! \param[in] ncmp Number of nodal components for this matrix block
  //! \param[in] basis Basis index associated with this block. If negative,
  //! this diagonal block is identically zero, but the dimension parameters are
  //! used for the associated off-diagonal blocks later.
  bool redim(size_t blkIndex, size_t nen, size_t ncmp = 1, char basis = 1);
  //! \brief Sets the dimension of an off-diagonal block sub-matrix.
  //! \param[in] blkIndex Sub-matrix block index
  //! \param[in] symm Symmetry flag (0: not, 1: yes, -1: skew-symmetric)
  //!
  //! \details This method must not be invoked until the \a redim method has
  //! been invoked for all the diagonal blocks.
  bool redimOffDiag(size_t blkIndex, char symm = 1);
  //! \brief Sets the dimension of the Newton matrix.
  //! \details This method must not be invoked until the \a redim method has
  //! been invoked for all the diagonal blocks.
  void redimNewtonMat();

  //! \brief Returns the element-level Newton matrix.
  virtual const Matrix& getNewtonMatrix() const;

  //! \brief Returns the element-level right-hand-side vector
  //! associated with the Newton matrix.
  virtual const Vector& getRHSVector() const;

private:
  //! \brief A struct with some key parameters for each block.
  struct Block
  {
    char   basis; //!< Basis associated with the block
    size_t nen;   //!< Number of element nodes for the basis of the block
    size_t ncmp;  //!< Number of nodal components for the block
    size_t idof;  //!< The first element DOF of the block
    //! \brief Default constructor.
    Block(size_t n = 0, size_t c = 0, char b = 0) : basis(b), nen(n), ncmp(c) {}
  };

  //! \brief A struct with some key parameters for each basis.
  struct Basis
  {
    size_t nen;  //!< Number of element nodes for the basis
    size_t ncmp; //!< Number of nodal components for the basis
    size_t idof; //!< First element dof for the basis
    //! \brief Default constructor.
    Basis() : nen(0), ncmp(0), idof(1) {}
  };

  std::vector<Block> blockInfo; //!< Block information
  std::vector<Basis> basisInfo; //!< Basis information
  std::vector<char>  symmetry;  //!< Symmetry flag for off-diagonal blocks
};

#endif
