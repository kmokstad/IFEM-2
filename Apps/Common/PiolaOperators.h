//==============================================================================
//!
//! \file PiolaOperators.h
//!
//! \date Apr 30 2024
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various weak, discrete Piola-mapped operators.
//!
//==============================================================================

#ifndef PIOLA_OPERATORS_H_
#define PIOLA_OPERATORS_H_

class FiniteElement;
class Tensor;
class Vec3;

#include "MatVec.h"
#include "EqualOrderOperators.h"

/*! \brief Common operators using Piola mapped discretizations.
 *  \details The operators use the block ordering used in the BlockElmMats class.
 */

class PiolaOperators
{
public:
  //! \brief Common weak operators using Piola mapped discretizations.
  class Weak
  {
  public:
    //! \brief Compute an advection term.
    //! \param[out] EM The element matrices to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] AC Advecting field
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] cnvForm Form of advection operator
    static void Advection(Matrices& EM, const FiniteElement& fe,
                          const Vec3& AC, double scale = 1.0,
                          WeakOperators::ConvectionForm cnvForm = WeakOperators::CONVECTIVE);

    //! \brief Compute a (nonlinear) convection term.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] U  Advecting field
    //! \param[in] dUdX Field gradient
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] form Which form of the convective form to use
    static void Convection(Matrices& EM, const FiniteElement& fe,
                           const Vec3& U, const Tensor& dUdX, double scale,
                           WeakOperators::ConvectionForm form = WeakOperators::CONVECTIVE);

    //! \brief Compute a gradient term.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    static void Gradient(Matrices& EM,
                         const FiniteElement& fe, double scale = 1.0);

    //! \brief Compute a laplacian.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] stress Whether to add extra stress formulation terms
    //! \details Only the upper blocks are added with the stress formulation
    static void Laplacian(Matrices& EM,
                          const FiniteElement& fe,
                          double scale, bool stress);

    //! \brief Compute a mass term.
    //! \param[out] EM The element matrices to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    static void Mass(Matrices& EM,
                     const FiniteElement& fe, double scale);

      //! \brief Compute a vector-source term.
      //! \param[out] EV The element vectors to add contribution to
      //! \param[in] fe The finite element to evaluate for
      //! \param[in] f Vector with contributions
      //! \param[in] scale Scaling factor for contribution
    static void Source(Vectors& EV, const FiniteElement& fe,
                       const Vec3& f, double scale);
  };

  //! \brief Common weak residual operators using div-compatible discretizations.
  class Residual
  {
  public:
    //! \brief Compute a convection term in a residual vector.
    //! \param EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] U Advected field
    //! \param[in] dUdX Advected field gradient
    //! \param[in] UC Advecting field
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] form Which form of the convective form to use
    static void Convection(Vectors& EV, const FiniteElement& fe,
                           const Vec3& U, const Tensor& dUdX,
                           const Vec3& UC, double scale,
                           WeakOperators::ConvectionForm form = WeakOperators::CONVECTIVE);

    //! \brief Compute a gradient term.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    static void Gradient(Vectors& EV, const FiniteElement& fe,
                         double scale = 1.0);

    //! \brief Compute a laplacian term in a residual vector.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] dUdX Current solution gradient
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] stress Whether to add extra stress formulation terms
    static void Laplacian(Vectors& EV, const FiniteElement& fe,
                          const Tensor& dUdX, double scale = 1.0,
                          bool stress = false);
  };

  //! \brief Add the full Piola operator to the blocks.
  //! \param EM Vector of block matrices
  //! \param fe Finite element information
  //! \param A Piola operator to add
  static void Copy(Matrices& EM,
                   const FiniteElement& fe,
                   const Matrix& A);

  //! \brief Add the full Piola vector to the blocks.
  //! \param EV Vector of block vectors
  //! \param fe Finite element information
  //! \param V Piola vector to add
  static void Copy(Vectors& EV,
                   const FiniteElement& fe,
                   const Vector& V);
};

#endif
