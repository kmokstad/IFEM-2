// $Id$
//==============================================================================
//!
//! \file Functions.h
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Specific function implementations.
//!
//==============================================================================

#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H

#include "Function.h"

class TensorFunc;
namespace tinyxml2 { class XMLElement; }

using VecTimeFunc = utl::Function<Real,Vec3>; //!< Vector-valued real function
using IntFunc     = utl::Function<int,Real>;  //!< Real-valued integer function


/*!
  \brief A scalar-valued constant function.
*/

class ConstantFunc : public ScalarFunc
{
  Real fval; //!< The function value

public:
  //! \brief Constructor initializing the function value.
  explicit ConstantFunc(Real v) : fval(v) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return fval == Real(0); }

protected:
  //! \brief Evaluates the constant function.
  Real evaluate(const Real&) const override { return fval; }
};


/*!
  \brief A scalar-valued linear function.
*/

class LinearFunc : public ScalarFunc
{
  using Point = std::pair<Real,Real>; //!< Convenience type

  std::vector<Point> fvals; //!< Values for piece-wise linear function
  Real               scale; //!< Scaling factor

public:
  //! \brief Constructor initializing the scaling parameter.
  explicit LinearFunc(Real s = Real(1)) : scale(s) {}

  //! \brief Constructor initializing piece-wise linear function values.
  //! \param[in] file Name of file to read the function values from
  //! \param[in] c Which column of the file to read function values from
  //! \param[in] s Optional scaling factor
  explicit LinearFunc(const char* file, int c = 2, Real s = Real(1));

  //! \brief Constructor initializing piece-wise linear function values.
  //! \param[in] x List of function argument values
  //! \param[in] y List of function values associated with the \a x values
  LinearFunc(const std::vector<Real>& x, const std::vector<Real>& y);

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override;
  //! \brief Returns whether the function is constant or not.
  bool isConstant() const override { return fvals.size() == 1; }

  //! \brief Returns the first-derivative of the function.
  Real deriv(Real x) const override;

protected:
  //! \brief Evaluates the function at \a x.
  Real evaluate(const Real& x) const override;

private:
  //! \brief Retuns the index of the first point after \a x.
  size_t locate(Real x) const;
};


/*!
  \brief A vector-valued linear function.
*/

class LinVecFunc : public VecTimeFunc
{
  using Point = std::pair<Real,Vec3>; //!< Convenience type

  std::vector<Point> fvals; //!< Values for piece-wise linear function

public:
  //! \brief Constructor initializing piece-wise linear function values.
  //! \param[in] file Name of file to read the function values from
  //! \param[in] c The column of the file to start reading function values from
  explicit LinVecFunc(const char* file, int c = 2);

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override;
  //! \brief Returns whether the function is constant or not.
  bool isConstant() const override { return fvals.size() == 1; }

protected:
  //! \brief Evaluates the function at \a x.
  Vec3 evaluate(const Real& x) const override;
};


/*!
  \brief A scalar-valued ramp function, linear up to \a xmax.
*/

class RampFunc : public ScalarFunc
{
  Real fval; //!< Max function value
  Real xmax; //!< The function is linear from \a x = 0 to \a x = \a xmax

public:
  //! \brief Constructor initializing the function parameters.
  explicit RampFunc(Real f = Real(1), Real x = Real(1)) : fval(f), xmax(x) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return fval == Real(0); }
  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return xmax <= Real(0); }

  //! \brief Returns the first-derivative of the function.
  Real deriv(Real x) const override;

protected:
  //! \brief Evaluates the function at \a x.
  Real evaluate(const Real& x) const override;
};


/*!
  \brief A scalar-valued dirac function.
*/

class DiracFunc : public ScalarFunc
{
  Real amp;  //!< The amplitude of the dirac function
  Real xmax; //!< Associated \a x value

public:
  //! \brief Constructor initializing the function parameters.
  explicit DiracFunc(Real a = Real(1), Real x = Real(0)) : amp(a), xmax(x) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return amp == Real(0); }

protected:
  //! \brief Evaluates the function at \a x.
  Real evaluate(const Real& x) const override;
};


/*!
  \brief A scalar-valued step function.
*/

class StepFunc : public ScalarFunc
{
  Real amp;  //!< The amplitude of the step function
  Real xmax; //!< Associated \a x value

public:
  //! \brief Constructor initializing the function parameters.
  explicit StepFunc(Real a, Real x = Real(0)) : amp(a), xmax(x) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return amp == Real(0); }

protected:
  //! \brief Evaluates the function at \a x.
  Real evaluate(const Real& x) const override;
};


/*!
  \brief A scalar-valued sinusoidal function.
*/

class SineFunc : public ScalarFunc
{
  Real scale; //!< Amplitude of the sine function
  Real freq;  //!< Angular frequency of the sine function
  Real phase; //!< Phase shift of the sine function

public:
  //! \brief Constructor initializing the function parameters.
  explicit SineFunc(Real s = Real(1), Real f = Real(1), Real p = Real(0))
    : scale(s), freq(f), phase(p) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return scale == Real(0); }
  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return this->isZero(); }

  //! \brief Returns the first-derivative of the function.
  Real deriv(Real x) const override;

protected:
  //! \brief Evaluates the function at \a x.
  Real evaluate(const Real& x) const override;
};


/*!
  \brief A scalar-valued spatial function, constant in space and time.
*/

class ConstFunc : public RealFunc
{
  Real fval; //!< The function value

public:
  //! \brief Constructor initializing the function value.
  explicit ConstFunc(Real v) : fval(v) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return fval == Real(0); }

protected:
  //! \brief Evaluates the constant function.
  Real evaluate(const Vec3&) const override { return fval; }
};


/*!
  \brief A scalar-valued spatial function, constant in space, varying in time.
*/

class ConstTimeFunc : public RealFunc
{
  const ScalarFunc* tfunc; //!< The time-dependent function value

public:
  //! \brief Constructor initializing the function value.
  explicit ConstTimeFunc(const ScalarFunc* f) : tfunc(f) {}
  //! \brief The destructor frees the time function.
  virtual ~ConstTimeFunc() { delete tfunc; }

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return tfunc->isZero(); }
  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return tfunc->isConstant(); }

  //! \brief Returns first-derivative of the function.
  Real deriv(const Vec3& X, int dir) const override;

protected:
  //! \brief Evaluates the time-varying function.
  Real evaluate(const Vec3& X) const override;
};


/*!
  \brief A scalar-valued spatial function, varying in space and time.
  \details The function value is defined as a product between one
  space-dependent component and one time-dependent component.
*/

class SpaceTimeFunc : public RealFunc
{
  const RealFunc*   sfunc; //!< The space-dependent term
  const ScalarFunc* tfunc; //!< The time-dependent term

public:
  //! \brief Constructor initializing the function terms.
  SpaceTimeFunc(const RealFunc* s, const ScalarFunc* t) : sfunc(s), tfunc(t) {}
  //! \brief The destructor frees the space and time functions.
  virtual ~SpaceTimeFunc() { delete sfunc; delete tfunc; }

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return sfunc->isZero() || tfunc->isZero(); }
  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return tfunc->isConstant(); }

  //! \brief Returns first-derivative of the function.
  Real deriv(const Vec3& X, int dir) const override;
  //! \brief Returns second-derivative of the function.
  Real dderiv(const Vec3& X, int dir1, int dir2) const override;

protected:
  //! \brief Evaluates the space-time function.
  Real evaluate(const Vec3& X) const override;
};


/*!
  \brief A scalar-valued spatial function, linear in \a x.
*/

class LinearXFunc : public RealFunc
{
  Real a; //!< The function derivative
  Real b; //!< The function value at \a x = 0

public:
  //! \brief Constructor initializing the function parameters.
  explicit LinearXFunc(Real A, Real B = Real(0)) : a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return a == Real(0) && b == Real(0); }

  //! \brief Returns first-derivative of the function.
  Real deriv(const Vec3&, int dir) const override;

protected:
  //! \brief Evaluates the linear function.
  Real evaluate(const Vec3& X) const override;
};


/*!
  \brief A scalar-valued spatial function, linear in \a y.
*/

class LinearYFunc : public RealFunc
{
  Real a; //!< The function derivative
  Real b; //!< The function value at \a y = 0

public:
  //! \brief Constructor initializing the function parameters.
  explicit LinearYFunc(Real A, Real B = Real(0)) : a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return a == Real(0) && b == Real(0); }

  //! \brief Returns first-derivative of the function.
  Real deriv(const Vec3&, int dir) const override;

protected:
  //! \brief Evaluates the linear function.
  Real evaluate(const Vec3& X) const override;
};


/*!
  \brief A scalar-valued spatial function, linear in \a z.
*/

class LinearZFunc : public RealFunc
{
  Real a; //!< The function derivative
  Real b; //!< The function value at \a z = 0

public:
  //! \brief Constructor initializing the function parameters.
  explicit LinearZFunc(Real A, Real B = Real(0)) : a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return a == Real(0) && b == Real(0); }

  //! \brief Returns first-derivative of the function.
  Real deriv(const Vec3&, int dir) const override;

protected:
  //! \brief Evaluates the linear function.
  Real evaluate(const Vec3& X) const override;
};


/*!
  \brief A scalar-valued spatial function, quadratic in \a x.
*/

class QuadraticXFunc : public RealFunc
{
  Real max; //!< Max value of function
  Real a;   //!< First root where function is zero
  Real b;   //!< Second root where function is zero

public:
  //! \brief Constructor initializing the function parameters.
  QuadraticXFunc(Real MAX, Real A, Real B) : max(MAX), a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return max == Real(0); }

  //! \brief Returns first-derivative of the function.
  Real deriv(const Vec3& X, int dir) const override;
  //! \brief Returns second-derivative of the function.
  Real dderiv(const Vec3&, int dir1, int dir2) const override;

protected:
  //! \brief Evaluates the quadratic function.
  Real evaluate(const Vec3& X) const override;
};


/*!
  \brief A scalar-valued spatial function, quadratic in \a y.
*/

class QuadraticYFunc : public RealFunc
{
  Real max; //!< Max value of function
  Real a;   //!< First root where function is zero
  Real b;   //!< Second root where function is zero

public:
  //! \brief Constructor initializing the function parameters.
  QuadraticYFunc(Real MAX, Real A, Real B) : max(MAX), a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return max == Real(0); }

  //! \brief Returns first-derivative of the function.
  Real deriv(const Vec3& X, int dir) const override;
  //! \brief Returns second-derivative of the function.
  Real dderiv(const Vec3&, int dir1, int dir2) const override;

protected:
  //! \brief Evaluates the quadratic function.
  Real evaluate(const Vec3& X) const override;
};


/*!
  \brief A scalar-valued spatial function, quadratic in \a z.
*/

class QuadraticZFunc : public RealFunc
{
  Real max; //!< Max value of function
  Real a;   //!< First root where function is zero
  Real b;   //!< Second root where function is zero

public:
  //! \brief Constructor initializing the function parameters.
  QuadraticZFunc(Real MAX, Real A, Real B) : max(MAX), a(A), b(B) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return max == Real(0); }

  //! \brief Returns first-derivative of the function.
  Real deriv(const Vec3& X, int dir) const override;
  //! \brief Returns second-derivative of the function.
  Real dderiv(const Vec3&, int dir1, int dir2) const override;

protected:
  //! \brief Evaluates the quadratic function.
  Real evaluate(const Vec3& X) const override;
};


/*!
  \brief A scalar-valued spatial function, defining a rotation about the Z-axis.
  \details The time component of the function argument multiplied with the
  function parameter \a A, is interpreted as the angle of rotation (in radians)
  about the global Z-axis passing through the point \a x0, \a y0.
  The function then returns the translation in either X- or Y-direction
  (depending on the \a retX argument to the constructor) of the global point
  { \a X.x, \a X.y } corresponding to this rotation.
  \note If the function is passed a Vec3 object as argument (and not a Vec4),
  it will always return zero.
*/

class LinearRotZFunc : public RealFunc
{
  bool rX; //!< Flag telling whether to return the X- (true) or Y-component
  Real A;  //!< Magnitude of the rotation
  Real x0; //!< Global X-coordinate of rotation centre
  Real y0; //!< Global Y-coordinate of rotation centre

public:
  //! \brief Constructor initializing the function parameters.
  LinearRotZFunc(bool retX, Real a, Real x_0 = Real(0), Real y_0 = Real(0))
    : rX(retX), A(a), x0(x_0), y0(y_0) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return A == Real(0); }
  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return this->isZero(); }

  //! \brief Returns first-derivative of the function.
  Real deriv(const Vec3&, int dir) const override;

protected:
  //! \brief Evaluates the rotation function.
  Real evaluate(const Vec3& X) const override;
};


/*!
  \brief A scalar-valued spatial function, step in \a x.
*/

class StepXFunc : public RealFunc
{
  Real fv; //!< The non-zero function value
  Real x0; //!< The function is zero for \a x < \a x0
  Real x1; //!< The function is zero for \a x > \a x1
  char d;  //!< Coordinate to step in (default X).

public:
  //! \brief Constructor initializing the function parameters.
  explicit StepXFunc(Real v, Real a = Real(0), Real b = Real(1), char dir = 'X')
    : fv(v), x0(a), x1(b), d(dir) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return fv == Real(0); }

protected:
  //! \brief Evaluates the step function.
  Real evaluate(const Vec3& X) const override;
};


/*!
  \brief A scalar-valued spatial function, step in \a x and \a y.
*/

class StepXYFunc : public RealFunc
{
  Real fv; //!< The non-zero function value
  Real x0; //!< The function is zero for \a x < \a x0
  Real y0; //!< The function is zero for \a y < \a y0
  Real x1; //!< The function is zero for \a x > \a x1
  Real y1; //!< The function is zero for \a y > \a y1

public:
  //! \brief Constructor initializing the function parameters.
  explicit StepXYFunc(Real v,
                      Real X1 = Real( 1), Real Y1 = Real( 1),
                      Real X0 = Real(-1), Real Y0 = Real(-1))
    : fv(v), x0(X0), y0(Y0), x1(X1), y1(Y1) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return fv == Real(0); }

protected:
  //! \brief Evaluates the step function.
  Real evaluate(const Vec3& X) const override;
};


/*!
  \brief A scalar-valued spatial function, linear interpolation.
*/

class Interpolate1D : public RealFunc
{
  LinearFunc lfunc; //!< The piece-wise linear function doing the interpolation
  int        dir;   //!< In which direction to perform the interpolation
  Real       time;  //!< Ramp-up time

public:
  //! \brief Constructor initializing the function parameters from a file.
  //! \param[in] file Name of file to read function data from
  //! \param[in] dir_ Coordinate direction of the spatial variation
  //! \param[in] col Which column of the file to read function values from
  //! \param[in] ramp Ramp-up time
  Interpolate1D(const char* file, int dir_, int col = 2, Real ramp = Real(0))
    : lfunc(file,col), dir(dir_), time(ramp) {}

  //! \brief Constructor initializing the function parameters from two lists.
  //! \param[in] xVal List of function argument values
  //! \param[in] yVal List of function values
  //! \param[in] dir_ Coordinate direction of the spatial variation
  //! \param[in] ramp Ramp-up time
  Interpolate1D(const std::vector<Real>& xVal, const std::vector<Real>& yVal,
                int dir_ = 1, Real ramp = Real(0))
    : lfunc(xVal,yVal), dir(dir_), time(ramp) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return lfunc.isZero(); }
  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return time <= Real(0); }

  //! \brief Returns first-derivative of the function.
  Real deriv(const Vec3&, int ddir) const override;

protected:
  //! \brief Evaluates the function by interpolation.
  Real evaluate(const Vec3& X) const override;
};


/*!
  \brief A vector-valued spatial function, constant in space and time.
*/

class ConstVecFunc : public VecFunc
{
  Vec3 fval; //!< The function value

public:
  //! \brief Constructor initializing the function value.
  explicit ConstVecFunc(const Vec3& v) : fval(v) {}

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return fval.isZero(0.0); }

protected:
  //! \brief Evaluates the constant function.
  Vec3 evaluate(const Vec3&) const override { return fval; }
};


namespace utl
{
  //! \brief Creates a scalar-valued function by parsing a character string.
  const RealFunc* parseRealFunc(char* cline, Real A = Real(1),
                                bool print = true);

  //! \brief Creates a time function by parsing a character string.
  //! \param[in] func Character string to parse function definition from
  //! \param[in] type %Function definition type flag
  //! \param[in] eps Domain increment for calculation of numerical derivative
  ScalarFunc* parseTimeFunc(const char* func,
                            const std::string& type = "expression",
                            Real eps = Real(1.0e-8));

  //! \brief Creates a vector-valued time function by parsing a char string.
  //! \param[in] func Character string to parse function definition from
  //! \param[in] type %Function definition type flag
  VecTimeFunc* parseVecTimeFunc(const char* func,
                                const std::string& type);

  //! \brief Creates a scalar-valued function by parsing a character string.
  //! \param[in] func Character string to parse function definition from
  //! \param[in] type %Function definition type flag
  //! \param[in] print If \e true, print function definition
  RealFunc* parseRealFunc(const std::string& func,
                          const std::string& type = "expression",
                          bool print = true);

  //! \brief Creates a scalar-valued function by parsing a character string.
  //! \param[in] function %Function expression
  //! \param[in] autodiff if \e true, auto-differentiation is enabled
  //!
  //! \details The implementation of this method is in the file ExprFunctions.C
  //! for encapsulation of the autodiff package.
  RealFunc* parseExprRealFunc(const std::string& function, bool autodiff);

  //! \brief Creates a scalar-valued int function by parsing a character string.
  //! \param[in] func Character string to parse function definition from
  //! \param[in] type %Function definition type flag
  IntFunc* parseIntFunc(const std::string& func,
                        const std::string& type = "expression");

  //! \brief Creates a vector-valued function by parsing a character string.
  //! \param[in] func Character string to parse function definition from
  //! \param[in] type %Function definition type flag
  //! \param[in] variables Optional variable definition for expression functions
  VecFunc* parseVecFunc(const std::string& func,
                        const std::string& type = "expression",
                        const std::string& variables = "");

  //! \brief Creates a tensor-valued function by parsing a character string.
  //! \param[in] func Character string to parse function definition from
  //! \param[in] type %Function definition type flag
  TensorFunc* parseTensorFunc(const std::string& func,
                              const std::string& type);

  //! \brief Creates a vector-valued function defining a surface traction.
  //! \param[in] func Character string to parse function definition from
  //! \param[in] type %Function definition type flag
  //! \param[in] dir Coordinate direction of the traction (0=normal direction)
  TractionFunc* parseTracFunc(const std::string& func,
                              const std::string& type = "expression",
                              int dir = 0);

  //! \brief Creates a vector-valued function defining a surface traction.
  //! \param[in] elem Pointer to XML-element to parse function definition from
  TractionFunc* parseTracFunc(const tinyxml2::XMLElement* elem);
}

#endif
