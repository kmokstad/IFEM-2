// $Id$
//==============================================================================
//!
//! \file Function.h
//!
//! \date Jul 1 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General functions with arbitrary argument and value type.
//!
//==============================================================================

#ifndef UTL_FUNCTION_H
#define UTL_FUNCTION_H

#include "matrixnd.h"
#include "Tensor.h"
#include "Vec3.h"
#include <cstddef>


namespace utl
{
  /*!
    \brief Base class for unary functions of arbitrary result and argument type.
  */

  template<class Arg, class Result>
  class Function
  {
  protected:
    //! \brief The constructor is protected to allow sub-class instances only.
    Function() {}

  public:
    //! \brief Empty destructor.
    virtual ~Function() {}

    //! \brief Returns whether the function is identically zero or not.
    virtual bool isZero() const { return false; }
    //! \brief Returns whether the function is time-independent or not.
    virtual bool isConstant() const { return true; }

  protected:
    //! \brief Evaluates the function for the argument \a x.
    virtual Result evaluate(const Arg& x) const = 0;

  public:
    //! \brief Operator returning the function value for the given argument.
    Result operator()(const Arg& x) const { return this->evaluate(x); }

    typedef Arg    Input;  //!< Input type
    typedef Result Output; //!< Output type
  };


  /*!
    \brief Base class for binary function of arbitrary result and argument type.
    \details The two arguments have to be of the same type.
  */

  template<class Arg, class Result>
  class Function2
  {
  protected:
    //! \brief The constructor is protected to allow sub-class instances only.
    Function2() {}

  public:
    //! \brief Empty destructor.
    virtual ~Function2() {}

    //! \brief Returns whether the function is identically zero or not.
    virtual bool isZero() const { return false; }

  protected:
    //! \brief Evaluates the function for the arguments \a x and \a y.
    virtual Result evaluate(const Arg& x, const Arg& y) const = 0;

  public:
    //! \brief Operator returning the function value for the given arguments.
    Result operator()(const Arg& x, const Arg& y) const
    { return this->evaluate(x,y); }

    typedef Arg    Input;  //!< Input type
    typedef Result Output; //!< Output type
  };


  /*!
    \brief Base class for unary spatial function of arbitrary result type.
    \details Includes interfaces for evaluation of first and second derivates.
  */

  template<class Result>
  class SpatialFunction : public Function<Vec3,Result>
  {
  protected:
    //! \brief The constructor is protected to allow sub-class instances only.
    explicit SpatialFunction(const Result& val) : zero(val) {}

  public:
    //! \brief Returns a first-derivative of the function.
    virtual Result deriv(const Vec3&, int) const { return zero; }
    //! \brief Returns a second-derivative of the function.
    virtual Result dderiv(const Vec3&, int, int) const { return zero; }

  protected:
    //! \brief Returns the gradient of the function as a 1D array.
    virtual std::vector<Real> evalGradient(const Vec3&) const { return {}; }
    //! \brief Returns the hessian of the function as a 1D array.
    virtual std::vector<Real> evalHessian(const Vec3&) const { return {}; }
    //! \brief Returns the time derivatives of the function as a 1D array.
    virtual std::vector<Real> evalTimeDerivative(const Vec3&) const { return {}; }

    Result zero; //!< Return value for default implementations of derivatives
  };
}


/*!
  \brief Scalar-valued unary function of a scalar value.
*/

class ScalarFunc : public utl::Function<Real,Real>
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  ScalarFunc() {}

public:
  //! \brief Returns the first-derivative of the function.
  virtual Real deriv(Real) const { return Real(0); }
};


/*!
  \brief Base class for unary spatial functions of arbitrary result type.
  \details Includes an interface for returning the function value as an array,
  and for assigning additional state-dependent parameter values.
*/

class FunctionBase
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  FunctionBase() : ncmp(1) {}

public:
  //! \brief Empty destructor.
  virtual ~FunctionBase() {}

  //! \brief Returns the function type flag.
  virtual unsigned char getType() const = 0;

  //! \brief Returns the function value as an array.
  virtual std::vector<Real> getValue(const Vec3&) const = 0;
  //! \brief Returns a representative scalar equivalent of the function value.
  virtual Real getScalarValue(const Vec3&) const = 0;

  //! \brief Returns the number of components of the return value.
  size_t dim() const { return ncmp; }

  //! \brief Sets the active patch.
  virtual bool initPatch(size_t) { return true; }

  //! \brief Checks if a specified point is within the function domain.
  virtual bool inDomain(const Vec3&) const { return true; }

  //! \brief Sets an additional parameter in the function.
  virtual void setParam(const std::string&, double) {}
  //! \brief Sets additional parameter values in the function.
  void setParam(const std::string& name, const Vec3& value)
  {
    this->setParam(name + "x", value.x);
    this->setParam(name + "y", value.y);
    this->setParam(name + "z", value.z);
  }

protected:
  size_t ncmp; //!< Number of components in the return value
};


/*!
  \brief Scalar-valued unary function of a spatial point.
*/

class RealFunc : public utl::SpatialFunction<Real>, public FunctionBase
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  RealFunc() : utl::SpatialFunction<Real>(Real(0)) {}

public:
  //! \brief Returns the function type flag.
  unsigned char getType() const override { return 1; }

  //! \brief Returns the function value as an array.
  std::vector<Real> getValue(const Vec3& X) const override
  {
    return std::vector<Real>(1,this->evaluate(X));
  }

  //! \brief Evaluates first derivatives of the function.
  virtual Vec3 gradient(const Vec3& X) const
  {
    Vec3 result;
    for (size_t d = 1; d <= 3; ++d)
      result[d-1] = this->deriv(X,d);

    return result;
  }

  //! \brief Evaluates second derivatives of the function.
  virtual SymmTensor hessian(const Vec3& X) const
  {
    SymmTensor result(3);
    for (size_t d1 = 1; d1 <= 3; ++d1)
      for (size_t d2 = d1; d2 <= 3; ++d2)
        result(d1,d2) = this->dderiv(X, d1, d2);

    return result;
  }

  //! \brief Returns a representative scalar equivalent of the function value.
  Real getScalarValue(const Vec3& X) const override { return this->evaluate(X); }

  //! \brief Returns the time derivative of the function.
  Real timeDerivative(const Vec3& X) const { return this->deriv(X,4); }
};


/*!
  \brief Vector-valued unary function of a spatial point.
*/

class VecFunc : public utl::SpatialFunction<Vec3>, public FunctionBase
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  explicit VecFunc(size_t n = 3) : utl::SpatialFunction<Vec3>(Vec3())
  {
    ncmp = n;
  }

public:
  //! \brief Returns the function type flag.
  unsigned char getType() const override { return 2; }

  //! \brief Returns the function value as an array.
  std::vector<Real> getValue(const Vec3& X) const override
  {
    return this->evaluate(X).vec(ncmp);
  }

  //! \brief Returns a representative scalar equivalent of the function value.
  Real getScalarValue(const Vec3& X) const override
  {
    return this->evaluate(X).length();
  }

  //! \brief Evaluates first derivatives of the function.
  Tensor gradient(const Vec3& X) const
  {
    Tensor result(ncmp);
    result = this->evalGradient(X);
    return result;
  }

  //! \brief Evaluates second derivatives of the function.
  utl::matrix3d<Real> hessian(const Vec3& X) const
  {
    utl::matrix3d<Real> result(ncmp,ncmp,ncmp);
    result.fill(this->evalHessian(X).data());
    return result;
  }

  //! \brief Evaluates time derivatives of the function.
  Vec3 timeDerivative(const Vec3& X) const { return this->evalTimeDerivative(X); }
};


/*!
  \brief Vector-valued binary function of a spatial point and normal vector.
*/

class TractionFunc : public utl::Function2<Vec3,Vec3>
{
public:
  //! \brief Returns whether the traction is always normal to the face or not.
  virtual bool isNormalPressure() const { return false; }

  //! \brief Returns the time-derivative of the function.
  virtual Vec3 timeDerivative(const Vec3&, const Vec3&) const { return Vec3(); }
};

#endif
