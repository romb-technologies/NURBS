
#ifndef NURBS_H
#define NURBS_H

#include "declarations.h"

namespace NURBS {


class Curve
{
public:
  ~Curve() = default;

  /*!
   * \brief Create a NURBS curve
   * \param points Nx2 matrix where each row is one of N control points that define the curve
   */
  Curve(Eigen::MatrixX2d points);

  /*!
   * \brief Create a NURBS curve
   * \param points A vector of control points that define the curve
   */
  Curve(const PointVector& points);

  Curve(const Curve& curve);
  Curve(Curve&&) = default;
  Curve& operator=(const Curve&);
  Curve& operator=(Curve&&) = default;


  /*!
   * \brief Get order of the curve (Nth order curve is described with N+1 points);
   * \return Order of curve
   */
  unsigned order() const;

  /*!
   * \brief Get a vector of control points
   * \return A vector of control points
   */
  PointVector controlPoints() const;

  /*!
   * \brief Get the control point at index idx
   * \param idx Index of chosen control point
   * \return A vector of control points
   */
  Point controlPoint(unsigned idx) const;
};

}

#endif // NURBS_H