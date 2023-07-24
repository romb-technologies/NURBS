#ifndef NURBS_H
#define NURBS_H

#include <map>
#include <memory>
#include <iostream>
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
   * \brief Raise the curve order by 1
   *
   * Curve will always retain its shape
   * \warning Resets cached data
   */
  void elevateOrder();

  /*!
   * \brief Lower the curve order by 1
   *
   * If current shape cannot be described by lower order, it will be best aproximation
   * \warning CAN THROW: Cannot be called for curves of 1st order
   * \warning Resets cached data
   */
  void lowerOrder();

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

  /*!
   * \brief Set the new coordinates to a control point
   * \param idx Index of chosen control point
   * \param point New control point
   */
  void setControlPoint(unsigned idx, const Point& point);

  /*!
   * \brief Get first and last control points
   * \return A pair of end points
   */
  std::pair<Point, Point> endPoints() const;

  /*!
   * \brief Reverse order of control points
   */
  void reverse();

  /*!
   * \brief Get a polyline representation of the curve as a vector of points on curve
   * \param flatness Error tolerance of approximation
   * \return A vector of polyline vertices
   */
  PointVector polyline(double flatness = 0.5) const;

  /*!
   * \brief Get the point on curve for a given t
   * \param t Curve parameter
   * \return Point on a curve for a given t
   */
  Point valueAt(double t) const;

  /*!
   * \brief Get the point on curve for a given t
   * \param t Curve parameter
   * \return Point on a curve for a given t
   */
  Point valueAt2(double t) const;

  /*!
   * \brief Get the point vector on curve for given parameters
   * \param t_vector Curve parameters
   * \return Matrix of points on a curve for given parameters
   */
  Eigen::MatrixX2d valueAt(const std::vector<double>& t_vector) const;

  Eigen::MatrixXd basisFunction2(int i, int k) const;

  /*!
   * \brief Get the bounding box of curve
   * \return Bounding box (if use_roots is false, returns the bounding box of control points)
   */
  BoundingBox boundingBox() const;

  /*!
   * \brief Get the derivative of a curve
   * \return Derivative curve
   */
  const Curve& derivative() const;

  /*!
   * \brief Get the nth derivative of a curve
   * \param n Desired number of derivative
   * \return Derivative curve
   * \warning double n cannot be zero
   */
  const Curve& derivative(unsigned n) const;

  /*!
   * \brief Get value of a derivative for a given t
   * \param t Curve parameter
   * \return Curve derivative at t
   */
  Vector derivativeAt(double t) const;

  /*!
   * \brief Get value of an nth derivative for a given t
   * \param n Desired number of derivative
   * \param t Curve parameter
   * \return nth curve derivative at t
   */
  Vector derivativeAt(unsigned n, double t) const;

  /*!
   * \brief Get curvature of the curve for a given t
   * \param t Curve parameter
   * \return Curvature of a curve for a given t
   */
  double curvatureAt(double t) const;

  /*!
   * \brief Get curvature derivative of the curve for a given t
   * \param t Curve parameter
   * \return Curvature derivative of a curve for a given t
   */
  double curvatureDerivativeAt(double t) const;

  /*!
   * \brief Get the tangent of the curve for a given t
   * \param t Curve parameter
   * \param normalize If the resulting tangent should be normalized
   * \return Tangent of a curve for a given t
   */
  Vector tangentAt(double t, bool normalize = true) const;

  /*!
   * \brief Get the normal of the curve for a given t
   * \param t Curve parameter
   * \param normalize If the resulting normal should be normalized
   * \return Normal of a curve for given t
   */
  Vector normalAt(double t, bool normalize = true) const;

  /*!
   * \brief Get the parameter t where curve is closest to given point
   * \param point Point to project on curve
   * \return double t
   */
  double projectPoint(const Point& point) const;

protected:
  /*!
   * \brief N x 2 matrix where each row corresponds to control Point
   * \warning Any changes made to control_points_ require a call to resetCache() funtion!
   */
  Eigen::MatrixX2d control_points_;

  using BasisFunctionsMap = std::map<unsigned, Eigen::MatrixXd>;

  /// Reset all privately cached data
  inline void resetCache();

private:
  /// Number of control points (order + 1)
  unsigned N_{}, p_{};

  mutable std::unique_ptr<const Curve> cached_derivative_;    /*! If generated, stores derivative for later use */
  mutable std::unique_ptr<std::vector<double>> cached_roots_; /*! If generated, stores roots for later use */
  mutable std::unique_ptr<BoundingBox> cached_bounding_box_;  /*! If generated, stores bounding box for later use */
  mutable std::unique_ptr<PointVector> cached_polyline_;      /*! If generated, stores polyline for later use */
  mutable double cached_polyline_flatness_{};                 /*! Flatness of cached polyline */


  mutable std::unique_ptr<Eigen::VectorXd> knot_vector_;
  mutable unsigned m_{};

  Eigen::VectorXd knotVector() const;
  int getKnotSpanIndex(double u, int p) const;
  Eigen::VectorXd getBasisFunctions(int i, double u, int p) const;
  Eigen::VectorXd getDerivativeBasisFunctions(int i, double u, int p, int n) const;
};

}

#endif // NURBS_H
