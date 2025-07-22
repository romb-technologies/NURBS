#ifndef NURBS_H
#define NURBS_H

#include <limits>
#include <map>
#include <memory>
#include <numeric>

#include "NURBS/declarations.h"
#include "NURBS/span.h"
#include "NURBS/utils.h"

#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/Polynomials>

namespace NURBS
{

/*!
 * \brief A NURBS curve class
 *
 * A class for storing and using any-order NURBS curve.
 * It uses private and static caching for storing often accessed data.
 * Private caching is used for data concerning individual curves, while
 * static caching is used for common data (coefficient matrices)
 */
class Curve
{
public:
  ~Curve() = default;

  /*!
   * \brief Create a NURBS curve of order \c p based on an array of control points.
   * \param points Nx2 matrix where each row is one of N control points that define the curve.
   * \param p Order of the curve, defaults to 3.
   *
   * Creates a new NURBS curve based on a Nx2 matrix of points.
   * Its weights are all set to 1 by default, \c N is automatically set to be the
   * number of points in the array. \c N must be higher or equal to \c p+1.
   * The knot vector is automatically set to consist of equally spaced knots.
   */
  Curve(Eigen::MatrixX2d points, int p = 3);

  /*!
   * \brief Create a NURBS curve of order \c p based on an array of control points.
   * \param points A std::vector of N control points that define the curve.
   * \param p Order of the curve, defaults to 3.
   *
   * Creates a new NURBS curve based on a std::vector of points.
   * Its weights are all set to 1 by default, \c N is automatically set to be the
   * number of points in the vector. \c N must be higher or equal to \c p+1.
   * The knot vector is automatically set to consist of equally spaced knots.
   */
  Curve(const PointVector& points, int p = 3);

  /*!
   * \brief Create a NURBS curve with defined weighted control points, knot vector and order.
   * \param points Nx3 matrix where each row is one of N weighted control points that define the curve.
   * \param knotvector Array of knots in ascending order, must be of length N+p+1.
   * \param p Order of the curve, defaults to 3.
   */
  Curve(Eigen::MatrixX3d wpoints, Eigen::ArrayXd knotvector, int p = 3);

  Curve(const Curve& curve);
  Curve(Curve&&) = default;

  Curve& operator=(const Curve&);
  Curve& operator=(Curve&&) = default;

  /*!
   * \brief Get the order of the curve;
   * \return Order of the curve
   */
  unsigned order() const;

  /*!
   * \brief Raise the curve order by \c t
   * \param t Number of orders, defaults to 1.
   *
   * Curve will always retain its shape
   * \warning Resets cached data
   */
  void elevateOrder(uint t = 1);

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
   * \brief Get the control point at index \c idx
   * \param idx Index of chosen control point
   * \return Control point
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
   * \brief Reverse order of control points in place
   */
  void reverse();

  /*!
   * \brief Get a polyline representation of the curve as a vector of points on curve
   * \param flatness Error tolerance of approximation
   * \return A vector of polyline vertices
   */
  PointVector polyline(double flatness = 0.5) const;

  /*!
   * \brief Get the point on this curve for a given \c t
   * \param t Curve parameter
   * \return Point on a curve for a given \c t
   */
  Point valueAt(double t) const;

  /*!
   * \brief Get the point vector on curve for given parameters
   * \param t_vector Curve parameters
   * \return Matrix of points on a curve for given parameters
   */
  Eigen::MatrixX2d valueAt(const std::vector<double>& t_vector) const;

  /*!
   * \brief Get the bounding box of curve
   * \return Bounding box (if use_roots is false, returns the bounding box of control points)
   */
  BoundingBox boundingBox(bool use_roots = true) const;

  /*!
   * \brief Get the derivative of a curve
   * \return Derivative curve
   */
  const Curve& derivative() const;

  /*!
   * \brief Get the \c n-th derivative of a curve
   * \param n Desired number of derivative
   * \return Derivative curve
   * \warning double n cannot be zero
   */
  const Curve& derivative(unsigned n) const;

  /*!
   * \brief Get value of a derivative for a given \c t
   * \param t Curve parameter
   * \return Curve derivative at \c t
   */
  Vector derivativeAt(double t) const;

  /*!
   * \brief Get value of an \c n-th derivative for a given \c t
   * \param n Desired number of derivative
   * \param t Curve parameter
   * \return \c n-th curve derivative at \c t
   */
  Vector derivativeAt(unsigned n, double t) const;

  /*!
   * \brief Get roots of the curve on both axes
   * \return A vector of parameters where curve passes through axes
   */
  std::vector<double> roots() const;

  /*!
   * \brief Get all extrema of the curve
   * \return A vector of parameters where extrema are
   */
  std::vector<double> extrema() const;

  /*!
   * \brief Get curvature of the curve for a given \c t
   * \param t Curve parameter
   * \return Curvature of a curve for a given \c t
   */
  double curvatureAt(double t) const;

  /*!
   * \brief Get curvature derivative of the curve for a given \c t
   * \param t Curve parameter
   * \return Curvature derivative of a curve for a given \c t
   */
  double curvatureDerivativeAt(double t) const;

  /*!
   * \brief Get the tangent of the curve for a given \c t
   * \param t Curve parameter
   * \param normalize If the resulting tangent should be normalized
   * \return Tangent of a curve for a given \c t
   */
  Vector tangentAt(double t, bool normalize = true) const;

  /*!
   * \brief Get the normal of the curve for a given \c t
   * \param t Curve parameter
   * \param normalize If the resulting normal should be normalized
   * \return Normal of a curve for given \c t
   */
  Vector normalAt(double t, bool normalize = true) const;

  /*!
   * \brief Get the parameter \c t where curve is closest to given point
   * \param point Point to project on curve
   * \return double \c t
   */
  double projectPoint(const Point& point) const;

  /*!
   * \brief Get the intersections between this curve and another curve
   * \param curve Second curve
   * \return Vector of points where the curves intersect
   */
  PointVector intersections(const Curve& other) const;

  /*!
   * \brief Get the weight of the control point at index \c idx
   * \param idx Weight index
   * \return Weight at index \c idx
   */
  double weight(int idx) const;

  /*!
   * \brief Set the weight of the control point at index \c idx
   * \param w New weight
   * \param idx Weight index
   */
  void setWeight(int idx, double value);

  /*!
   * \brief Get the weight vector of the curve
   * \return Vector of weights
   */
  Eigen::VectorXd weights() const;

  /*!
   * \brief Get the knot vector of the curve
   * \return Array of knots
   */
  Eigen::ArrayXd knotVector() const;

  /*!
   * \brief Set the knot at index \c idx
   * \param value New knot value
   * \param idx Knot index
   */
  void setKnot(int idx, double value);

  /*!
   * \brief Get the knot at index \c idx
   * \param idx Knot index
   * \return Knot value
   */
  double knot(int idx) const;

  /*!
   * \brief Append a new control point to the end of this curve
   * \param point New point
   */
  void appendPoint(Point point);

  /*!
   * \brief Insert a new knot into this curve at parameter \c t
   * \param t New knot
   * \param s Multiplicity of knot
   * \param r Number of insertions
   */
  void insertKnot(double t, int r);

  /*!
   * \brief Split this curve into two new curves at parameter \c t
   * \param t Where to split the curve
   * \return Pair of new curves
   */
  std::pair<Curve, Curve> splitCurve(double t) const;

  /*!
   * \brief Create a series of Bezier curves that is equivalent to this curve
   * \return Vector of Bezier curves
   */
  std::vector<Curve> piecewiseBezier() const;

  /*!
   * \brief Evaluate the basis functions of this curve at parameter \c t
   * \param t Curve parameter
   * \return Vector of basis function values
   */
  Eigen::VectorXd getBasisFunctionsAt(double t) const;

  /*!
   * \brief Get the index of the span in which parameter \c t lies
   * \param t Curve parameter
   * \return Span index
   */
  int getKnotSpanIndex(double t) const;

  /*!
   * \brief Evaluate the length of this curve
   * \return Curve length
   */
  double length() const;

  /*!
   * \brief Evaluate the length of this curve from its starting point
   * up to the point at parameter \c t
   * \return Curve length up to \c t;
   */
  double length(double t) const;

  /*!
   * \brief Remove the \c ix-th knot of this curve \c k times
   * \param ix Index of knot to be removed
   * \param k Number of times to remove the knot (must be less or equal to knot's multiplicity, defaults to 1)
   */
  void removeKnot(int ix, int k = 1);

  /*!
   * \brief Remove the \c ix-th control point of this curve
   * \param ix Index of control point to be removed
   */
  void removeControlPoint(int ix);

  /*!
   * \brief Join this curve with another curve to create a new curve
   * \param other Reference to second curve
   * \return New curve that is a result of joining the two curves
   */
  Curve join(Curve& other);

  /*!
   * \brief Make this curve continuous with \c source_curve
   * \param source_curve
   * \return New curve that is a result of joining the two curves
   */
  void applyContinuity(const Curve& source_curve, const std::vector<double>& beta_coeffs);

protected:
  /*!
   * \brief N x 3 matrix where each row corresponds to a weighted control Point and its weight
   * \warning Any changes made to weighted_control_points_ require a call to resetCache() funtion!
   */
  Eigen::MatrixX3d weighted_control_points_;

  /*!
   * \brief Reset all privately cached data
   */
  inline void resetCache();

private:
  /// Number of control points
  unsigned N_{};
  /// Order of curve
  unsigned p_{};

  mutable std::optional<std::vector<double>> cached_roots_; /*! If generated, stores roots for later use */
  mutable std::optional<BoundingBox> cached_bounding_box_;  /*! If generated, stores bounding box for later use */
  mutable std::optional<PointVector> cached_polyline_;      /*! If generated, stores polyline for later use */
  mutable double cached_polyline_flatness_{};               /*! Flatness of cached polyline */

  /// Knot vector
  Eigen::ArrayXd T_;
  /// Knot spans
  mutable std::vector<Span> spans_;

  Span& getKnotSpan(double t) const;
  int getKnotMultiplicity(double t) const;
  void normalizeKnotVector();
};

} // namespace NURBS

#endif // NURBS_H
