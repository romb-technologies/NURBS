#ifndef SPAN_H
#define SPAN_H

#include <memory>

#include "declarations.h"
#include "utils.h"

#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/Polynomials>

namespace NURBS
{

/*!
 * \brief A single knot span in a NURBS curve
 *
 * A class that represents a single interval between two knots in a knot vector.
 * Used for more intuitive and efficient handling of NURBS curves,
 * by definition they are a segment of an existing NURBS curve and not a
 * separate entity from the curve.
 * It uses private caching for storing often accessed data.
 */
class Span
{
public:
  ~Span() = default;
  Span(Eigen::Ref<Eigen::MatrixX3d> wpoints_, Eigen::Ref<Eigen::ArrayXd> knot_v, uint p);

  Span(Eigen::Ref<const Eigen::MatrixXd> basis_func, Eigen::Ref<const Eigen::MatrixXd> vbf,
       Eigen::Ref<const Eigen::RowVectorXd> wbf, double start, double end);

  /// Start knot of this span.
  inline double start() const { return start_; }
  /// End knot of this span.
  inline double end() const { return end_; }
  /*!
   * \brief Check if the curve parameter \c t is contained in this span.
   * \param t Curve parameter
   * \return True if this span contains \c t
   */
  bool contains(double t) const;
  /*!
   * \brief Update the knot span's basis function
   * \warning Must be called every time this span's
   *  knot vector segment, or both knot vector and control points are changed.
   */
  void update(Eigen::Ref<Eigen::ArrayXd> knots, Eigen::Ref<Eigen::MatrixX3d> wpoints);
  /*!
   * \brief Update the knot span's \c cached_vbf and \c cached_wbf
   * \warning Must be called every time this knot span's
   *  control points and/or weights are changed. If both knot vector and control points/weights
   *  are changed use the \c update method.
   */
  void updateControlPoints(Eigen::Ref<Eigen::MatrixX3d> wpoints);
  /*!
   * \brief Retrieve this knot span's basis function in matrix form.
   * \return Basis function
   */
  Eigen::MatrixXd basisFunction() const;
  /*!
   * \brief Evaluate the polyline representation of this knot span.
   * \return A vector of polyline vertices
   */
  PointVector polyline(double flatness = 0.5) const;
  /*!
   * \brief Get the point on this knot span for a given t
   * \param u Span parameter
   * \return Point on this span for a given t
   *
   * The value \c u ranges from 0.0 to 1.0 regardless of the span's
   * position on the curve, where 0.0 is the beginning of the span
   * and 1.0 is the end of the span.
   */
  Point valueAt(double u) const;

  /*!
   * \brief Get value of an nth derivative for a given t
   * \param n Desired number of derivative
   * \param u Span parameter (ranges from 0.0 to 1.0, not start_t to end_t)
   * \return nth span derivative at u
   */
  Point derivativeAt(int n, double u) const;

  /*!
   * \brief Get value of a derivative for a given t
   * \param t Curve parameter
   * \return Curve derivative at t
   */
  Point derivativeAt(double u) const;
  /// Reset all privately cached data
  void resetCache();
  /*!
   * \brief Evaluate the length of this curve
   * \return Curve length
   */
  double length() const;
  /*!
   * \brief Evaluate the length of this curve from its starting point up to the point at parameter t
   * \return Curve length up to t;
   */
  double length(double t) const;
  /*!
   * \brief Retrieve the cached weights times basis function
   * \return Vector
   */
  Eigen::RowVectorXd cachedWBF() const;
  /*!
   * \brief Retrieve the cached weighted control points
   * times basis function
   * \return Matrix
   */
  Eigen::MatrixXd cachedVBF() const;

  /*!
   * \brief Get the bounding box of this knot span
   * \return Bounding box
   */
  BoundingBox boundingBox() const;
  std::vector<double> extrema() const;
  PointVector intersections(const Span& other) const;

private:
  mutable std::optional<double> cached_length_;
  mutable std::optional<PointVector> cached_polyline_;
  mutable std::optional<Eigen::VectorXd> cached_chebyshev_coeffs_; /*!  If generated, stores chebyshev coefficients
                                                                        for calculating the length of the curve */
  mutable std::optional<BoundingBox> cached_bounding_box_;

  /// This knot span's basis function
  Eigen::MatrixXd basis_function_;
  /// Weights times basis function
  Eigen::RowVectorXd cached_wbf_;
  /// Weighted control points times basis function
  Eigen::MatrixXd cached_vbf_;
  /*!
   * \brief Span order
   * \warning Must always match curve order
   */
  const uint p_;
  double start_, end_;

  /*!
   * \brief Get the bounding box of this knot span's implicit control points
   * \warning This bounding box is bigger than the bounds of the actual curve and
   * should not be used for precise calculations.
   * \return Bounding box
   */
  static inline BoundingBox fastBoundingBox(const Eigen::MatrixXd& vbf, const Eigen::RowVectorXd& wbf,
                                            const Eigen::MatrixXd& inverse_basis_function);
};
} // namespace NURBS

#endif // SPAN_H
