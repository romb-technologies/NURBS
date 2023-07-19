#include "NURBS/nurbs.h"

#include <numeric>

#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/Polynomials>

using namespace NURBS;

struct _PolynomialRoots : public std::vector<double>
{
  explicit _PolynomialRoots(unsigned reserve) { std::vector<double>::reserve(reserve); }
  void clear(){}         // no-op so that PolynomialSolver::RealRoots() doesn't clear it
  void push_back(double t) // only allow valid roots
  {
    if (t >= 0 && t <= 1)
      std::vector<double>::push_back(t);
  }
};

inline unsigned _exp2(unsigned exp)
{
  return 1 << exp;
}

inline double _pow(double base, unsigned exp)
{
  double result = exp & 1 ? base : 1;
  while (exp >>= 1)
  {
    base *= base;
    if (exp & 1)
      result *= base;
  }
  return result;
}

inline Eigen::RowVectorXd _powSeries(double base, unsigned exp)
{
  Eigen::RowVectorXd power_series(exp);
  power_series(0) = 1;
  for (unsigned k = 1; k < exp; k++)
    power_series(k) = power_series(k - 1) * base;
  return power_series;
}

inline Eigen::VectorXd _trimZeroes(const Eigen::VectorXd& vec)
{
  auto idx = vec.size();
  while (idx && std::abs(vec(idx - 1)) < _epsilon)
    --idx;
  return vec.head(idx);
}

///// Curve::Curve

Curve::Curve(Eigen::MatrixX2d points) : control_points_(std::move(points)), N_(control_points_.rows()) {}

Curve::Curve(const PointVector& points)
    : control_points_(Eigen::Index(points.size()), Eigen::Index(2)), N_(points.size())
{
  for (unsigned k = 0; k < N_; k++)
    control_points_.row(k) = points[k];
}

Curve::Curve(const Curve& curve) : Curve(curve.control_points_) {}

Curve& Curve::operator=(const Curve& curve)
{
  control_points_ = curve.control_points_;
  resetCache();
  return *this;
}

PointVector Curve::controlPoints() const
{
  PointVector points(N_);
  for (unsigned k = 0; k < N_; k++)
    points[k] = control_points_.row(k);
  return points;
}

Point Curve::controlPoint(unsigned idx) const { return control_points_.row(idx); }

void Curve::setControlPoint(unsigned idx, const Point& point)
{
  control_points_.row(idx) = point;
  resetCache();
}

PointVector Curve::polyline(double flatness) const
{
  if (!cached_polyline_ || cached_polyline_flatness_ != flatness)
  {
      cached_polyline_flatness_ = flatness;
      cached_polyline_ = std::make_unique<PointVector>();
      for(float t=0; t<=1; t+=0.1) {
          cached_polyline_->emplace_back(valueAt(t));
      }
  }
  return *cached_polyline_;
}

Point Curve::valueAt(double t) const
{
  if (N_ == 0)
    return {0, 0};
  return (_powSeries(t, N_) * bernsteinCoeffs(N_) * control_points_).transpose();
}

BoundingBox Curve::boundingBox() const
{
  if (!cached_bounding_box_)
  {
    cached_bounding_box_ = std::make_unique<BoundingBox>(Point(control_points_.col(0).minCoeff(), control_points_.col(1).minCoeff()),
                                                         Point(control_points_.col(0).maxCoeff(), control_points_.col(1).maxCoeff()));
  }
  return *cached_bounding_box_;
}

void Curve::resetCache()
{
  N_ = control_points_.rows();
  cached_derivative_.reset();
  cached_roots_.reset();
  cached_bounding_box_.reset();
  cached_polyline_.reset();
}

int Curve::getKnotSpanIndex(u, p) {
    if(u == knot_vector_[N_+1])
        return N_;
    int low = p;
    int high = N_ + 1;
    int mid = (low + high)/2;
    while (u < knot_vector_[mid] || u >= knot_vector_[mid+1])
    {
        if (u < knot_vector_[mid]) high = mid;
        else low = mid;
        mid = (low+high)/2;
    {
    return(mid);
}

