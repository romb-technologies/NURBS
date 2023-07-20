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

Curve::Curve(Eigen::MatrixX2d points)
    : control_points_(std::move(points))
    , N_(control_points_.rows())
    , p_(2) {}

Curve::Curve(const PointVector& points)
    : control_points_(Eigen::Index(points.size()), Eigen::Index(2))
    , N_(points.size())
    , p_(2)
{
  for (unsigned k = 0; k < N_; k++)
    control_points_.row(k) = points[k];
}

Curve::Curve(const Curve& curve)
    : Curve(curve.control_points_) {}

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

std::pair<Point, Point> Curve::endPoints() const
{
    return {control_points_.row(0), control_points_.row(N_ - 1)};
}

PointVector Curve::polyline(double flatness) const
{
  if (!cached_polyline_ || cached_polyline_flatness_ != flatness)
  {
      cached_polyline_flatness_ = flatness;
      cached_polyline_ = std::make_unique<PointVector>();
      for(float t=0.0; t<=1; t+=0.01) {
          cached_polyline_->emplace_back(valueAt(t));
      }
  }
  return *cached_polyline_;
}

Point Curve::valueAt(double t) const
{
  if (N_ == 0)
    return {0, 0};
  int i = getKnotSpanIndex(t, p_);
  Eigen::VectorXd funs = getBasisFunctions(i, t, p_);
  Point c(0, 0);
  for (uint j=0; j<=p_; j++) {
      c += funs[j] * control_points_.row(i+j-p_);
  }
  return c;
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

Eigen::VectorXd Curve::knotVector() const
{
    if (!knot_vector_ || m_ != N_+p_+1) {
        m_ = N_ + p_ + 1;
        knot_vector_ = std::make_unique<Eigen::VectorXd>(m_);
        auto& knot_vector = *knot_vector_;

        for (uint i=0; i<p_+1; i++) {
            knot_vector(i) = 0;
        }
        int knots = m_-2*(p_+1);
        double interval = 1.0 / (N_ - 2);
        for (uint i=0; i<knots; i++) {
            knot_vector(p_+1+i) = (i+1) * interval;
        }
        for (uint i=m_-(p_+1); i<m_; i++) {
            knot_vector(i) = 1;
        }
    }
    return *knot_vector_;
}

int Curve::getKnotSpanIndex(double u, int p) const {
    if(u == knotVector()[N_])
        return N_ - 1;
    int low = p;
    int high = N_ + 1;
    int mid = (low + high)/2;
    while (u < (*knot_vector_)[mid] || u >= (*knot_vector_)[mid+1])
    {
        if (u < (*knot_vector_)[mid]) high = mid;
        else low = mid;
        mid = (low+high)/2;
    }
    return(mid);
}

Eigen::VectorXd Curve::getBasisFunctions(int i, double u, int p) const {
    Eigen::VectorXd N(p+1), left(p+1), right(p+1);
    N[0] = 1.0;
    for (int j=1; j<=p; j++)
    {
        left[j] = u - (*knot_vector_)[i+1-j];
        right[j] = (*knot_vector_)[i+j] - u;
        double saved = 0.0;
        for (int r=0; r<j; r++)
        {
            double temp = N[r]/(right[r+1]+left[j-r]);
            N[r] = saved+right[r+1]*temp;
            saved = left[j-r]*temp;
        }
        N[j] = saved;
    }
    return N;
}

