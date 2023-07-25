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
  Eigen::RowVectorXd power_series(exp+1);
  power_series(0) = 1;
  for (unsigned k = 1; k <= exp; k++)
    power_series(k) = power_series(k - 1) * base;
  return power_series;
}

inline Eigen::RowVectorXd _powSeriesDerivative(double base, unsigned exp, unsigned drv)
{
  Eigen::RowVectorXd power_series = Eigen::RowVectorXd::Ones(exp+1);
  for (int i=0; i<drv-1; i++) {
    for (int j=0; j<=exp; j++) {
        power_series(j) *= j-i;
    }
  }
  for (unsigned k = drv; k <= exp; k++)
    power_series(k) *= pow(base, k-drv+1);
  return power_series;
}

inline Eigen::MatrixXd _slice(Eigen::MatrixXd base, int start, int end)
{
    Eigen::MatrixXd out(end-start, base.cols());
    for (int i=start; i<end; i++) {
        out.row(i-start) = base.row(i);
    }
    return out;
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
  *knot_vector_ = *curve.knot_vector_;
  p_ = curve.p_;
  resetCache();
  return *this;
}

unsigned Curve::order() const {
    return p_;
}

void Curve::elevateOrder() {
    p_++;
    resetCache();
}

void Curve::lowerOrder() {
    if (p_>0) p_--;
    resetCache();
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

void Curve::reverse() {
    control_points_ = control_points_.colwise().reverse().eval();
    resetCache();
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
  auto& kv = *knot_vector_;
  double u = (t-kv(i))/(kv(i+1)-kv(i));

  return _powSeries(u, p_) * basisFunction2(i, p_+1) * _slice(control_points_, i-p_, i+1);
}


Eigen::MatrixX2d Curve::valueAt(const std::vector<double>& t_vector) const
{
    Eigen::MatrixX2d c(t_vector.size(), 2);
    if (N_ == 0)
        return {0, 0};
    for (int k=0; k<t_vector.size(); k++) {
        int i = getKnotSpanIndex(t_vector[k], p_);
        Eigen::VectorXd funs = getBasisFunctions(i, t_vector[k], p_);
        Point c_(0, 0);
        for (uint j=0; j<=p_; j++) {
            c_ += funs[j] * control_points_.row(i+j-p_);
        }
        c.row(k) = c;
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

const Curve& Curve::derivative() const
{

}


const Curve& Curve::derivative(unsigned n) const {

}

Vector Curve::derivativeAt(unsigned n, double t) const
{
  if (N_ == 0)
    return {0, 0};

  int i = getKnotSpanIndex(t, p_);
  auto& kv = *knot_vector_;
  double u = (t-kv(i))/(kv(i+1)-kv(i));

  return _powSeriesDerivative(u, p_, 2) * basisFunction2(i, p_+1) * _slice(control_points_, i-p_, i+1);
}

Vector Curve::derivativeAt(double t) const
{
    return derivativeAt(2, t);
}

double Curve::curvatureAt(double t) const
{
  Vector d1 = derivativeAt(t);
  Vector d2 = derivativeAt(2, t);

  return (d1.x() * d2.y() - d1.y() * d2.x()) / _pow(d1.norm(), 3);
}

double Curve::curvatureDerivativeAt(double t) const
{
  Vector d1 = derivativeAt(t);
  Vector d2 = derivativeAt(2, t);
  Vector d3 = derivativeAt(3, t);

  return (d1.x() * d3.y() - d1.y() * d3.x()) / _pow(d1.norm(), 3) -
         3 * d1.dot(d2) * (d1.x() * d2.y() - d1.y() * d2.x()) / _pow(d1.norm(), 5);
}

Vector Curve::tangentAt(double t, bool normalize) const {
    Vector p(derivativeAt(t));
    if (normalize && p.norm() > 0)
      p.normalize();
    return p;
}


Vector Curve::normalAt(double t, bool normalize) const
{
  Vector tangent = tangentAt(t, normalize);
  return {-tangent.y(), tangent.x()};
}

//todo: napravit bolje
double Curve::projectPoint(const Point &point) const
{
    std::pair<double, double> min_point(0.0, (point - valueAt(0.0)).norm());

    for (double t=0.0; t<=1.0; t+=0.01) {
        double dist = (point - valueAt(t)).norm();
        min_point = dist < min_point.second ? std::make_pair(t, dist) : min_point;
    }
    for (double t = std::max(min_point.first-0.01, 0.0);
         t < min_point.first+0.01 && t<=1.0;
         t += 0.0001) {
        double dist = (point - valueAt(t)).norm();
        min_point = dist < min_point.second ? std::make_pair(t, dist) : min_point;
    }

    return min_point.first;
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
        // Default is uniform knot vector
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
    if(u <= knotVector()[0])
        return p;
    if(u >= knotVector()[N_])
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

Eigen::MatrixXd Curve::basisFunction2(int i, int k) const {
    if (k==1) {
        Eigen::VectorXd m(1); m<<1;
        return m;
    }
    else {
        Eigen::MatrixXd m(k-1, k-1);
        m = basisFunction2(i, k-1);
        auto& t = *knot_vector_;

        Eigen::MatrixXd m1(k, k-1), m2 = Eigen::MatrixXd::Zero(k-1, k), m3(k, k-1), m4 = Eigen::MatrixXd::Zero(k-1, k);

        m1 << m, Eigen::MatrixXd::Zero(1, k-1);
        m3 << Eigen::MatrixXd::Zero(1, k-1), m;

        for (int j=0; j<k-1; j++) {
            int temp = i-(k-2-j);

            m2(j, j) = 1 - (
                        t[i]-t[temp])
                        /(t[temp+k-1]-t[temp]
                    );
            m2(j, j+1) = (t[i]-t[temp])
                        /(t[temp+k-1]-t[temp]);

            m4(j, j) = -(t[i+1]-t[i])
                        /(t[temp+k-1]-t[temp]
                    );
            m4(j, j+1) = (t[i+1]-t[i])
                        /(t[temp+k-1]-t[temp]);
        }
        Eigen::MatrixXd out(k, k);
        return (m1*m2) + (m3*m4);
    }
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

Eigen::VectorXd Curve::getDerivativeBasisFunctions(int i, double u, int p, int n) const {
    Eigen::MatrixXd ndu(p+1, p+1), ders(p+1, p+1), a(p+1, p+1);
    Eigen::VectorXd left(p+1), right(p+1);
    ndu(0, 0) = 1.0;

    for (int j=1; j<=p; j++) {
        left[j] = u - (*knot_vector_)(i+1-j);
        right[j] = (*knot_vector_)(i+j) - u;
        double saved = 0.0;
        for (int r=0; r<j; r++) {
            ndu(j, r) = right(r+1) + left(j-r);
            double temp = ndu(r, j-1) / ndu(j, r);

            ndu(r, j) = saved + right(r+1) * temp;
            saved = left(j-r) * temp;
        }
        ndu(j, j) = saved;
    }

    for (int j=0; j<=p; j++) {
        ders(0, j) = ndu(j, p);
    }
    for (int r=0; r<=p; r++) {
        unsigned s1(0), s2(1), j1, j2;
        a(0, 0) = 1.0;
        for (int k=1; k<=n; k++) {
            double d = 0.0;
            int rk(r-k), pk(p-k);
            if (r >= k) {
                a(s2, 0) = a (s1, 0) / ndu(pk+1, rk);
                d = a(s2, 0) * ndu(rk, pk);
            }
            if (rk >= -1) j1 = 1;
                else j1 = -rk;
            if (r-1 <= pk) j2 = k-1;
                else j2 = p-r;
            for (int j=j1; j<=j2; j++) {
                a(s2, j) = (a(s1, j)-a(s1, j-1))/ndu(pk+1, rk+j);
                d += a(s2, j)*ndu(rk+j, pk);
            }
            if (r <= pk) {
                a(s2, k) = -a(s1, k-1)/ndu(pk+1, r);
                d += a(s2, k)*ndu(r, pk);
            }
            ders(k, r) = d;
            unsigned j=s1; s1=s2; s2=j;
        }
    }

    unsigned r = p;
    for (int k=1; k<=n; k++) {
        for (int j=0; j<=p; j++) ders(k, j) *= r;
        r *= (p-k);
    }

    return ders;
}

