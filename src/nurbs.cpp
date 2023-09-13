#include "NURBS/nurbs.h"

#include <numeric>
#include <limits>

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
  for (uint i=0; i<drv; i++) {
    for (uint j=0; j<=exp; j++) {
        power_series(j) *= j-i;
    }
  }
  for (unsigned k = drv; k <= exp; k++)
    power_series(k) *= pow(base, k-drv);
  return power_series;
}


inline Eigen::VectorXd _trimZeroes(const Eigen::VectorXd& vec)
{
  auto idx = vec.size();
  while (idx && std::abs(vec(idx - 1)) < _epsilon)
    --idx;
  return vec.head(idx);
}

Eigen::VectorXd _multiplyPolynomials(const Eigen::VectorXd &poly1, const Eigen::VectorXd &poly2)
{
    Eigen::VectorXd result = Eigen::VectorXd::Zero(poly1.size() + poly2.size() - 1);
    for (int i = 0; i < poly1.size(); i++)
        for (int j = 0; j < poly2.size(); j++)
            result[i + j] += poly1[i] * poly2[j];

    return result;
}

///// Curve::Curve

Curve::Curve(Eigen::MatrixX2d points, int p)
    : N_(points.rows())
    , p_(p)
    , T_(N_+p_+1)
{
    weighted_control_points_ = Eigen::MatrixX3d(N_, 3);
    weighted_control_points_.leftCols<2>() = std::move(points);
    weighted_control_points_.rightCols<1>() = Eigen::VectorXd::Ones(N_);

    // stara metoda
    uint m = N_+p_+1;
    for (uint i=0; i<p_+1; i++) {
        T_[i] = 0;
    }
    uint knots = m-2*(p_+1);
    double interval = 1.0 / (N_-p_);
    for (uint i=0; i<knots; i++) {
        T_[p_+1+i] = (i+1) * interval;
    }
    for (uint i=m-(p_+1); i<m; i++) {
        T_[i] = 1;
    }

    //spans
    for (uint i=0; i<N_-p_; i++) {
        spans.emplace_back(new Span(weighted_control_points_.middleRows(i, p_+1),
                           T_.segment(i+1, 2*p_),
                           T_(i+p_), T_(i+p_+1), p_));
    }

}

Curve::Curve(const PointVector& points)
    : N_(points.size())
    , p_(2)
    , T_(N_+p_+1)
{
    uint m = N_+p_+1;
    for (unsigned k = 0; k < N_; k++)
        weighted_control_points_.row(k).head(2) = points[k];
    weighted_control_points_.rightCols<1>() = Eigen::VectorXd::Ones(N_);

    for (uint i=0; i<p_+1; i++) {
        T_[i] = 0;
    }
    uint knots = m-2*(p_+1);
    double interval = 1.0 / (N_ - 2);
    for (uint i=0; i<knots; i++) {
        T_[p_+1+i] = (i+1) * interval;
    }
    for (uint i=m-(p_+1); i<m; i++) {
        T_[i] = 1;
    }

    //spans
    for (uint i=0; i<N_-p_; i++) {
        spans.emplace_back(new Span(weighted_control_points_.middleRows(i, p_+1),
                           T_.segment(i+1, 2*p_),
                           T_(i+p_), T_(i+p_+1), p_));
    }
}

Curve::Curve(Eigen::MatrixX3d wpoints, Eigen::ArrayXd knotvector, int p)
    : N_(wpoints.rows())
    , p_(p)
    , T_(N_+p_+1)
    , weighted_control_points_(wpoints)
{
    T_ << knotvector;
    for (int i=0; i<N_+p_+1 - knotvector.rows(); i++)
        T_ << T_.tail<1>();

    normalizeKnotVector();

    //spans
    for (uint i=0; i<N_-p_; i++) {
        spans.emplace_back(new Span(weighted_control_points_.middleRows(i, p_+1),
                           T_.segment(i+1, 2*p_),
                           T_(i+p_), T_(i+p_+1), p_));
    }

}

Curve::Curve(const Curve& curve)
    : Curve(curve.weighted_control_points_, curve.T_, curve.p_) {}

Curve& Curve::operator=(const Curve& curve)
{
  weighted_control_points_ = curve.weighted_control_points_;
  T_ = curve.T_;
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
    points[k] = controlPoint(k);
  return points;
}

Point Curve::controlPoint(unsigned idx) const
{
    Eigen::Vector3d pt = weighted_control_points_.row(idx);
    return pt.head(2) / pt(2);
}

void Curve::setControlPoint(unsigned idx, const Point& point)
{
  weighted_control_points_.row(idx).head(2) = point*weighted_control_points_(idx, 2);

  for (int i=std::max<int>(0, idx-p_); i<=idx && i<spans.size(); i++) {
      spans[i]->updateControlPoints();
  }

  resetCache();
}

std::pair<Point, Point> Curve::endPoints() const
{
    return {controlPoint(0), controlPoint(N_ - 1)};
}

void Curve::reverse() {
    weighted_control_points_ = weighted_control_points_.colwise().reverse().eval();
    resetCache();
}

PointVector Curve::polyline(double flatness) const
{
    if (!cached_polyline_ || cached_polyline_flatness_ != flatness)
    {
        cached_polyline_flatness_ = flatness;
        cached_polyline_ = std::make_unique<PointVector>();
        for(double t = T_(p_); t < T_(N_) + 0.005; t+=0.01) {
            cached_polyline_->emplace_back(valueAt(t));
        }
    }
    return *cached_polyline_;
}


Point Curve::valueAt(double t) const
{
    if (N_ == 0)
        return {0, 0};

    Span *sp = getKnotSpan(t);
    double u = (t - sp->start_t)/(sp->end_t - sp->start_t);

    Eigen::RowVectorXd pw = _powSeries(u, p_);
    return (pw * sp->cached_v_bf) / pw.dot(sp->cached_w_bf);
}


Eigen::MatrixX2d Curve::valueAt(const std::vector<double>& t_vector) const
{
    Eigen::MatrixXd out(t_vector.size(), 2);
    for (unsigned k = 0; k < t_vector.size(); k++)
        out.row(k) = valueAt(t_vector[k]);
    return out;
}


BoundingBox Curve::boundingBox() const
{
    if (!cached_bounding_box_)
    {
      auto extremes = valueAt(extrema());
      extremes.conservativeResize(extremes.rows() + 2, Eigen::NoChange);
      extremes.row(extremes.rows() - 1) = controlPoint(0);
      extremes.row(extremes.rows() - 2) = controlPoint(N_ - 1);

      cached_bounding_box_ = std::make_unique<BoundingBox>(Point(extremes.col(0).minCoeff(), extremes.col(1).minCoeff()),
                                                           Point(extremes.col(0).maxCoeff(), extremes.col(1).maxCoeff()));
    }
    return *cached_bounding_box_;
}

const Curve& Curve::derivative() const
{

}


const Curve& Curve::derivative(unsigned n) const
{

}

Vector Curve::derivativeAt(unsigned n, double t) const
{
  if (N_ == 0)
    return {0, 0};

  Span *sp = getKnotSpan(t);
  double u = (t - sp->start_t)/(sp->end_t - sp->start_t);

  // temporary solution

  Eigen::RowVectorXd pw   = _powSeries(u, p_);
  Eigen::RowVectorXd pwd1 = _powSeriesDerivative(u, p_, 1);
  Eigen::RowVectorXd pwd2 = _powSeriesDerivative(u, p_, 2);
  Eigen::RowVectorXd pwd3 = _powSeriesDerivative(u, p_, 3);

  Eigen::MatrixX2d r = sp->cached_v_bf;
  Eigen::VectorXd s = sp->cached_w_bf;

  Eigen::RowVector2d ru = pw * r;
  double su = pw.dot(s);

  // derivatives of 1/S(u) in point t
  double d1su = - pwd1.dot(s)/pow(su, 2);

  double d2su = - pwd2.dot(s) / pow(su, 2)
                + 2*pow(pwd1.dot(s), 2)/pow(su, 3);

  double d3su = - (pwd3.dot(s) / pow(su, 2))
                + 4*(pwd1.dot(s)*pwd2.dot(s)/pow(su, 3))
                - 6*(pow(pwd1.dot(s), 3)/pow(su, 4));

  switch (n) {
  case 1:
      return ru * d1su
             + (pwd1 * r) / su;
  case 2:
      return ru * d2su
             + 2*(pwd1 * r)*d1su
             + (pwd2 * r)/su;
  case 3:
      return ru * d3su
             + 3*(pwd1 * r)*d2su
             + 3*(pwd2 * r)*d1su
             + (pwd3 * r)/su;
  default:
      return valueAt(t);
  }
}

Vector Curve::derivativeAt(double t) const
{
    return derivativeAt(1, t);
}

std::vector<double> Curve::roots() const
{
  if (!cached_roots_)
  {
    cached_roots_ = std::make_unique<std::vector<double>>();
    if (N_ > 1)
    {
        Eigen::PolynomialSolver<double, Eigen::Dynamic> poly_solver;
        for (int i=0; i<spans.size(); i++)
        {
            Eigen::MatrixXd bezier_polynomial = spans[i]->cached_v_bf;

            auto trimmed_x = _trimZeroes(bezier_polynomial.col(0));
            auto trimmed_y = _trimZeroes(bezier_polynomial.col(1));

            _PolynomialRoots roots(trimmed_x.size() + trimmed_y.size());
            if (trimmed_x.size() > 1)
            {
                poly_solver.compute(trimmed_x);
                poly_solver.realRoots(roots);
            }
            if (trimmed_y.size() > 1)
            {
                poly_solver.compute(trimmed_y);
                poly_solver.realRoots(roots);
            }
            for (int j=0; j<trimmed_x.size() + trimmed_y.size(); j++)
                cached_roots_->emplace_back(roots[i]);
        }
    }
  }
  return *cached_roots_;
}

std::vector<double> Curve::extrema() const
{
    std::vector<double> extr;
    if (N_ > 1)
    {
        Eigen::PolynomialSolver<double, Eigen::Dynamic> poly_solver;
        for (int i=0; i<spans.size(); i++)
        {
            Span *sp = spans[i];

            // d/du R(u)
            Eigen::MatrixX2d p1 = Eigen::MatrixXd::Zero(p_+1, 2);
            p1.topRows(p_) = (sp->cached_v_bf.array().colwise() * _powSeriesDerivative(1, p_, 1).transpose().array()).bottomRows(p_);

            // d/du S(u)
            Eigen::RowVectorXd pb = Eigen::VectorXd::Zero(p_+1);
            pb.head(p_) = (sp->cached_w_bf.array() * _powSeriesDerivative(1, p_, 1).array()).tail(p_);

            Eigen::MatrixX2d poly(2*p_+1, 2);
            poly.col(0) = - _multiplyPolynomials(pb, sp->cached_v_bf.col(0)) + _multiplyPolynomials(p1.col(0), sp->cached_w_bf);
            poly.col(1) = - _multiplyPolynomials(pb, sp->cached_v_bf.col(1)) + _multiplyPolynomials(p1.col(1), sp->cached_w_bf);

            auto trimmed_x = _trimZeroes(poly.col(0));
            auto trimmed_y = _trimZeroes(poly.col(1));

            _PolynomialRoots roots(trimmed_x.size() + trimmed_y.size());
            if (trimmed_x.size() > 1)
            {
                poly_solver.compute(trimmed_x);
                poly_solver.realRoots(roots);
            }
            if (trimmed_y.size() > 1)
            {
                poly_solver.compute(trimmed_y);
                poly_solver.realRoots(roots);
            }
            for (int j=0; j<trimmed_x.size() + trimmed_y.size(); j++)
                if (roots[j]>=0.0 && roots[j]<=1.0)
                    extr.emplace_back(roots[j]*(sp->end_t - sp->start_t) + sp->start_t);
        }
    }
    return extr;
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

Vector Curve::tangentAt(double t, bool normalize) const
{
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

double Curve::projectPoint(const Point &point) const
{
    std::pair<double, double> min_point(0.0, (point - valueAt(0.0)).norm());

    for (int i=0; i<spans.size(); i++) {
        Span *sp = spans[i];

        if (sp->start_t == sp->end_t)
            continue;

        Eigen::MatrixX2d
                p1 = Eigen::MatrixXd::Zero(p_+1, 2);

        p1.topRows(p_) = (sp->cached_v_bf.array().colwise() * _powSeriesDerivative(1, p_, 1).transpose().array()).bottomRows(p_);

        Eigen::RowVectorXd pb = Eigen::VectorXd::Zero(p_+1);
        pb.head(p_) = (sp->cached_w_bf.array() * _powSeriesDerivative(1, p_, 1).array()).tail(p_);

        Eigen::MatrixX2d left = sp->cached_v_bf - (point*sp->cached_w_bf).transpose();

        Eigen::MatrixX2d right(2*p_+1, 2);
        right.col(0) = - _multiplyPolynomials(pb, sp->cached_v_bf.col(0)) + _multiplyPolynomials(sp->cached_w_bf, p1.col(0));
        right.col(1) = - _multiplyPolynomials(pb, sp->cached_v_bf.col(1)) + _multiplyPolynomials(sp->cached_w_bf, p1.col(1));

        Eigen::VectorXd poly = _multiplyPolynomials(left.col(0), right.col(0)) + _multiplyPolynomials(left.col(1), right.col(1));

        std::vector<double> candidates;
        Eigen::PolynomialSolver<double, Eigen::Dynamic> poly_solver;
        poly_solver.compute(_trimZeroes(poly));
        poly_solver.realRoots(candidates);
        candidates.emplace_back(1.0);

        for (int i=0; i<candidates.size(); i++) {
            double t = candidates[i] * (sp->end_t - sp->start_t) + sp->start_t;
            if (t>=0 && t<=1) {
                double dist = (point - valueAt(t)).norm();
                min_point = dist < min_point.second ? std::make_pair(t, dist) : min_point;
            }
        }
    }


    return min_point.first;
}

//todo: fix
PointVector Curve::intersections(const Curve& curve) const
{
  PointVector intersections;
  auto addIntersection = [&intersections](Point new_point) {
    // check if not already found, and add new point
    if (std::none_of(intersections.begin(), intersections.end(),
                     [&new_point](const Point& point) { return (point - new_point).norm() < _epsilon; }))
      intersections.emplace_back(std::move(new_point));
  };

  std::vector<std::pair<Curve, Curve>> subcurve_pairs;

  if (this != &curve)
    subcurve_pairs.emplace_back(*this, *this);
  else
  {
    // for self intersections divide curve into subcurves at extrema
    auto t = extrema();
    std::sort(t.begin(), t.end());
    std::vector<Curve> subcurves;
    subcurves.emplace_back(*this);
    for (unsigned k = 0; k < t.size(); k++)
    {
      Curve new_curve = std::move(subcurves.back());
      subcurves.pop_back();
      subcurves.emplace_back(new_curve.splitCurve(t[k] - _epsilon / 2).first);
      subcurves.emplace_back(new_curve.splitCurve(t[k] - _epsilon / 2).second);

#if __cpp_init_captures
      std::for_each(t.begin() + k + 1,
                    t.end(),
                    [t = t[k]](double& x) {
                        x = (x - t) / (1 - t);
      });
#else
      std::for_each(t.begin() + k + 1, t.end(), [&t, k](double& x) { x = (x - t[k]) / (1 - t[k]); });
#endif
    }

    // create all pairs of subcurves
    for (unsigned k = 0; k < subcurves.size(); k++)
      for (unsigned i = k + 1; i < subcurves.size(); i++)
        subcurve_pairs.emplace_back(subcurves[k], subcurves[i]);
  }

  while (!subcurve_pairs.empty())
  {
#if __cpp_structured_bindings
    auto [cp_a, cp_b] = std::move(subcurve_pairs.back());
#else
    Eigen::MatrixX2d cp_a, cp_b;
    std::tie(cp_a, cp_b) = std::move(subcurve_pairs.back());
#endif
    subcurve_pairs.pop_back();

    BoundingBox bbox1(cp_a.boundingBox());
    BoundingBox bbox2(cp_b.boundingBox());

    if (!bbox1.intersects(bbox2))
      ; // no intersection
    else if (bbox1.diagonal().norm() < _epsilon)
      addIntersection(bbox1.center());
    else if (bbox2.diagonal().norm() < _epsilon)
      addIntersection(bbox2.center());
    else
    {
      // intersection exists, but segments are still too large
      // - divide both segments in half
      // - insert all combinations for next iteration
      // - last pair is one where both subcurves have smallest t ranges
      auto subcurve_a = cp_a.splitCurve(0.5);
      auto subcurve_b = cp_b.splitCurve(0.5);
      subcurve_pairs.emplace_back(subcurve_a.first, subcurve_b.first);
      subcurve_pairs.emplace_back(subcurve_a.second, std::move(subcurve_b.first));
      subcurve_pairs.emplace_back(std::move(subcurve_a.first), subcurve_b.second);
      subcurve_pairs.emplace_back(std::move(subcurve_a.second), std::move(subcurve_b.second));
    }
  }

  return intersections;
}


void Curve::resetCache()
{
  N_ = weighted_control_points_.rows();
  cached_derivative_.reset();
  cached_roots_.reset();
  cached_bounding_box_.reset();
  cached_polyline_.reset();
}

Eigen::ArrayXd Curve::knotVector() const
{
    return T_;
}

void Curve::setKnot(int idx, double value) {
    if (idx > 0)
        value = std::max(value, T_(idx-1));
    else
        value = std::max(value, 0.0);
    if (idx < T_.rows()-1)
        value = std::min(value, T_(idx+1));
    else
        value = std::min(value, 1.0);
    T_(idx) = value;
    resetCache();

    for (int i=0; i<spans.size(); i++) {
        spans[i]->update();
    }
}

double Curve::knot(int idx) const
{
    return T_(idx);
}

Eigen::VectorXd Curve::weights() const
{
    return weighted_control_points_.col(2);
}

double Curve::weight(int idx) const
{
    return weighted_control_points_(idx, 2);
}

void Curve::setWeight(double w, unsigned idx)
{
    weighted_control_points_.row(idx) *= w / weighted_control_points_(idx, 2);

    for (int i=std::max<int>(0, idx-p_); i<=idx && i<spans.size(); i++) {
        spans[i]->updateControlPoints();
    }

    resetCache();
}

int Curve::getKnotSpanIndex(double t) const
{
    for (int span=p_; span<N_; span++) {
        if (T_(span+1) > t)
            return span;
    }
    return N_-1;
}

int Curve::getKnotMultiplicity(double t) const {
    return (T_ == t).count();
}

void Curve::insertKnot(double t, int r)
{
    int s = getKnotMultiplicity(t);
    int k = getKnotSpanIndex(t);
    Span* sp = spans[k-p_];
    r = std::min(r, int(p_+1-s));

    int mp = T_.rows();
    int nq = N_ + r;

    // create new knot vector
    Eigen::ArrayXd T_new(mp+r);
    T_new.head(k+1) = T_.head(k+1);
    T_new.segment(k+1, r) = Eigen::ArrayXd::Ones(r) * t;
    T_new.tail(mp-k-1) = T_.tail(mp-k-1);

    // save unaltered control points
    Eigen::MatrixX3d wpoints_new(nq, 3);
    wpoints_new.topRows(k-p_+1) = weighted_control_points_.topRows(k-p_+1);
    wpoints_new.bottomRows(N_-k+s) = weighted_control_points_.bottomRows(N_-k+s);

    //setup new control points
    Eigen::MatrixX3d wpoints_segment(sp->wpoints.topRows(p_-s+1));
    Eigen::VectorXd t_segment(sp->knots);

    for (int j=1; j<=r && j+s<=p_; j++) /* Insert the knot r times */
    {
        for (int i=0; i<=(int)p_-j-s; i++)
        {
            double alpha = (t-t_segment(i+j-1))/(t_segment(i+p_)-t_segment(i+j-1));
            wpoints_segment.row(i) = alpha*wpoints_segment.row(i+1) + (1.0-alpha)*wpoints_segment.row(i);
        }
        wpoints_new.row(k-(int)p_+j) = wpoints_segment.row(0);
        wpoints_new.row(k+r-j-s) = wpoints_segment.row((int)p_-j-s);
    }

    // load remaining control points
    for (int i=k-p_+r+1; i<k-s; i++)
        wpoints_new.row(i) = wpoints_segment.row(i-k+p_-r);

    T_ = T_new;
    weighted_control_points_ = wpoints_new;

    for (int i=0; i<r; i++) {
        spans.emplace_back(new Span(weighted_control_points_.middleRows(nq-(p_+1), p_+1),
                           T_.segment(nq-p_, 2*p_),
                           T_(nq-1), T_(nq), p_));
    }

    // reassign spans
    for (int i=0; i<spans.size(); i++) {
        new (&(spans[i]->wpoints)) Eigen::Ref<Eigen::MatrixX3d> {weighted_control_points_.middleRows(i, p_+1)};
        new (&(spans[i]->knots)) Eigen::Ref<Eigen::VectorXd> {T_.segment(i+1, 2*p_)};
        spans[i]->update();
    }

    resetCache();

}

void Curve::appendPoint(Point point)
{
    weighted_control_points_.conservativeResize(N_+1, 3);
    weighted_control_points_.row(N_).head(2) = point;
    weighted_control_points_(N_, 2) = 1.0;

    T_.conservativeResize(N_+p_+2);
    T_.head(N_+p_+1) = T_.head(N_+p_+1)*(N_-p_)/(N_-p_+1);
    T_.tail(p_+1) = 1;

    for (int i=0; i<spans.size(); i++) {
        new (&(spans[i]->wpoints)) Eigen::Ref<Eigen::MatrixX3d> {weighted_control_points_.middleRows(i, p_+1)};
        new (&(spans[i]->knots)) Eigen::Ref<Eigen::VectorXd> {T_.segment(i+1, 2*p_)};
        spans[i]->update();
    }

    spans.emplace_back(new Span(weighted_control_points_.middleRows(N_-p_, p_+1),
                       T_.segment(N_-p_+1, 2*p_),
                       T_(N_), T_(N_+1), p_));

    resetCache();
}

std::pair<Curve, Curve> Curve::splitCurve(double t) const
{
    Curve split(*this);

    split.insertKnot(t, p_+1);
    int k = split.getKnotSpanIndex(t) - (p_+1);

    int n1 = k+1, n2 = split.N_ - n1;
    Curve c1(split.weighted_control_points_.topRows(n1),
             split.T_.head(n1+p_+1),
             p_);
    Curve c2(split.weighted_control_points_.bottomRows(n2),
             split.T_.tail(n2+p_+1),
             p_);
    return {c1, c2};
}

std::vector<Curve> Curve::piecewiseBezier() const
{
    Curve temp(*this);
    std::vector<Curve> out;
    while (temp.N_ > p_+1) {
        auto pair = temp.splitCurve(temp.T_(p_+1));
        out.emplace_back(pair.first);
        temp = pair.second;
    }
    out.emplace_back(temp);
    return out;
}

void Curve::normalizeKnotVector()
{
    T_ -= T_(0);
    T_ /= T_(T_.rows()-1);
    for (int i=0; i<spans.size(); i++) {
        spans[i]->update();
    }
}

Span *Curve::getKnotSpan(double t) const
{
    for (int i=0; i<spans.size(); i++) {
        if (spans[i]->contains(t)) return spans[i];
    }
    return spans.back();
}

Eigen::VectorXd Curve::getBasisFunctionsAt(double t) const
{
    Span *sp = getKnotSpan(t);
    double u = (t - sp->start_t)/(sp->end_t - sp->start_t);
    return _powSeries(u, p_) * sp->basis_function_;
}

double Curve::length(double t) const
{
  if (t < 0.0 || t > 1.0)
    throw std::logic_error{"Length can only be calculated for t within [0.0, 1.0] range."};

  auto evaluate_chebyshev = [](double t, const Eigen::VectorXd& coeff) {
    t = 2 * t - 1;
    double tn{t}, tn_1{1}, res{coeff(0) + coeff(1) * t};
    for (unsigned k = 2; k < coeff.size(); k++)
    {
      std::swap(tn_1, tn);
      tn = 2 * t * tn_1 - tn;
      res += coeff(k) * tn;
    }
    return res;
  };

  if (!cached_chebyshev_coeffs_)
  {
    constexpr unsigned START_LOG_N = 10;
    unsigned log_n = START_LOG_N - 1;
    unsigned n = _exp2(START_LOG_N - 1);

    Eigen::VectorXd derivative_cache(2 * n + 1);
    auto updateDerivativeCache = [this, &derivative_cache](double n) {
      derivative_cache.conservativeResize(n + 1);
      derivative_cache.tail(n / 2) =
          ((1 + Eigen::cos(Eigen::ArrayXd::LinSpaced(n / 2, 1, n - 1) * M_PI / n)) / 2).unaryExpr([this](double t) {
            return derivativeAt(t).norm();
          });
    };

    derivative_cache.head(2) << derivativeAt(1.0).norm(), derivativeAt(0.0).norm();
    for (unsigned k = 2; k <= n; k *= 2)
      updateDerivativeCache(k);

    Eigen::VectorXd chebyshev;
    Eigen::FFT<double> fft;
    Eigen::VectorXcd fft_out;
    do
    {
      n *= 2;
      log_n++;
      updateDerivativeCache(n);

      unsigned N = 2 * n;
      Eigen::VectorXd coeff(N);
      coeff(0) = derivative_cache(0);
      coeff(n) = derivative_cache(1);

      for (unsigned k = 1; k <= log_n; k++)
      {
        auto lin_spaced = Eigen::ArrayXi::LinSpaced(_exp2(k - 1), 0, _exp2(k - 1) - 1);
        auto index_c = _exp2(log_n + 1 - (k + 1)) + lin_spaced * _exp2(log_n + 1 - k);
        auto index_dc = _exp2(k - 1) + 1 + lin_spaced;
        // TODO: make use of slicing & indexing in Eigen3.4
        // coeff(index_c) = coeff(N - index_c) = derivative_cache(index_dc) / n;
        for (unsigned i = 0; i < lin_spaced.size(); i++)
          coeff(index_c(i)) = coeff(N - index_c(i)) = derivative_cache(index_dc(i)) / n;
      }

      fft.fwd(fft_out, coeff);
      chebyshev = (fft_out.real().head(n - 1) - fft_out.real().segment(2, n - 1)).array() /
                  Eigen::ArrayXd::LinSpaced(n - 1, 4, 4 * (n - 1));
    } while (std::fabs(chebyshev.tail<1>()[0]) > _epsilon * 1e-2);

    unsigned cut = 0;
    while (std::fabs(chebyshev(cut)) > _epsilon * 1e-2)
      cut++;
    cached_chebyshev_coeffs_ = std::make_unique<Eigen::VectorXd>(cut + 1);
    (*cached_chebyshev_coeffs_) << 0, chebyshev.head(cut);
    (*cached_chebyshev_coeffs_)(0) = -evaluate_chebyshev(0, *cached_chebyshev_coeffs_);
  }
  return evaluate_chebyshev(t, *cached_chebyshev_coeffs_);
}

double Curve::length() const {
  return length(1.0);
}

