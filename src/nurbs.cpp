#include "NURBS/nurbs.h"

using namespace NURBS;

///// Curve::Curve

Curve::Curve(Eigen::MatrixX2d points, int p) : N_(points.rows()), p_(p), T_(N_ + p_ + 1)
{
  weighted_control_points_ = Eigen::MatrixX3d(N_, 3);
  weighted_control_points_.leftCols<2>() = std::move(points);
  weighted_control_points_.rightCols<1>() = Eigen::VectorXd::Ones(N_);

  // stara metoda
  uint m = N_ + p_ + 1;
  for (uint i = 0; i < p_ + 1; i++)
  {
    T_(i) = 0;
  }
  uint knots = m - 2 * (p_ + 1);
  double interval = 1.0 / (N_ - p_);
  for (uint i = 0; i < knots; i++)
  {
    T_(p_ + 1 + i) = (i + 1) * interval;
  }
  for (uint i = m - (p_ + 1); i < m; i++)
  {
    T_(i) = 1;
  }

  // spans
  for (uint i = 0; i < N_ - p_; i++)
  {
    spans_.emplace_back(Span(weighted_control_points_.middleRows(i, p_ + 1), T_.segment(i + 1, 2 * p_), p_));
  }
}

Curve::Curve(const PointVector& points, int p) : N_(points.size()), p_(p), T_(N_ + p_ + 1)
{
  uint m = N_ + p_ + 1;

  weighted_control_points_ = Eigen::MatrixX3d(N_, 3);
  for (unsigned k = 0; k < N_; k++)
    weighted_control_points_.row(k).head(2) = points[k];
  weighted_control_points_.rightCols<1>() = Eigen::VectorXd::Ones(N_);

  for (uint i = 0; i < p_ + 1; i++)
  {
    T_(i) = 0;
  }
  uint knots = m - 2 * (p_ + 1);
  double interval = 1.0 / (N_ - 2);
  for (uint i = 0; i < knots; i++)
  {
    T_(p_ + 1 + i) = (i + 1) * interval;
  }
  for (uint i = m - (p_ + 1); i < m; i++)
  {
    T_(i) = 1;
  }

  // spans
  for (uint i = 0; i < N_ - p_; i++)
  {
    spans_.emplace_back(Span(weighted_control_points_.middleRows(i, p_ + 1), T_.segment(i + 1, 2 * p_), p_));
  }
}

Curve::Curve(Eigen::MatrixX3d wpoints, Eigen::ArrayXd knotvector, int p)
    : N_(wpoints.rows()), p_(p), T_(N_ + p_ + 1), weighted_control_points_(wpoints)
{
  T_ << knotvector;
  for (int i = 0; i < N_ + p_ + 1 - knotvector.rows(); i++)
    T_ << T_.tail<1>();

  normalizeKnotVector();

  // spans
  for (uint i = 0; i < N_ - p_; i++)
  {
    spans_.emplace_back(Span(weighted_control_points_.middleRows(i, p_ + 1), T_.segment(i + 1, 2 * p_), p_));
  }
}

Curve::Curve(const Curve& curve) : Curve(curve.weighted_control_points_, curve.T_, curve.p_) {}

Curve& Curve::operator=(const Curve& curve)
{
  weighted_control_points_ = curve.weighted_control_points_;
  T_ = curve.T_;
  p_ = curve.p_;
  resetCache();
  return *this;
}

unsigned Curve::order() const { return p_; }

int Bin(int n, int k)
{
  if (k == 0 || k == n)
    return 1;
  return Bin(n - 1, k - 1) + Bin(n - 1, k);
}

// (wake me up) wake me up inside (i can't wake up) wake me up inside (save me)
void Curve::elevateOrder(uint t)
{

  uint new_p = p_ + t;

  Eigen::MatrixXd bezalfs(new_p + 1, p_ + 1);
  Eigen::VectorXd alfs(p_ - 1);
  Eigen::MatrixX3d bpts(p_ + 1, 3), ebpts(p_ + t + 1, 3), Nextbpts(p_ + 1, 3);

  Eigen::MatrixX3d new_wpoints(N_ + 2 * t, 3);
  Eigen::ArrayXd new_t(2 * N_ + new_p + 1);

  /* Compute Bezier degree elevation coefficients */
  bezalfs(0, 0) = bezalfs(new_p, p_) = 1.0;
  for (int i = 1; i <= new_p / 2; i++)
  {
    double inv = 1.0 / Bin(new_p, i);
    int mpi = std::min((int)p_, i);
    for (int j = std::max(i - (int)t, 0); j <= mpi; j++)
      bezalfs(i, j) = inv * Bin(p_, j) * Bin(t, i - j);
  }
  for (int i = new_p / 2 + 1; i <= new_p - 1; i++)
  {
    int mpi = std::min((int)p_, i);
    for (int j = std::max(i - (int)t, 0); j <= mpi; j++)
      bezalfs(i, j) = bezalfs(new_p - i, p_ - j);
  }
  int new_m = new_p + 1, kind = new_p + 1, a = p_, r = -1, b = p_ + 1, cind = 1;
  double ua = T_(0), ub = T_(N_ + p_);
  new_wpoints.row(0) = weighted_control_points_.row(0);
  for (int i = 0; i <= new_p; i++)
    new_t(i) = ua;
  for (int i = 0; i <= p_; i++)
    bpts.row(i) = weighted_control_points_.row(i);

  while (b < N_ + p_ + 1)
  {
    int mul = getKnotMultiplicity(T_(b));
    b += mul - 1;
    new_m += mul + t;
    new_t.conservativeResize(new_m);
    ub = T_(b);
    int oldr = r;
    r = p_ - mul;
    int lbz, rbz;

    if (oldr > 0)
      lbz = (oldr + 2) / 2;
    else
      lbz = 1;

    if (r > 0)
      rbz = new_p - (r + 1) / 2;
    else
      rbz = new_p;
    if (r > 0)
    { /* Insert knot to get Bezier segment */
      double numer = ub - ua;
      for (int k = p_; k > mul; k--)
        alfs(k - mul - 1) = numer / (T_(a + k) - ua);
      for (int j = 1; j <= r; j++)
      {
        int save = r - j;
        int s = mul + j;
        for (int k = p_; k >= s; k--)
        {
          bpts.row(k) = alfs(k - s) * bpts.row(k) + (1.0 - alfs(k - s)) * bpts.row(k - 1);
        }
        Nextbpts.row(save) = bpts.row(p_);
      }
    } /* End of "insert knot" */

    for (int i = lbz; i <= new_p; i++)
    /* Degree elevate Bezier */
    { /* Only points lbz, ... ,ph are used below */
      ebpts.row(i) = Eigen::Vector3d(0.0, 0.0, 0.0);
      int mpi = std::min(int(p_), i);
      for (int j = std::max(0, i - (int)t); j <= mpi; j++)
        ebpts.row(i) = ebpts.row(i) + bezalfs(i, j) * bpts.row(j);
    } /* End of degree elevating Bezier */

    if (oldr > 1)
    { /* Must remove knot u=U[a] oldr times */
      int first = kind - 2;
      int last = kind;
      double den = ub - ua;
      double bet = (ub - new_t(kind - 1)) / den;
      for (int tr = 1; tr < oldr; tr++)
      { /* Knot removal loop */
        int i = first;
        int j = last;
        int kj = j - kind + 1;
        while (j - i > tr) /* Loop and compute the new */
        {                  /* control points for one removal step */
          if (i < cind)
          {
            double alf = (ub - new_t(i)) / (ua - new_t(i));
            new_wpoints.row(i) = alf * new_wpoints.row(i) + (1.0 - alf) * new_wpoints.row(i - 1);
          }
          if (j >= lbz)
          {
            if (j - tr <= kind - new_p + oldr)
            {
              double gam = (ub - new_t(j - tr)) / den;
              ebpts.row(kj) = gam * ebpts.row(kj) + (1.0 - gam) * ebpts.row(kj + 1);
            }
            else
            {
              ebpts.row(kj) = bet * ebpts.row(kj) + (1.0 - bet) * ebpts.row(kj + 1);
            }
          }
          i++;
          j--;
          kj--;
        }
        first--;
        last++;
      }
    } /* End of knot removal */

    if (a != p_)
      for (int i = 0; i < new_p - oldr; i++)
      {
        new_t(kind) = ua;
        kind = kind + 1;
      }

    new_t(kind) = ua;
    for (int j = lbz; j <= rbz; j++)
    {
      new_wpoints.row(cind) = ebpts.row(j);
      new_wpoints(cind, 2) = 1.0;
      cind++;
    }
    if (b < N_ + p_ + 1) // set up for next pass through loop
    {
      for (int j = 0; j < r; j++)
        bpts.row(j) = Nextbpts.row(j);
      for (int j = r; j <= p_; j++)
        bpts.row(j) = weighted_control_points_.row(b - p_ + j);
      a = b++;
      ua = ub;
    }
  }
  for (int i = 0; i <= new_p; i++)
    new_t(kind + i) = ub;

  N_ = new_m - new_p - 1;
  weighted_control_points_ = new_wpoints.topRows(N_);
  p_ = new_p;
  T_ = new_t.head(new_m);

  spans_.clear();
  for (uint i = 0; i < N_ - p_; i++)
  {
    spans_.emplace_back(Span(weighted_control_points_.middleRows(i, p_ + 1), T_.segment(i + 1, 2 * p_), p_));
  }

  resetCache();
}

void Curve::lowerOrder()
{
  uint new_p = p_ - 1, new_m = new_p, kind = new_p + 1;
  int r = -1, a = p_, b = p_ + 1, cind = 1, m = N_ + p_ + 1;

  Eigen::MatrixXd bezalfs(new_p + 1, p_ + 1);
  Eigen::VectorXd alphas(p_ - 1);
  Eigen::MatrixX3d bpts(p_ + 1, 3), rbpts(p_, 3), // (p - num. of orders + 1)
      Nextbpts(p_ + 1, 3);

  Eigen::MatrixX3d new_wpoints(N_, 3);
  Eigen::ArrayXd new_t(N_ + p_);

  new_wpoints.row(0) = weighted_control_points_.row(0);
  new_t.head(new_p + 1) = T_(0);
  bpts.topRows(p_ + 1) = weighted_control_points_.topRows(p_ + 1);

  int lbz, rbz;
  while (b < m)
  {
    int mul = getKnotMultiplicity(T_(b));
    while (b < m - 1 && std::fabs(T_(b) - T_(b + 1)) < _epsilon)
      b = b + 1;
    new_m += mul + 1;
    int oldr = r;
    r = p_ - mul;
    if (oldr > 0)
      lbz = (oldr + 2) / 2;
    else
      lbz = 1;
    if (r > 0)
    {
      double numer = T_(b) - T_(a);
      for (int k = p_; k > mul; k--)
        alphas(k - mul - 1) = numer / (T_(a + k) - T_(a));
      for (int j = 1; j <= r; j++)
      {
        int save = r - j;
        int s = mul + j;
        for (int k = p_; k >= s; k--)
          bpts.row(k) = alphas(k - s) * bpts.row(k) + (1.0 - alphas(k - s)) * bpts.row(k - 1);
        Nextbpts.row(save) = bpts.row(p_);
      }
    }
    /* Degree reduce Bezier segment */
    // rbpts = BezDegreeReduce(bpts);
    Eigen::MatrixXd elevate_order_coeffs = Eigen::MatrixXd::Zero(p_ + 1, p_);
    elevate_order_coeffs.diagonal().setLinSpaced(1, 1 - (p_ + 1) / p_);
    elevate_order_coeffs.diagonal(-1).setLinSpaced(1. / p_, 1);

    Eigen::MatrixXd lower_order_coeffs =
        (elevate_order_coeffs.transpose() * elevate_order_coeffs).inverse() * elevate_order_coeffs.transpose();

    rbpts = lower_order_coeffs * bpts;

    if (oldr > 0)
    {
      int first = kind;
      int last = kind;
      int i, j;
      for (int k = 0; k < oldr; k++)
      {
        i = first;
        j = last;
        int kj = j - kind;
        while (j - i > k)
        {
          double alfa = (T_(a) - new_t(i - 1)) / (T_(b) - new_t(i - 1));
          double beta = (T_(a) - new_t(j - k - 1)) / (T_(b) - new_t(j - k - 1));
          new_wpoints.row(i - 1) -= (1.0 - alfa) * new_wpoints.row(i - 2) / alfa;
          rbpts.row(kj) = (rbpts.row(kj) - beta * rbpts.row(kj + 1)) / (1.0 - beta);
          i++;
          j--;
          kj--;
        }
        first--;
        last++;
      } /* End for (k=O; k<oldr; k++) loop */
      cind = i - 1;
    } /* End if (oldr > 0) */
    /* Load knot vector and control points */
    if (a != p_)
      new_t.segment(kind, new_p - oldr) = T_(a);
    for (int i = lbz; i <= new_p; i++)
    {
      new_wpoints.row(cind++) = rbpts.row(i);
    }
    /* Set up for next pass through */
    if (b < m)
    {
      int i;
      for (i = 0; i < r; i++)
        bpts.row(i) = Nextbpts.row(i);
      bpts.row(i) = weighted_control_points_.row(b - p_ + i);
      for (int i = r; i <= p_; i++)
        b++;
      a = b;
    }
  } /* End of while (b < m) loop */
  for (int i = 0; i <= new_p; i++)
    new_t(kind + i) = T_(b);

  weighted_control_points_ = new_wpoints;
  N_ = new_m - new_p - 1;
  p_ = new_p;
  T_ = new_t.head(new_m);

  spans_.clear();
  for (uint i = 0; i < N_ - p_; i++)
  {
    spans_.emplace_back(Span(weighted_control_points_.middleRows(i, p_ + 1), T_.segment(i + 1, 2 * p_), p_));
  }

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
  weighted_control_points_.row(idx).head(2) = point * weighted_control_points_(idx, 2);

  for (int i = std::max<int>(0, idx - p_); i <= idx && i < spans_.size(); i++)
  {
    spans_[i].updateControlPoints(weighted_control_points_.middleRows(i, p_ + 1));
  }

  resetCache();
}

std::pair<Point, Point> Curve::endPoints() const { return {valueAt(T_(p_ + 1)), valueAt(T_(N_))}; }

void Curve::reverse()
{
  weighted_control_points_ = weighted_control_points_.colwise().reverse().eval();
  resetCache();
}

PointVector Curve::polyline(double flatness) const
{
  if (!cached_polyline_)
  {
    cached_polyline_ = PointVector();
    for (int i = 0; i < spans_.size(); i++)
    {
      PointVector poly = spans_[i].polyline();
      cached_polyline_->insert(cached_polyline_->end(), poly.begin(), poly.end());
    }
  }
  return cached_polyline_.value();
}

Point Curve::valueAt(double t) const
{
  if (N_ == 0)
    return {0, 0};
  const Span& sp = getKnotSpan(t);
  double u = (t - sp.start()) / (sp.end() - sp.start());
  return sp.valueAt(u);
}

Eigen::MatrixX2d Curve::valueAt(const std::vector<double>& t_vector) const
{
  Eigen::MatrixXd out(t_vector.size(), 2);
  for (unsigned k = 0; k < t_vector.size(); k++)
    out.row(k) = valueAt(t_vector[k]);
  return out;
}

BoundingBox Curve::boundingBox(bool use_roots) const
{
  if (use_roots)
  {
    if (!cached_bounding_box_)
    {
      auto extremes = valueAt(extrema());
      extremes.conservativeResize(extremes.rows() + 2, Eigen::NoChange);
      extremes.row(extremes.rows() - 1) = endPoints().first;
      extremes.row(extremes.rows() - 2) = endPoints().second;

      cached_bounding_box_ = BoundingBox(Point(extremes.col(0).minCoeff(), extremes.col(1).minCoeff()),
                                         Point(extremes.col(0).maxCoeff(), extremes.col(1).maxCoeff()));
    }
    return *cached_bounding_box_;
  }
  else
  {
    Eigen::MatrixX2d unweighted_pts = weighted_control_points_.leftCols(2).array().colwise() / weights().array();
    return BoundingBox(Point(unweighted_pts.col(0).minCoeff(), unweighted_pts.col(1).minCoeff()),
                       Point(unweighted_pts.col(0).maxCoeff(), unweighted_pts.col(1).maxCoeff()));
  }
}

const Curve& Curve::derivative() const {}

const Curve& Curve::derivative(unsigned n) const {}

Vector Curve::derivativeAt(unsigned n, double t) const
{
  if (N_ == 0)
    return {0, 0};

  const Span& sp = getKnotSpan(t);
  double u = (t - sp.start()) / (sp.end() - sp.start());

  return sp.derivativeAt(n, u);
}

Vector Curve::derivativeAt(double t) const { return derivativeAt(1, t); }

std::vector<double> Curve::roots() const
{
  if (!cached_roots_)
  {
    cached_roots_ = std::vector<double>();
    if (N_ > 1)
    {
      Eigen::PolynomialSolver<double, Eigen::Dynamic> poly_solver;
      for (const auto& span : spans_)
      {
        Eigen::MatrixXd bezier_polynomial = span.cachedVBF();

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
        for (int j = 0; j < trimmed_x.size() + trimmed_y.size(); j++)
          cached_roots_->emplace_back(roots[j]);
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
    for (const auto& span : spans_)
    {
      std::vector<double> extr_sp = span.extrema();
      for (auto ex : extr_sp)
      {
        extr.emplace_back(ex * (span.end() - span.start()) + span.start());
      }
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

double Curve::projectPoint(const Point& point) const
{
  std::pair<double, double> min_point(0.0, (point - valueAt(0.0)).norm());

  for (int i = 0; i < spans_.size(); i++)
  {
    const Span& sp = spans_[i];

    if (std::fabs(sp.start() - sp.end()) < _epsilon)
      continue;

    Eigen::MatrixX2d p1 = Eigen::MatrixXd::Zero(p_ + 1, 2);

    p1.topRows(p_) =
        (sp.cachedVBF().array().colwise() * _powSeriesDerivative(1, p_, 1).transpose().array()).bottomRows(p_);

    Eigen::RowVectorXd pb = Eigen::VectorXd::Zero(p_ + 1);
    pb.head(p_) = (sp.cachedWBF().array() * _powSeriesDerivative(1, p_, 1).array()).tail(p_);

    Eigen::MatrixX2d left = sp.cachedVBF() - (point * sp.cachedWBF()).transpose();

    Eigen::MatrixX2d right(2 * p_ + 1, 2);
    right.col(0) = -_multiplyPolynomials(pb, sp.cachedVBF().col(0)) + _multiplyPolynomials(sp.cachedWBF(), p1.col(0));
    right.col(1) = -_multiplyPolynomials(pb, sp.cachedVBF().col(1)) + _multiplyPolynomials(sp.cachedWBF(), p1.col(1));

    Eigen::VectorXd poly =
        _multiplyPolynomials(left.col(0), right.col(0)) + _multiplyPolynomials(left.col(1), right.col(1));

    std::vector<double> candidates;
    Eigen::PolynomialSolver<double, Eigen::Dynamic> poly_solver;
    auto trim = _trimZeroes(poly);
    if (trim.size() > 0)
      poly_solver.compute(trim);
    poly_solver.realRoots(candidates);
    candidates.emplace_back(1.0);

    for (int i = 0; i < candidates.size(); i++)
    {
      double t = candidates[i] * (sp.end() - sp.start()) + sp.start();
      if (t >= 0 && t <= 1)
      {
        double dist = (point - valueAt(t)).norm();
        min_point = dist < min_point.second ? std::make_pair(t, dist) : min_point;
      }
    }
  }

  return min_point.first;
}

PointVector Curve::intersections(const Curve& other) const
{
  PointVector intersections;
  for (int i = 0; i < spans_.size(); i++)
    for (int j = 0; j < other.spans_.size(); j++)
    {
      if (this == &other && j < i)
        continue;
      PointVector ints = spans_[i].intersections(other.spans_[j]);
      intersections.insert(intersections.end(), std::make_move_iterator(ints.begin()),
                           std::make_move_iterator(ints.end()));
    }
  return intersections;
}

void Curve::resetCache()
{
  N_ = weighted_control_points_.rows();
  cached_roots_.reset();
  cached_bounding_box_.reset();
  cached_polyline_.reset();
  for (int i = 0; i < spans_.size(); i++)
  {
    spans_[i].resetCache();
  }
}

Eigen::ArrayXd Curve::knotVector() const { return T_; }

void Curve::setKnot(int idx, double value)
{
  if (idx < 0 || idx >= T_.rows())
    return;

  if (idx > 0)
    value = std::max(value, T_(idx - 1));
  else
    value = std::max(value, 0.0);
  if (idx < T_.rows() - 1)
    value = std::min(value, T_(idx + 1));
  else
    value = std::min(value, 1.0);
  T_(idx) = value;
  resetCache();

  // todo: update only affected knots
  for (int i = 0; i < spans_.size(); i++)
  {
    spans_[i].update(T_.segment(i + 1, 2 * p_), weighted_control_points_.middleRows(i, p_ + 1));
  }
}

double Curve::knot(int idx) const { return T_(idx); }

Eigen::VectorXd Curve::weights() const { return weighted_control_points_.col(2); }

double Curve::weight(int idx) const { return weighted_control_points_(idx, 2); }

void Curve::setWeight(int idx, double value)
{
  weighted_control_points_.row(idx) *= value / weighted_control_points_(idx, 2);

  for (int i = std::max<int>(0, idx - p_); i <= idx && i < spans_.size(); i++)
  {
    spans_[i].updateControlPoints(weighted_control_points_.middleRows(i, p_ + 1));
  }

  resetCache();
}

int Curve::getKnotSpanIndex(double t) const
{
  for (int span = p_; span < N_; span++)
  {
    if (T_(span + 1) > t)
      return span;
  }
  return N_ - 1;
}

int Curve::getKnotMultiplicity(double t) const { return ((T_ - t).abs() < _epsilon).count(); }

void Curve::insertKnot(double t, int r)
{
  int s = getKnotMultiplicity(t);
  int k = getKnotSpanIndex(t);
  const Span& sp = spans_[k - p_];
  r = std::min(r, int(p_ + 1 - s));

  int mp = T_.rows();
  int nq = N_ + r;

  // create new knot vector
  Eigen::ArrayXd T_new(mp + r);
  T_new.head(k + 1) = T_.head(k + 1);
  T_new.segment(k + 1, r) = Eigen::ArrayXd::Ones(r) * t;
  T_new.tail(mp - k - 1) = T_.tail(mp - k - 1);

  // save unaltered control points
  Eigen::MatrixX3d wpoints_new(nq, 3);
  wpoints_new.topRows(k - p_ + 1) = weighted_control_points_.topRows(k - p_ + 1);
  wpoints_new.bottomRows(N_ - k + s) = weighted_control_points_.bottomRows(N_ - k + s);

  // setup new control points
  Eigen::MatrixX3d wpoints_segment(weighted_control_points_.middleRows(k - p_, p_ - s + 1));
  Eigen::VectorXd t_segment(T_.segment(k - p_ + 1, 2 * p_));

  for (int j = 1; j <= r && j + s <= p_; j++) /* Insert the knot r times */
  {
    for (int i = 0; i <= (int)p_ - j - s; i++)
    {
      double alpha = (t - t_segment(i + j - 1)) / (t_segment(i + p_) - t_segment(i + j - 1));
      wpoints_segment.row(i) = alpha * wpoints_segment.row(i + 1) + (1.0 - alpha) * wpoints_segment.row(i);
    }
    wpoints_new.row(k - (int)p_ + j) = wpoints_segment.row(0);
    wpoints_new.row(k + r - j - s) = wpoints_segment.row((int)p_ - j - s);
  }

  // load remaining control points
  for (int i = k - p_ + r + 1; i < k - s; i++)
    wpoints_new.row(i) = wpoints_segment.row(i - k + p_ - r);

  T_ = T_new;
  weighted_control_points_ = wpoints_new;

  for (int i = 0; i < r; i++)
  {
    spans_.emplace_back(
        Span(weighted_control_points_.middleRows(nq - (p_ + 1), p_ + 1), T_.segment(nq - p_, 2 * p_), p_));
  }

  // reassign spans
  for (int i = 0; i < spans_.size(); i++)
  {
    spans_[i].update(T_.segment(i + 1, 2 * p_), weighted_control_points_.middleRows(i, p_ + 1));
  }

  resetCache();
}

void Curve::appendPoint(Point point)
{
  weighted_control_points_.conservativeResize(N_ + 1, 3);
  weighted_control_points_.row(N_).head(2) = point;
  weighted_control_points_(N_, 2) = 1.0;

  T_.conservativeResize(N_ + p_ + 2);
  T_.head(N_ + p_ + 1) = T_.head(N_ + p_ + 1) * (N_ - p_) / (N_ - p_ + 1);
  T_.tail(p_ + 1) = 1;

  for (int i = 0; i < spans_.size(); i++)
  {
    spans_[i].update(T_.segment(i + 1, 2 * p_), weighted_control_points_.middleRows(i, p_ + 1));
  }

  spans_.emplace_back(Span(weighted_control_points_.middleRows(N_ - p_, p_ + 1), T_.segment(N_ - p_ + 1, 2 * p_), p_));

  resetCache();
}

std::pair<Curve, Curve> Curve::splitCurve(double t) const
{
  Curve split(*this);

  split.insertKnot(t, p_ + 1);
  int k = split.getKnotSpanIndex(t) - (p_ + 1);

  int n1 = k + 1, n2 = split.N_ - n1;
  Curve c1(split.weighted_control_points_.topRows(n1), split.T_.head(n1 + p_ + 1), p_);
  Curve c2(split.weighted_control_points_.bottomRows(n2), split.T_.tail(n2 + p_ + 1), p_);
  return {c1, c2};
}

std::vector<Curve> Curve::piecewiseBezier() const
{
  Curve temp(*this);
  std::vector<Curve> out;
  while (temp.N_ > p_ + 1)
  {
    auto pair = temp.splitCurve(temp.T_(p_ + 1));
    out.emplace_back(pair.first);
    temp = pair.second;
  }
  out.emplace_back(temp);
  return out;
}

void Curve::normalizeKnotVector()
{
  T_ -= T_(0);
  T_ /= T_(T_.rows() - 1);
  for (int i = 0; i < spans_.size(); i++)
  {
    spans_[i].update(T_.segment(i + 1, 2 * p_), weighted_control_points_.middleRows(i, p_ + 1));
  }
}

Span& Curve::getKnotSpan(double t) const
{
  for (int i = 0; i < spans_.size(); i++)
  {
    if (spans_[i].contains(t))
      return spans_[i];
  }
  return spans_.back();
}

Eigen::VectorXd Curve::getBasisFunctionsAt(double t) const
{
  const Span& sp = getKnotSpan(t);
  double u = (t - sp.start()) / (sp.end() - sp.start());
  return _powSeries(u, p_) * sp.basisFunction();
}

double Curve::length(double t) const
{
  // analytic
  if (t < _epsilon)
    return 0.0;

  int ix = 0;
  double len = 0.0;

  for (Span& sp : spans_)
  {
    len += sp.length();
  }

  return len;
}

double Curve::length() const { return length(1.0); }

void Curve::removeKnot(int ix, int k)
{
  while (T_(ix) >= T_(ix + 1))
    ix++;
  int s = getKnotMultiplicity(T_(ix));

  int first = ix - p_ + 1;
  int last = ix - s - 1;
  Eigen::MatrixX3d Pt;

  int t, i, j;
  Pt = Eigen::MatrixX3d(T_.rows(), 3);
  for (t = 0; t < k && t < s; t++)
  {

    first--;
    last++;
    int off = first - 1;
    Pt.row(0) = weighted_control_points_.row(off);
    Pt.row(last + 1 - off) = weighted_control_points_.row(last + 1);

    i = first;
    j = last;
    while (j - i > t)
    {
      double alpha_i = (T_(ix) - T_(i)) / (T_(i + p_ + 1 + t) - T_(i));
      double alpha_j = (T_(ix) - T_(j - t)) / (T_(j + p_ + 1) - T_(j - t));

      Pt.row(i - off) = (weighted_control_points_.row(i) - (1 - alpha_i) * Pt.row(i - off - 1)) / alpha_i;
      Pt.row(j - off) = (weighted_control_points_.row(j) - alpha_j * Pt.row(j - off + 1)) / (1 - alpha_j);

      i++;
      j--;
    }

    i = first;
    j = last;
    while (j - i > t)
    {
      weighted_control_points_.row(i) = Pt.row(i - off);
      weighted_control_points_.row(j) = Pt.row(j - off);
      i++;
      j--;
    }
  }
  i = (2 * ix - s - p_) / 2;
  j = i;
  for (int m = 1; m < t; m++)
  {
    if (m % 2 == 1)
      i++;
    else
      j--;
  }

  for (int kn = ix + 1; kn < N_ + p_ + 1; kn++)
  {
    T_(kn - t) = T_(kn);
  }
  for (int m = i + 1; m < N_; m++)
  {
    weighted_control_points_.row(j++) = weighted_control_points_.row(m);
  }

  weighted_control_points_.conservativeResize(N_ - t, 3);
  T_.conservativeResize(N_ + p_ + 1 - t);

  for (int m = 0; m < t; m++)
  {

    spans_.pop_back();
  }
  // reassign spans
  for (int i = 0; i < spans_.size(); i++)
  {
    spans_[i].update(T_.segment(i + 1, 2 * p_), weighted_control_points_.middleRows(i, p_ + 1));
  }

  resetCache();
}

Curve Curve::join(Curve& other)
{
  // todo: elevate order
  if (p_ != other.p_)
    return *this;

  auto ends1 = endPoints();
  auto ends2 = other.endPoints();

  if ((ends1.second - ends2.second).norm() < (ends1.second - ends2.first).norm())
    other.reverse();
  if ((ends1.first - ends2.first).norm() < (ends1.second - ends2.first).norm())
    this->reverse();
  if ((ends1.first - ends2.second).norm() < (ends1.first - ends2.first).norm())
  {
    this->reverse();
    other.reverse();
  }

  Eigen::MatrixX3d points(N_ + other.N_, 3);
  points << weighted_control_points_, other.weighted_control_points_;

  Eigen::ArrayXd knots(points.rows() + p_ + 1);
  knots << T_.head(N_ + p_ - 1), other.T_.tail(other.N_ + p_ - 1) + T_(N_ + p_);

  return Curve(points.leftCols(2), p_);
}

void Curve::applyContinuity(const Curve& source_curve, const std::vector<double>& beta_coeffs)
{
  unsigned c_order = std::min((uint)beta_coeffs.size(), p_);

  // pascal triangle matrix (binomial coefficients) - rowwise
  Eigen::MatrixXd pascal_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  pascal_matrix.row(0).setOnes();
  for (unsigned k = 1; k <= c_order; k++)
    for (unsigned i = 1; i <= k; i++)
      pascal_matrix(i, k) = pascal_matrix(i - 1, k - 1) + pascal_matrix(i, k - 1);

  // inverse of pascal matrix, i.e., pascal matrix with alternating signs - colwise
  Eigen::MatrixXd pascal_alternating_matrix = pascal_matrix.transpose().inverse();

  // https://en.wikipedia.org/wiki/Bell_polynomials -> equivalent to equations of geometric continuity
  Eigen::MatrixXd bell_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  bell_matrix(0, c_order) = 1;
  for (unsigned k = 0; k < c_order; k++)
    bell_matrix.block(1, c_order - k - 1, k + 1, 1) =
        bell_matrix.block(0, c_order - k, k + 1, k + 1) *
        pascal_matrix.block(0, k, k + 1, 1)
            .cwiseProduct(Eigen::Map<const Eigen::MatrixXd>(beta_coeffs.data(), k + 1, 1));

  // diagonal: (N-1)! / (N-k-1)!
  Eigen::MatrixXd factorial_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  factorial_matrix(0, 0) = 1;
  for (unsigned k = 1; k <= c_order; k++)
    factorial_matrix(k, k) = factorial_matrix(k - 1, k - 1) * (N_ - k);

  // derivatives of given curve
  Eigen::Matrix2Xd derivatives(Eigen::Index(2), Eigen::Index(c_order + 1));
  for (unsigned k = 0; k < c_order + 1; k++)
    derivatives.col(k) = source_curve.derivativeAt(k, 1.0);

  // based on the beta coefficients and geometric continuity equations, calculate new derivatives
  Eigen::MatrixXd new_derivatives = (derivatives * bell_matrix).rowwise().reverse().transpose();

  weighted_control_points_.topRows(c_order + 1).leftCols(2) =
      (factorial_matrix * pascal_alternating_matrix).inverse() * new_derivatives;
  resetCache();
}
