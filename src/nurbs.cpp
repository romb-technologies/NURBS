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

Eigen::VectorXd _multiplyPolynomials(const Eigen::VectorXd& first, const Eigen::VectorXd& second)
{
    Eigen::VectorXd out(first.rows() + second.rows() - 1);

    Eigen::MatrixXd mul = first * second.transpose();
    mul.colwise().reverseInPlace();

    for (int i=-mul.rows()+1; i<mul.cols(); i++)
        out(i+mul.rows()-1) = mul.diagonal(i).sum();
    return out;
}

///// Curve::Span

Curve::Span::Span(Eigen::Ref<Eigen::MatrixX3d> wpoints,
                  Eigen::Ref<Eigen::ArrayXd> knot_v, double start, double end, uint p) :
    wpoints(wpoints), knots(knot_v), start_t(start), end_t(end), p_(p)
{
    update();
}

bool Curve::Span::contains(double t) const
{
    if (start_t==end_t)
        return false;
    else
        return (t >= start_t) && (t <= end_t);
}

Eigen::MatrixXd Curve::Span::getBasisFunction() const
{
    return basis_function_;
}

void Curve::Span::update()
{
    // update start and end
    start_t = knots(p_-1);
    end_t = knots(p_);

    // generate basis function
    Eigen::MatrixXd m(1, 1); m<<1;

    static const int i = p_-1;
    for (int k=2; k<=p_+1; k++)
    {
        Eigen::MatrixXd m1(k, k-1),
                m2 = Eigen::MatrixXd::Zero(k-1, k),
                m3(k, k-1),
                m4 = Eigen::MatrixXd::Zero(k-1, k);

        m1 << m, Eigen::MatrixXd::Zero(1, k-1);
        m3 << Eigen::MatrixXd::Zero(1, k-1), m;

        Eigen::ArrayXd
                d0 = Eigen::ArrayXd::Constant(k-1, start_t),
                d1 = Eigen::ArrayXd::Constant(k-1, end_t - start_t),
                ddwn = Eigen::ArrayXd::Zero(k-1);

        ddwn = knots.segment(i+1, k-1) - knots.segment(i-k+2, k-1);
        d0 -= knots.segment(i-k+2, k-1);

        d0 /= ddwn;
        d1 /= ddwn;

        m2.diagonal() = 1 - d0;
        m2.diagonal(1) = d0;

        m4.diagonal() = -d1;
        m4.diagonal(1) = d1;

        m = (m1*m2)+(m3*m4);
    }
    basis_function_ = m;

    updateControlPoints();
}

void Curve::Span::updateControlPoints()
{
    // generate w_bf, v_bf
    Eigen::MatrixX3d test = wpoints;
    cached_v_bf = basis_function_ * wpoints.leftCols<2>();
    cached_w_bf = basis_function_ * wpoints.col(2);
}


///// Curve::Curve

Curve::Curve(Eigen::MatrixX2d points)
    : N_(points.rows())
    , p_(3)
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

Curve::Curve(const Curve& curve)
    : Curve(curve.weighted_control_points_.leftCols<2>()) {}

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

            Eigen::RowVectorXd temp(p_+1);
            temp << 1, 0, 0, 0;

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

double Curve::knot(int idx)
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
    int span;
    for (span=p_; span<N_; span++) {
        if (t < T_(span+1))
            return span-p_;
    }
    return N_-1;
}

void Curve::insertKnot(double t, int s, int r)
{
    int mp = T_.rows();
    int nq = N_ + r;
    int k = getKnotSpanIndex(t);

    // load new knot vector
    Eigen::ArrayXd UQ(mp+s*r);
    for (int i=0; i<=k+p_; i++) UQ(i)=T_(i);
    for (int i=1; i<=r; i++) UQ(k+p_+i)=t;
    for (int i=k+p_+1; i<mp; i++) UQ(i+r) = T_(i);

    Eigen::VectorXd test = UQ;

    // save unaltered control points
    Eigen::MatrixX3d PQ(nq, 3), Rw(p_+1, 3);
    for (uint i=0; i<=k; i++) PQ.row(i) = weighted_control_points_.row(i);
    for (uint i=k+p_-s; i<N_; i++) PQ.row(i+r) = weighted_control_points_.row(i);
    for (uint i=0; i<=p_-s; i++) Rw.row(i) = weighted_control_points_.row(k+i);

    int L = 1;
    for (int j=1; j<=r; j++) /* Insert the knot r times */
    {
        L = k+j;
        for (uint i=0; i<=p_-j-s; i++)
        {
            double alpha = (t-T_(L+i))/(T_(i+k+p_+1)-T_(L+i));
            Rw.row(i) = alpha*Rw.row(i+1) + (1.0-alpha)*Rw.row(i);
        }
        PQ.row(L) = Rw.row(0);
        PQ.row(k+p_+r-j-s) = Rw.row(p_-j-s);
    }

    // load remaining control points
    for(int i=L+1; i<k+p_-s; i++)
        PQ.row(i) = Rw.row(i-L);

    T_ = UQ;
    weighted_control_points_ = PQ;

    spans.emplace_back(new Span(weighted_control_points_.middleRows(N_-1+r-p_, p_+1),
                       T_.segment(N_+r-p_, 2*p_),
                       T_(N_+r-1), T_(N_+r), p_));

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

Curve::Span *Curve::getKnotSpan(double t) const
{
    for (int i=0; i<spans.size(); i++) {
        if (spans[i]->contains(t)) return spans[i];
    }
    return spans.back();
}
