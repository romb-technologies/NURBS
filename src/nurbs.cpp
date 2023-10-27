#include "NURBS/nurbs.h"

#include <numeric>
#include <limits>

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
        spans_.emplace_back(Span(weighted_control_points_.middleRows(i, p_+1),
                           T_.segment(i+1, 2*p_),
                           T_(i+p_), T_(i+p_+1), p_));
    }

}

Curve::Curve(const PointVector& points, int p)
    : N_(points.size())
    , p_(p)
    , T_(N_+p_+1)
{
    uint m = N_+p_+1;

    weighted_control_points_ = Eigen::MatrixX3d(N_, 3);
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
        spans_.emplace_back(Span(weighted_control_points_.middleRows(i, p_+1),
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
        spans_.emplace_back(Span(weighted_control_points_.middleRows(i, p_+1),
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

int Bin(int n, int k) {
    if (k == 0 || k == n)
       return 1;
    return Bin(n - 1, k - 1) + Bin(n - 1, k);
}

// (wake me up) wake me up inside (i can't wake up) wake me up inside (save me)
void Curve::elevateOrder(uint t) {

    uint new_p = p_+t;
    
    Eigen::MatrixXd bezalfs(new_p+1, p_+1);
    Eigen::VectorXd alfs(p_-1);
    Eigen::MatrixX3d
            bpts(p_+1, 3),
            ebpts(p_+t+1, 3),
            Nextbpts(p_+1, 3);

    Eigen::MatrixX3d new_wpoints(2*N_, 3);
    Eigen::ArrayXd new_t(2*N_+new_p+1);
    
    /* Compute Bezier degree elevation coefficients */
    bezalfs(0, 0) = bezalfs(new_p, p_) = 1.0;
    for (uint i=1; i<=new_p/2; i++)
    {
        double inv = 1.0/Bin(new_p,i);
        uint mpi = std::min(p_,i);
        for (uint j=std::max(i-t, 0U); j<=mpi; j++)
            bezalfs(i, j) = inv*Bin(p_,j)*Bin(t,i-j);
    }
    for (uint i=new_p/2+1; i<=new_p-1; i++)
    {
        uint mpi = std::min(p_,i);
        for (uint j=std::max(0U, i-t); j<=mpi; j++)
            bezalfs(i, j) = bezalfs(new_p-i, p_-j);
    }
    int new_m = new_p+1, kind = new_p+1, a = p_, r = -1,
            b = p_+1, cind = 1;
    double ua = T_(0), ub = T_(N_+p_);
    new_wpoints.row(0) = weighted_control_points_.row(0);
    for (int i=0; i<=new_p; i++)
        new_t(i) = ua;
    for (int i=0; i<=p_; i++)
        bpts.row(i) = weighted_control_points_.row(i);

    while (b < N_+p_+1)
    {
        int mul = getKnotMultiplicity(T_(b));
        b += mul - 1;
        new_m += mul+t;
        double ub = T_(b);
        int oldr = r;
        r = p_-mul;
        int lbz, rbz;

        if (oldr > 0) lbz = (oldr+2)/2;
        else lbz = 1;

        if (r > 0)  rbz = new_p-(r+1)/2;
        else rbz = new_p;
        if (r > 0)
        { /* Insert knot to get Bezier segment */
            double numer = ub-ua;
            for (int k=p_; k>mul; k--)
                alfs(k-mul-1) = numer/(T_(a+k)-ua);
            for (int j=1; j<=r; j++)
            {
                int save = r-j;
                int s = mul+j;
                for (int k=p_; k>=s; k--)
                {
                    bpts.row(k) = alfs(k-s)*bpts.row(k) +
                            (1.0-alfs(k-s))*bpts.row(k-1) ;
                }
                Nextbpts.row(save) = bpts.row(p_);
            }
        } /* End of "insert knot" */

        for (uint i=lbz; i<=new_p; i++)
            /* Degree elevate Bezier */
        { /* Only points lbz, ... ,ph are used below */
            ebpts.row(i) = Eigen::Vector3d(0.0, 0.0, 0.0);
            int mpi = std::min(p_,i);
            for (int j=std::max(0U, i-t); j<=mpi; j++)
                ebpts.row(i) = ebpts.row(i) + bezalfs(i, j)*bpts.row(j);
        } /* End of degree elevating Bezier */

        if (oldr > 1)
        { /* Must remove knot u=U[a] oldr times */
            int first = kind-2;
            int last = kind;
            int den = ub-ua;
            int bet = (ub-new_t(kind-1))/den;
            for (int tr=1; tr<oldr; tr++)
            { /* Knot removal loop */
                int i = first;
                int j = last;
                int kj = j-kind+1;
                while (j-i > tr) /* Loop and compute the new */
                { /* control points for one removal step */
                    if (i < cind)
                    {
                        double alf = (ub-new_t(i))/(ua-new_t(i));
                        new_wpoints.row(i) = alf*new_wpoints.row(i)
                                + (1.0-alf)*new_wpoints.row(i-1);
                    }
                    if(j >= lbz)
                    {
                        if( j-tr <= kind-new_p+oldr )
                        {
                            double gam = (ub-new_t(j-tr))/den;
                            ebpts.row(kj) = gam*ebpts.row(kj)+(1.0-gam)*ebpts.row(kj+1);
                        }
                        else
                        {
                            ebpts.row(kj) = bet*ebpts.row(kj)+(1.0-bet)*ebpts.row(kj+1);
                        }
                    }
                    i++; j--; kj--;

                }
                first--;
                last++;
            }
        } /* End of knot removal */

        if (a != p_)
            for (int i=0; i<new_p-oldr; i++) {
                new_t(kind) = ua;
                kind = kind+1;
            }

        new_t(kind) = ua;
        for (int j=lbz; j<=rbz; j++)
        {
            new_wpoints.row(cind) = ebpts.row(j);
            new_wpoints(cind, 2) = 1.0;
            cind++;
        }
        if (b < N_+p_+1) // set up for next pass through loop
        {
            for (int j=0; j<r; j++) bpts.row(j) = Nextbpts.row(j);
            for (int j=r; j<=p_; j++) bpts.row(j) = weighted_control_points_.row(b-p_+j);
            a = b++; ua = ub;
        }
    }
    for (int i=0; i<=new_p; i++)
        new_t(kind+i) = ub;

    int new_n = new_m-new_p-1;
    Eigen::VectorXd test_t = new_t.matrix();

    weighted_control_points_ = new_wpoints.topRows(new_n);
    N_ = new_n;
    p_ = new_p;
    T_ = new_t.head(new_m);

    spans_.clear();
    for (uint i=0; i<N_-p_; i++) {
        spans_.emplace_back(Span(weighted_control_points_.middleRows(i, p_+1),
                                 T_.segment(i+1, 2*p_),
                                 T_(i+p_), T_(i+p_+1), p_));
    }

    resetCache();
}

void Curve::lowerOrder() {
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

  for (int i=std::max<int>(0, idx-p_); i<=idx && i<spans_.size(); i++) {
      spans_[i].updateControlPoints();
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
    if (!cached_polyline_)
    {
        cached_polyline_ = std::make_unique<PointVector>();
        for (int i=0; i<spans_.size(); i++) {
            PointVector poly = spans_[i].polyline();
            cached_polyline_->insert(cached_polyline_->end(), poly.begin(), poly.end());
        }
    }
    return *cached_polyline_;
}


Point Curve::valueAt(double t) const
{
    if (N_ == 0)
        return {0, 0};
    const Span &sp = getKnotSpan(t);
    double u = (t - sp.start_t_)/(sp.end_t_ - sp.start_t_);
    return sp.valueAt(u);
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

  const Span &sp = getKnotSpan(t);
  double u = (t - sp.start_t_)/(sp.end_t_ - sp.start_t_);

  return sp.derivativeAt(n, u);
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
        for (int i=0; i<spans_.size(); i++)
        {
            Eigen::MatrixXd bezier_polynomial = spans_[i].cached_vbf_;

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
        for (int i=0; i<spans_.size(); i++)
        {
            const Span &sp = spans_[i];

            // d/du R(u)
            Eigen::MatrixX2d p1 = Eigen::MatrixXd::Zero(p_+1, 2);
            p1.topRows(p_) = (sp.cached_vbf_.array().colwise() * _powSeriesDerivative(1, p_, 1).transpose().array()).bottomRows(p_);

            // d/du S(u)
            Eigen::RowVectorXd pb = Eigen::VectorXd::Zero(p_+1);
            pb.head(p_) = (sp.cached_wbf_.array() * _powSeriesDerivative(1, p_, 1).array()).tail(p_);

            Eigen::MatrixX2d poly(2*p_+1, 2);
            poly.col(0) = - _multiplyPolynomials(pb, sp.cached_vbf_.col(0)) + _multiplyPolynomials(p1.col(0), sp.cached_wbf_);
            poly.col(1) = - _multiplyPolynomials(pb, sp.cached_vbf_.col(1)) + _multiplyPolynomials(p1.col(1), sp.cached_wbf_);

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
                    extr.emplace_back(roots[j]*(sp.end_t_ - sp.start_t_) + sp.start_t_);
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

    for (int i=0; i<spans_.size(); i++) {
        const Span &sp = spans_[i];

        if (sp.start_t_ == sp.end_t_)
            continue;

        Eigen::MatrixX2d
                p1 = Eigen::MatrixXd::Zero(p_+1, 2);

        p1.topRows(p_) = (sp.cached_vbf_.array().colwise() * _powSeriesDerivative(1, p_, 1).transpose().array()).bottomRows(p_);

        Eigen::RowVectorXd pb = Eigen::VectorXd::Zero(p_+1);
        pb.head(p_) = (sp.cached_wbf_.array() * _powSeriesDerivative(1, p_, 1).array()).tail(p_);

        Eigen::MatrixX2d left = sp.cached_vbf_ - (point*sp.cached_wbf_).transpose();

        Eigen::MatrixX2d right(2*p_+1, 2);
        right.col(0) = - _multiplyPolynomials(pb, sp.cached_vbf_.col(0)) + _multiplyPolynomials(sp.cached_wbf_, p1.col(0));
        right.col(1) = - _multiplyPolynomials(pb, sp.cached_vbf_.col(1)) + _multiplyPolynomials(sp.cached_wbf_, p1.col(1));

        Eigen::VectorXd poly = _multiplyPolynomials(left.col(0), right.col(0)) + _multiplyPolynomials(left.col(1), right.col(1));

        std::vector<double> candidates;
        Eigen::PolynomialSolver<double, Eigen::Dynamic> poly_solver;
        poly_solver.compute(_trimZeroes(poly));
        poly_solver.realRoots(candidates);
        candidates.emplace_back(1.0);

        for (int i=0; i<candidates.size(); i++) {
            double t = candidates[i] * (sp.end_t_ - sp.start_t_) + sp.start_t_;
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
  for (int i=0; i<spans_.size(); i++) {
      spans_[i].resetCache();
  }
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

    for (int i=0; i<spans_.size(); i++) {
        spans_[i].update();
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

void Curve::setWeight(int idx, double value)
{
    weighted_control_points_.row(idx) *= value / weighted_control_points_(idx, 2);

    for (int i=std::max<int>(0, idx-p_); i<=idx && i<spans_.size(); i++) {
        spans_[i].updateControlPoints();
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
    const Span &sp = spans_[k-p_];
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
    Eigen::MatrixX3d wpoints_segment(sp.wpoints_.topRows(p_-s+1));
    Eigen::VectorXd t_segment(sp.knots);

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
        spans_.emplace_back(Span(weighted_control_points_.middleRows(nq-(p_+1), p_+1),
                           T_.segment(nq-p_, 2*p_),
                           T_(nq-1), T_(nq), p_));
    }

    // reassign spans
    for (int i=0; i<spans_.size(); i++) {
        new (&(spans_[i].wpoints_)) Eigen::Ref<Eigen::MatrixX3d> {weighted_control_points_.middleRows(i, p_+1)};
        new (&(spans_[i].knots)) Eigen::Ref<Eigen::VectorXd> {T_.segment(i+1, 2*p_)};
        spans_[i].update();
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

    for (int i=0; i<spans_.size(); i++) {
        new (&(spans_[i].wpoints_)) Eigen::Ref<Eigen::MatrixX3d> {weighted_control_points_.middleRows(i, p_+1)};
        new (&(spans_[i].knots)) Eigen::Ref<Eigen::VectorXd> {T_.segment(i+1, 2*p_)};
        spans_[i].update();
    }

    spans_.emplace_back(Span(weighted_control_points_.middleRows(N_-p_, p_+1),
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
    for (int i=0; i<spans_.size(); i++) {
        spans_[i].update();
    }
}

Span &Curve::getKnotSpan(double t) const
{
    for (int i=0; i<spans_.size(); i++) {
        if (spans_[i].contains(t)) return spans_[i];
    }
    return spans_.back();
}

Eigen::VectorXd Curve::getBasisFunctionsAt(double t) const
{
    const Span &sp = getKnotSpan(t);
    double u = (t - sp.start_t_)/(sp.end_t_ - sp.start_t_);
    return _powSeries(u, p_) * sp.getBasisFunction();
}

double Curve::length(double t) const
{
    // analytic
    if (t == 0.0) return 0.0;

    int ix = 0;
    double len = 0.0;

    for (Span& sp: spans_) {
        len += sp.length();
    }

    return len;
}

double Curve::length() const
{
    // polyline
    double out = 0.0;
    PointVector poly = polyline();
    for(int i=0; i<poly.size()-1; i++) {
        out += std::sqrt(pow(poly[i](0)-poly[i+1](0), 2) + pow(poly[i](1)-poly[i+1](1), 2));
    }
    return out;
    //    return length(1.0);
}

void Curve::removeKnot(int ix, int k)
{
    while(T_(ix) >= T_(ix+1)) ix++;
    int s = getKnotMultiplicity(T_(ix));

    int first = ix-p_+1;
    int last = ix-s-1;
    Eigen::MatrixX3d Pt;

    int t, i, j;
    Pt = Eigen::MatrixX3d(T_.rows(), 3);
    for (t=0; t<k && t<s; t++) {

        first--; last++;
        int off = first-1;
        Pt.row(0) = weighted_control_points_.row(off);
        Pt.row(last+1-off) = weighted_control_points_.row(last+1);

        i=first; j=last;
        while(j-i > t) {
            double alpha_i = (T_(ix)-T_(i))/(T_(i+p_+1+t)-T_(i));
            double alpha_j = (T_(ix)-T_(j-t))/(T_(j+p_+1)-T_(j-t));

            Pt.row(i-off) = (weighted_control_points_.row(i) - (1-alpha_i)*Pt.row(i-off-1)) / alpha_i;
            Pt.row(j-off) = (weighted_control_points_.row(j) - alpha_j*Pt.row(j-off+1)) / (1-alpha_j);

            i++; j--;
        }

        i=first; j=last;
        while(j-i > t) {
            weighted_control_points_.row(i) = Pt.row(i-off);
            weighted_control_points_.row(j) = Pt.row(j-off);
            i++; j--;
        }

    }
    i = (2*ix-s-p_)/2; j = i;
    for (int m=1; m<t; m++){
        if (m%2 == 1)
            i++;
        else
            j--;
    }

    for (int kn=ix+1; kn<N_+p_+1; kn++) {
        T_(kn-t) = T_(kn);
    }
    for (int m=i+1; m<N_; m++) {
        weighted_control_points_.row(j++) = weighted_control_points_.row(m);
    }

    weighted_control_points_.conservativeResize(N_-t, 3);
    T_.conservativeResize(N_+p_+1-t);

    for (int m=0; m<t; m++) {

        spans_.pop_back();

    }
    // reassign spans
    for (int i=0; i<spans_.size(); i++) {
        new (&(spans_[i].wpoints_)) Eigen::Ref<Eigen::MatrixX3d> {weighted_control_points_.middleRows(i, p_+1)};
        new (&(spans_[i].knots)) Eigen::Ref<Eigen::VectorXd> {T_.segment(i+1, 2*p_)};
        spans_[i].update();
    }

    resetCache();
}

Curve Curve::join(Curve &other)
{
    // todo: elevate order
    if (p_ != other.p_) return *this;

    auto ends1 = endPoints();
    auto ends2 = other.endPoints();

    if (dist(ends1.second, ends2.second) < dist(ends1.second, ends2.first))
        other.reverse();
    if (dist(ends1.first, ends2.first) < dist(ends1.second, ends2.first))
        this->reverse();
    if (dist(ends1.first, ends2.second) < dist(ends1.first, ends2.first)) {
        this->reverse();
        other.reverse();
    }

    Eigen::MatrixX3d points(N_ + other.N_, 3);
    points << weighted_control_points_, other.weighted_control_points_;

    Eigen::ArrayXd knots(points.rows() + p_ + 1);
    knots << T_.head(N_+p_-1), other.T_.tail(other.N_+p_-1) + T_(N_+p_);

    return Curve(points.leftCols(2), p_);
}

void Curve::applyContinuity(const Curve& source_curve, const std::vector<double>& beta_coeffs)
{
  unsigned c_order = beta_coeffs.size();

  Eigen::MatrixXd pascal_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  Eigen::MatrixXd pascal_alterating_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  pascal_alterating_matrix.diagonal(-1).setLinSpaced(-1, -static_cast<int>(c_order));
  pascal_alterating_matrix = pascal_alterating_matrix.exp();
  pascal_matrix = pascal_alterating_matrix.cwiseAbs().transpose();

  Eigen::MatrixXd bell_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  bell_matrix(0, c_order) = 1;

  for (unsigned k = 0; k < c_order; k++)
    bell_matrix.block(1, c_order - k - 1, k + 1, 1) =
        bell_matrix.block(0, c_order - k, k + 1, k + 1) *
        pascal_matrix.block(0, k, k + 1, 1)
            .cwiseProduct(Eigen::Map<const Eigen::MatrixXd>(beta_coeffs.data(), k + 1, 1));

  Eigen::MatrixXd factorial_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));

  factorial_matrix.diagonal() = Eigen::ArrayXd::LinSpaced(c_order + 1, 0, c_order).unaryExpr([this](unsigned k) {
    // (N-1)! / (N-k-1)! = e^(ln(N-1)! - ln(N-k-1)!)
    return std::exp(std::lgamma(N_) - std::lgamma(N_ - k));
  });

  Eigen::Matrix2Xd derivatives(Eigen::Index(2), Eigen::Index(c_order + 1));
  for (unsigned k = 0; k < c_order + 1; k++)
    derivatives.col(k) = source_curve.derivativeAt(k, 1.0);

  Eigen::MatrixXd derivatives_wanted = (derivatives * bell_matrix).rowwise().reverse().transpose();

  weighted_control_points_.topRows(c_order + 1).leftCols(2) = (factorial_matrix * pascal_alterating_matrix).inverse() * derivatives_wanted;
  resetCache();
}






